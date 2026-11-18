"""
Particles released on a cell face they exit through immediately.

A particle released on the face it's moving toward exits at once, so its
travel time to that face is zero. Pollock's method computes the travel
time as log(v_exit / v) / (dv/dx), where v is interpolated from the face
velocities at the particle's position. On the face, v should equal the
face velocity, and the log of their ratio be exactly zero. But a compiler
may reassociate the interpolation (e.g. ifx's default floating point model,
or gfortran -ffast-math), so v differs from the face velocity in the last
bit, giving a tiny positive or negative travel time. A negative one moves
the particle back in time, and could reorder or drop events reported at
the release time (e.g. a tracking time).

Particles are released at the midpoints of the faces of interior cells in
a heterogeneous flow field (lognormal K), so face velocities vary. There
is a tracking time at the release time. Checks:

    - no event is reported at a negative time,
    - every particle reports the tracking time at the release time,
    - an exit at (or within roundoff of) the release time is at exactly 0.

With standard IEEE arithmetic (e.g. gfortran without -ffast-math) the
interpolation is exact on the face, so this passes with or without the
fix. It guards builds with value-unsafe floating point optimizations.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtzdt"
nrow = ncol = 10
top, botm = 1.0, [0.0]
porosity = 0.1
GRIDS = ["dis", "disv"]
cases = [f"{simname}_{grid}" for grid in GRIDS]
RELEASE, EXIT, USERTIME = 0, 1, 5


def hk():
    rng = np.random.default_rng(42)
    return np.exp(rng.normal(0.0, 1.0, (1, nrow, ncol)))


def add_dis(model, grid):
    gwf = isinstance(model, flopy.mf6.ModflowGwf)
    if grid == "dis":
        cls = flopy.mf6.ModflowGwfdis if gwf else flopy.mf6.ModflowPrtdis
        cls(
            model,
            nlay=1,
            nrow=nrow,
            ncol=ncol,
            delr=1.0,
            delc=1.0,
            top=top,
            botm=botm,
        )
    else:
        cls = flopy.mf6.ModflowGwfdisv if gwf else flopy.mf6.ModflowPrtdisv
        cls(
            model,
            **get_disv_kwargs(
                nlay=1,
                nrow=nrow,
                ncol=ncol,
                delr=1.0,
                delc=1.0,
                tp=top,
                botm=botm,
            ),
        )


def cellid(grid, i, j):
    return (0, i, j) if grid == "dis" else (0, i * ncol + j)


def release_points(grid):
    """Midpoints of the 4 faces of each interior cell."""
    pts = []
    for i in range(1, nrow - 1):
        for j in range(1, ncol - 1):
            xc, yc = j + 0.5, nrow - i - 0.5
            for dx, dy in ((0.5, 0.0), (-0.5, 0.0), (0.0, 0.5), (0.0, -0.5)):
                pts.append([len(pts), cellid(grid, i, j), xc + dx, yc + dy, 0.5])
    return pts


def build_models(test, grid):
    mf6 = test.targets["mf6"]
    ws = Path(test.workspace)

    gwf_sim = flopy.mf6.MFSimulation(sim_name=test.name, exe_name=mf6, sim_ws=ws)
    flopy.mf6.ModflowTdis(gwf_sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwfname = get_model_name(test.name, "gwf")
    gwf = flopy.mf6.ModflowGwf(gwf_sim, modelname=gwfname, save_flows=True)
    add_dis(gwf, grid)
    k = hk() if grid == "dis" else hk().reshape(1, -1)
    flopy.mf6.ModflowGwfnpf(
        gwf, k=k, save_specific_discharge=True, save_saturation=True
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [cellid(grid, 0, 0), top],
            [cellid(grid, nrow - 1, ncol - 1), 0.0],
        ],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(gwf_sim)

    prt_ws = ws / "prt"
    prt_sim = flopy.mf6.MFSimulation(sim_name=test.name, exe_name=mf6, sim_ws=prt_ws)
    flopy.mf6.ModflowTdis(prt_sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    prtname = get_model_name(test.name, "prt")
    prt = flopy.mf6.ModflowPrt(prt_sim, modelname=prtname)
    add_dis(prt, grid)
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)
    pts = release_points(grid)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(pts),
        packagedata=pts,
        nreleasetimes=1,
        releasetimes=[(0.0,)],
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
        ntracktimes=1,
        tracktimes=[(0.0,)],
    )
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", ws / f"{gwfname}.hds"),
            ("GWFBUDGET", ws / f"{gwfname}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(prt_sim, pname="ems", filename=f"{prtname}.ems")
    prt_sim.register_solution_package(ems, [prt.name])
    return gwf_sim, prt_sim


def check_output(test, grid):
    prtname = get_model_name(test.name, "prt")
    trk = pd.read_csv(Path(test.workspace) / "prt" / f"{prtname}.trk.csv")
    assert trk.irpt.nunique() == len(release_points(grid))

    neg = trk[trk.t < 0]
    assert neg.empty, f"events at negative times:\n{neg}"

    lost, tiny, immediate = [], [], 0
    for irpt, p in trk.groupby("irpt"):
        if not ((p.ireason == USERTIME) & (p.t == 0.0)).any():
            lost.append(irpt)
        exits = p[p.ireason == EXIT]
        if exits.empty:
            continue
        t = exits.t.iloc[0]
        if t == 0.0:
            immediate += 1
        elif abs(t) < 1e-9:
            tiny.append((irpt, t))
    assert not lost, f"tracking time at release not reported: {lost}"
    assert not tiny, f"exits within roundoff of, but not at, release: {tiny}"
    # sanity check: about half the particles start on an outflow face
    assert immediate > len(release_points(grid)) // 4


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    grid = name.split("_")[-1]
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, grid),
        check=lambda t: check_output(t, grid),
        targets=targets,
        compare=None,
    )
    test.run()
