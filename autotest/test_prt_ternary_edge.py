"""
Particles tracked with the ternary method along internal subcell edges.

The ternary method splits each cell into triangular subcells sharing the
cell center. On a uniform grid with a CHD in opposite corners, flow runs
along the grid's diagonal (x + y = 10), so a particle released on it, e.g.
at a cell center, travels exactly along the internal edge between two
subcells, from the cell center to the cell corner.

The exit point at the corner could land just beyond the end of the cell
face due to roundoff, and the adjacent (internal) subcell face, along
which the normal velocity is zero, reports no exit either. The particle
was then terminated with status 9 (no exit from the subcell), typically
when restarting at the beginning of a later time step.

Particles are released at each interior cell center along the diagonal
and at points along the internal edges, with several time steps, so the
solution restarts at varying positions along each edge. Each particle
should travel to the CHD cell in the bottom right corner.
"""

from pathlib import Path

import flopy
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtedge"
cases = [simname]
nrow = ncol = 10
top, botm = 1.0, [0.0]
perlen, nstp = 1.0, 20
porosity = 0.1
TERM_NO_EXITS, TERM_NO_EXITS_SUB = 5, 9


def release_points():
    # cell (i, i) has center (i + 0.5, 9.5 - i), on the diagonal x + y = 10.
    # also release partway between the center and the corners, on the
    # internal subcell edges.
    pts = []
    for i in range(1, nrow - 1):
        icpl = i * ncol + i
        for f in (0.0, 0.3, -0.3):
            pts.append((icpl, i + 0.5 + f, 9.5 - i - f))
    return pts


def build_models(test):
    mf6 = test.targets["mf6"]
    disv_kwargs = get_disv_kwargs(
        nlay=1, nrow=nrow, ncol=ncol, delr=1.0, delc=1.0, tp=top, botm=botm
    )

    gwf_sim = flopy.mf6.MFSimulation(
        sim_name=test.name, exe_name=mf6, sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(gwf_sim, nper=1, perioddata=[(perlen, nstp, 1.0)])
    gwfname = get_model_name(test.name, "gwf")
    gwf = flopy.mf6.ModflowGwf(gwf_sim, modelname=gwfname, save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **disv_kwargs)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0), top], [(0, nrow * ncol - 1), 0.0]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(gwf_sim)

    prt_ws = Path(test.workspace) / "prt"
    prt_sim = flopy.mf6.MFSimulation(sim_name=test.name, exe_name=mf6, sim_ws=prt_ws)
    flopy.mf6.ModflowTdis(prt_sim, nper=1, perioddata=[(perlen, nstp, 1.0)])
    prtname = get_model_name(test.name, "prt")
    prt = flopy.mf6.ModflowPrt(prt_sim, modelname=prtname)
    flopy.mf6.ModflowPrtdisv(prt, **disv_kwargs)
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(release_points()),
        packagedata=[
            [irpt, (0, icpl), x, y, 0.5]
            for irpt, (icpl, x, y) in enumerate(release_points())
        ],
        perioddata={0: [("FIRST",)]},
        dev_forceternary=True,
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
    )
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", Path(test.workspace) / f"{gwfname}.hds"),
            ("GWFBUDGET", Path(test.workspace) / f"{gwfname}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(prt_sim, pname="ems", filename=f"{prtname}.ems")
    prt_sim.register_solution_package(ems, [prt.name])
    return gwf_sim, prt_sim


def check_output(test):
    prtname = get_model_name(test.name, "prt")
    trk = pd.read_csv(Path(test.workspace) / "prt" / f"{prtname}.trk.csv")
    assert not (trk.istatus == TERM_NO_EXITS_SUB).any(), (
        "particle(s) terminated in a subcell with no exit face:\n"
        f"{trk[trk.istatus == TERM_NO_EXITS_SUB]}"
    )
    last = trk.groupby("irpt").tail(1)
    assert len(last) == len(release_points())
    # every particle ends up in the CHD cell in the bottom right corner
    assert (last.icell == nrow * ncol).all(), last
    assert (last.istatus == TERM_NO_EXITS).all(), last


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build_models,
        check=check_output,
        targets=targets,
        compare=None,
    )
    test.run()
