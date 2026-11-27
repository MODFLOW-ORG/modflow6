"""
Tests particle tracking through a rectangular DISV cell whose
east side is divided among several smaller neighbors. Flow is
west to east. Particles released uniformly across the cell's
west face should split among the east neighbors in proportion
to the flow each neighbor receives.

The quad-refined Pollock method applies only to sides divided
once at the midpoint. Sides divided into thirds, or once off
center, are tracked with the ternary method.
"""

import itertools

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

simname = "prtqhang"
cases = {
    f"{simname}mid": [1.5],
    f"{simname}off": [1.0],
    f"{simname}thd": [1.0, 2.0],
}
nreleasepts = 60


def get_gridprops(splits):
    """
    West strip cell [0,1]x[0,3], center cell [1,4]x[0,3], and
    east cells [4,5]x[y0,y1] for each interval between splits.
    """
    ys = [0.0, *splits, 3.0]
    verts = [(0.0, 0.0), (1.0, 0.0), (0.0, 3.0), (1.0, 3.0)]
    verts += [(4.0, y) for y in ys] + [(5.0, y) for y in ys]
    iv4 = {y: 4 + i for i, y in enumerate(ys)}
    iv5 = {y: 4 + len(ys) + i for i, y in enumerate(ys)}
    cells = [
        [0, 2, 3, 1],
        [1, 3, *[iv4[y] for y in reversed(ys)]],
    ]
    for y0, y1 in itertools.pairwise(ys):
        cells.append([iv4[y0], iv4[y1], iv5[y1], iv5[y0]])
    cell2d = []
    for i, iverts in enumerate(cells):
        xc, yc = np.mean([verts[iv] for iv in iverts], axis=0)
        cell2d.append([i, xc, yc, len(iverts), *iverts])
    return {
        "nlay": 1,
        "ncpl": len(cells),
        "nvert": len(verts),
        "top": 1.0,
        "botm": [0.0],
        "vertices": [[i, x, y] for i, (x, y) in enumerate(verts)],
        "cell2d": cell2d,
    }


def build_models(idx, test):
    name = test.name
    gwf_name = f"{name}_gwf"
    prt_name = f"{name}_prt"
    gridprops = get_gridprops(cases[name])
    ncpl = gridprops["ncpl"]

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **gridprops)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)
    chd = [[(0, 0), 1.0]] + [[(0, i), 0.5] for i in range(2, ncpl)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{gwf_name}.bud", saverecord=[("BUDGET", "ALL")]
    )
    ims = flopy.mf6.ModflowIms(sim)
    sim.register_solution_package(ims, [gwf_name])

    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdisv(prt, **gridprops)
    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)
    dy = 3.0 / nreleasepts
    rpts = [[i, 0, 1, 1.001, (i + 0.5) * dy, 0.5] for i in range(nreleasepts)]
    flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=nreleasepts,
        packagedata=rpts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(prt, trackcsv_filerecord=f"{prt_name}.trk.csv")
    flopy.mf6.ModflowGwfprt(sim, exgmnamea=gwf_name, exgmnameb=prt_name)
    ems = flopy.mf6.ModflowEms(sim)
    sim.register_solution_package(ems, [prt_name])
    return sim


def check_output(idx, test):
    name = test.name
    sim = test.sims[0]
    gwf = sim.get_model(f"{name}_gwf")
    ncpl = gwf.modelgrid.ncpl
    east = range(3, ncpl + 1)

    # flow into each east cell, from the CHD budget
    chd = gwf.output.budget().get_data(text="CHD")[0]
    qeast = {int(n): -q for n, q in zip(chd["node"], chd["q"]) if n in east}
    qfrac = np.array([qeast[n] for n in east]) / sum(qeast.values())

    # each particle enters an east cell within that cell's extent
    pls = pd.read_csv(test.workspace / f"{name}_prt.trk.csv")
    entries = pls[pls.icell.isin(east)].groupby("irpt").first()
    assert len(entries) == nreleasepts
    ys = [0.0, *cases[name], 3.0]
    ylo = entries.icell.map(lambda n: ys[n - 3])
    yhi = entries.icell.map(lambda n: ys[n - 2])
    assert entries.y.between(ylo, yhi).all(), entries[["icell", "y"]]

    # particles split among east cells in proportion to flow
    pfrac = np.array([(entries.icell == n).sum() for n in east]) / nreleasepts
    assert np.allclose(pfrac, qfrac, atol=2.0 / nreleasepts), (pfrac, qfrac)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()
