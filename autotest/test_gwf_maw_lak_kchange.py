"""
Test MAW and LAK connection conductance when aquifer hydraulic
conductivity changes during a simulation.

Steady flow from a well or lake at a constant head, connected to the
first cell of a row of cells, to a constant head at the last cell. TVK
changes K at the connected cell. After the change, the heads and the
well or lake flows must match a model with that K given in NPF.

Case 1 (kc-maw): MAW, THIEM conductance equation; K changes at the start
    of period 2.
Case 2 (kc-maws): MAW, SKIN conductance equation; K changes at the start
    of period 2.
Case 3 (kc-mawts): MAW, THIEM conductance equation; K from a time series
    that changes at the second time step of period 2.
Case 4 (kc-lak): LAK, horizontal connection without a lakebed; K changes
    at the start of period 2.
Case 5 (kc-lakv): LAK above the cell, vertical connection with a lakebed
    leakance; K changes at the start of period 2.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["kc-maw", "kc-maws", "kc-mawts", "kc-lak", "kc-lakv"]
timeseries = [False, False, True, False, False]

nlay, nrow, ncol = 1, 1, 10
delr, delc = 10.0, 10.0
top, botm = 10.0, 0.0
k0, k1 = 10.0, 100.0
cell = (0, 0, 0)


def build_sim(ws, name, exe, tvk):
    idx = cases.index(name)
    gwfname = f"gwf-{name}"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 2, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    flopy.mf6.ModflowIms(
        sim, inner_dvclose=1e-10, outer_dvclose=1e-10, linear_acceleration="BICGSTAB"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    k = np.full((nlay, nrow, ncol), k0)
    if not tvk:
        # reference: the K after the change given in NPF
        k[cell] = k1
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    if tvk and timeseries[idx]:
        tvkpkg = flopy.mf6.ModflowUtltvk(
            npf, perioddata={1: [(cell, "K", "kts")]}, filename=f"{gwfname}.tvk"
        )
        tvkpkg.ts.initialize(
            filename=f"{gwfname}.tvk.ts",
            timeseries=[(0.0, k0), (1.5, k1), (2.0, k1)],
            time_series_namerecord="kts",
            interpolation_methodrecord="stepwise",
        )
    elif tvk:
        npf.tvk.initialize(filename=f"{gwfname}.tvk", perioddata={1: [(cell, "K", k1)]})
    flopy.mf6.ModflowGwfic(gwf, strt=1.5)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, ncol - 1), 1.0]])

    if name.startswith("kc-maw"):
        condeqn = "SKIN" if name == "kc-maws" else "THIEM"
        flopy.mf6.ModflowGwfmaw(
            gwf,
            save_flows=True,
            packagedata=[(0, 0.1, botm, 2.0, condeqn, 1)],
            connectiondata=[(0, 0, cell, top, botm, 1.0, 0.2)],
            perioddata={0: [(0, "status", "constant"), (0, "well_head", 2.0)]},
        )
    elif name == "kc-lak":
        flopy.mf6.ModflowGwflak(
            gwf,
            save_flows=True,
            nlakes=1,
            packagedata=[(0, 2.0, 1)],
            connectiondata=[
                (0, 0, cell, "horizontal", "none", botm, top, delr / 2, delc)
            ],
            perioddata={0: [(0, "status", "constant"), (0, "stage", 2.0)]},
        )
    else:
        flopy.mf6.ModflowGwflak(
            gwf,
            save_flows=True,
            nlakes=1,
            packagedata=[(0, 12.0, 1)],
            connectiondata=[(0, 0, cell, "vertical", 1.0, 0.0, 0.0, 0.0, 0.0)],
            perioddata={0: [(0, "status", "constant"), (0, "stage", 12.0)]},
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, cases[idx], "mf6", tvk=True), None


def check_output(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"

    # reference: the K after the change given in NPF; steady flow, so the
    # time steps after the change are compared (period 2, or its second
    # time step for the time series case)
    ws = test.workspace / "ref"
    sim = build_sim(ws, name, test.targets["mf6"], tvk=False)
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference simulation failed"
    steps = [2] if timeseries[idx] else [1, 2]

    h = flopy.utils.HeadFile(test.workspace / f"{gwfname}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{gwfname}.hds").get_alldata()
    assert np.abs(h[2] - h[0]).max() > 1e-3, "TVK change not applied"
    for i in steps:
        assert np.allclose(h[i], href[i], atol=1e-8), (
            f"time step {i + 1} heads differ from NPF reference: "
            f"max diff {np.abs(h[i] - href[i]).max()}"
        )

    # well or lake flow to the aquifer
    text = "MAW" if name.startswith("kc-maw") else "LAK"
    q = flopy.utils.CellBudgetFile(
        test.workspace / f"{gwfname}.cbc", precision="double"
    ).get_data(text=text)
    qref = flopy.utils.CellBudgetFile(
        ws / f"{gwfname}.cbc", precision="double"
    ).get_data(text=text)
    for i in steps:
        assert np.allclose(q[i]["q"], qref[i]["q"], atol=1e-8), (
            f"time step {i + 1} {text} flows differ from NPF reference: "
            f"max diff {np.abs(q[i]['q'] - qref[i]['q']).max()}"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
