"""
Test that EVTA with READASARRAYS and AUXMULTNAME correctly persists the
auxiliary multiplier array across stress periods when it is not re-read.

The multiplier is specified only in period 1. The ET rate is updated in
period 2 without re-specifying the multiplier. The test verifies that the
multiplier from period 1 continues to be applied.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["evta_auxmult01"]

# spatial discretization
nlay, nrow, ncol = 1, 5, 5
delr = delc = 1000.0
top = 10.0
botm = [0.0]
strt = 9.0
k = 10.0

# ET parameters
et_surface = 10.0
et_rate = 0.005
et_depth = 5.0
mult_val = 0.5


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nper = 3

    # ET rate changes in period 2; multiplier only in period 1
    rate1 = np.full((nrow, ncol), et_rate)
    rate2 = np.full((nrow, ncol), et_rate * 2.0)
    evt_rate = {0: rate1, 1: rate2}
    evt_aux = {0: [np.full((nrow, ncol), mult_val)]}

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-9, inner_dvclose=1e-9)
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, k=k)

    # constant head on one side
    chd_spd = [[(0, i, ncol - 1), strt] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

    # evta with auxmultname
    flopy.mf6.ModflowGwfevta(
        gwf,
        readasarrays=True,
        surface=et_surface,
        rate=evt_rate,
        depth=et_depth,
        auxiliary="ET_MULT",
        auxmultname="ET_MULT",
        aux=evt_aux,
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    fpth = test.workspace / f"{name}.cbc"
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    evt_budgets = cobj.get_data(text="EVT")

    # Key checks:
    # 1. ET is non-zero in all periods (multiplier persists)
    # 2. Period 2 ET > period 1 ET (rate doubled, mult same)
    # 3. Period 3 ET == period 2 ET (both persist)
    for iper, evt_data in enumerate(evt_budgets):
        q = evt_data["q"]
        et_out = -q.min()
        assert et_out > 0.0, (
            f"Period {iper + 1}: ET is zero - AUXMULTNAME not persisting across periods"
        )

    et1 = -evt_budgets[0]["q"].sum()
    et2 = -evt_budgets[1]["q"].sum()
    et3 = -evt_budgets[2]["q"].sum()

    assert et2 > et1, f"Period 2 ET ({et2}) should exceed period 1 ({et1})"
    assert np.isclose(et3, et2, rtol=1e-6), (
        f"Period 3 ET ({et3}) should equal period 2 ({et2})"
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
