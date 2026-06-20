"""
Test that WELG (grid array package) with AUXMULTNAME correctly persists
the auxiliary multiplier array across stress periods when it is not re-read.

The multiplier is specified only in period 1. The well rate is updated in
period 2 without re-specifying the multiplier. The test verifies that the
multiplier from period 1 continues to be applied.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["welg_auxmult01"]

# spatial discretization
nlay, nrow, ncol = 1, 5, 5
delr = delc = 1000.0
top = 10.0
botm = [0.0]
strt = 5.0
k = 10.0

DNODATA = 3.0e30


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nper = 3

    # well rate arrays - injection in cell (0,2,2) only
    q1 = np.full((nlay, nrow, ncol), DNODATA)
    q1[0, 2, 2] = 100.0
    q2 = np.full((nlay, nrow, ncol), DNODATA)
    q2[0, 2, 2] = 200.0
    wel_q = {0: q1, 1: q2}

    # multiplier array - only in period 1
    mult = np.full((nlay, nrow, ncol), DNODATA)
    mult[0, 2, 2] = 0.5
    wel_aux = {0: [mult]}

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

    # constant head on boundary to balance injection
    chd_spd = [[(0, i, ncol - 1), strt] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

    flopy.mf6.ModflowGwfwelg(
        gwf,
        maxbound=1,
        q=wel_q,
        auxiliary="WEL_MULT",
        auxmultname="WEL_MULT",
        aux=wel_aux,
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
    wel_budgets = cobj.get_data(text="WEL")

    # period 1: q=100 * mult=0.5 = 50
    # period 2: q=200 * mult=0.5 = 100 (mult persists)
    # period 3: same as period 2 (no new data)
    expected = [50.0, 100.0, 100.0]

    for iper, (wel_data, exp_q) in enumerate(zip(wel_budgets, expected)):
        q = wel_data["q"]
        actual_max = q.max()
        msg = f"Period {iper + 1}: expected wel flux {exp_q}, got {actual_max}"
        assert np.isclose(actual_max, exp_q, rtol=1e-6), msg
        assert actual_max > 0.0, (
            f"Period {iper + 1}: well rate is zero - "
            f"AUXMULTNAME not persisting across periods"
        )


@pytest.mark.developmode
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
