"""
Test that RCHA with READASARRAYS and AUXMULTNAME correctly persists the
auxiliary multiplier array across stress periods when it is not re-read.

The multiplier is specified only in period 1. Recharge is updated in period 2
without re-specifying the multiplier. The test verifies that the multiplier
from period 1 continues to be applied in period 2 (i.e., recharge is not zero).
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["rcha_auxmult01"]

# temporal discretization
nper = 3
perlen = [1.0, 1.0, 1.0]
nstp = [1, 1, 1]
tsmult = [1.0, 1.0, 1.0]

# spatial discretization
nlay, nrow, ncol = 1, 5, 5
delr = delc = 1000.0
top = 10.0
botm = [0.0]
strt = 5.0

# npf
k = 10.0

# recharge - provided in periods 1 and 2
rech_rate = 0.001
rech1 = np.full((nrow, ncol), rech_rate)
rech2 = np.full((nrow, ncol), rech_rate * 2.0)
rech = {0: rech1, 1: rech2}

# auxiliary multiplier - only provided in period 1
rch_mult = np.full((nrow, ncol), 0.5)
rch_aux = {0: [rch_mult]}


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)

    tdis_rc = [(p, n, t) for p, n, t in zip(perlen, nstp, tsmult)]
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc)

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

    # constant head on one side to allow flow
    chd_spd = [[(0, i, ncol - 1), strt] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

    # rcha with auxmultname - multiplier only in period 1
    flopy.mf6.ModflowGwfrcha(
        gwf,
        readasarrays=True,
        recharge=rech,
        auxiliary="RCH_MULT",
        auxmultname="RCH_MULT",
        aux=rch_aux,
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

    # get recharge for each stress period
    rch_budgets = cobj.get_data(text="RCH")

    # period 1: recharge * mult = 0.001 * 0.5 = 0.0005
    # period 2: recharge * mult = 0.002 * 0.5 = 0.001 (mult persists)
    # period 3: same as period 2 (no new data)
    expected_q = [
        rech_rate * 0.5 * delr * delc,
        rech_rate * 2.0 * 0.5 * delr * delc,
        rech_rate * 2.0 * 0.5 * delr * delc,
    ]

    for iper, (rch_data, exp_flux) in enumerate(zip(rch_budgets, expected_q)):
        q = rch_data["q"]
        actual_max = q.max()
        msg = f"Period {iper + 1}: expected rch flux {exp_flux}, got {actual_max}"
        assert np.isclose(actual_max, exp_flux, rtol=1e-6), msg
        # verify no zero recharge (the bug symptom)
        if iper > 0:
            assert actual_max > 0.0, (
                f"Period {iper + 1}: recharge is zero - "
                f"AUXMULTNAME not persisting across periods"
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
