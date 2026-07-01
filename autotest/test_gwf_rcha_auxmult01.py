"""
Test that RCHA with READASARRAYS and AUXMULTNAME correctly persists the
auxiliary multiplier array across stress periods when it is not re-read,
and that the multiplier is correctly overwritten when re-specified.

Case 1 (rcha_auxmult01): Multiplier set in period 1 only, recharge updated
in period 2. Verifies multiplier persists.

Case 2 (rcha_auxmult02): Multiplier set in period 1, then overridden with a
new value in period 3. Verifies override takes effect.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["rcha_auxmult01", "rcha_auxmult02"]

# spatial discretization
nlay, nrow, ncol = 1, 5, 5
delr = delc = 1000.0
top = 10.0
botm = [0.0]
strt = 5.0
k = 10.0
rech_rate = 0.001


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    if idx == 0:
        # Case 1: multiplier persists (not re-specified after period 1)
        nper = 3
        rech = {
            0: np.full((nrow, ncol), rech_rate),
            1: np.full((nrow, ncol), rech_rate * 2.0),
        }
        rch_aux = {0: [np.full((nrow, ncol), 0.5)]}
    else:
        # Case 2: multiplier overridden in period 3
        nper = 4
        rech = {0: np.full((nrow, ncol), rech_rate)}
        rch_aux = {0: [np.full((nrow, ncol), 0.5)], 2: [np.full((nrow, ncol), 2.0)]}

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

    chd_spd = [[(0, i, ncol - 1), strt] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

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
    rch_budgets = cobj.get_data(text="RCH")
    area = delr * delc

    if idx == 0:
        # Case 1: mult=0.5 persists across all periods
        expected = [
            rech_rate * 0.5 * area,  # period 1
            rech_rate * 2.0 * 0.5 * area,  # period 2: new rech, old mult
            rech_rate * 2.0 * 0.5 * area,  # period 3: both persist
        ]
    else:
        # Case 2: mult=0.5 in periods 1-2, mult=2.0 in periods 3-4
        expected = [
            rech_rate * 0.5 * area,  # period 1
            rech_rate * 0.5 * area,  # period 2: both persist
            rech_rate * 2.0 * area,  # period 3: mult overridden to 2.0
            rech_rate * 2.0 * area,  # period 4: override persists
        ]

    for iper, (rch_data, exp_flux) in enumerate(zip(rch_budgets, expected)):
        actual_max = rch_data["q"].max()
        msg = (
            f"Case {idx + 1}, Period {iper + 1}: "
            f"expected rch flux {exp_flux}, got {actual_max}"
        )
        assert np.isclose(actual_max, exp_flux, rtol=1e-6), msg


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
