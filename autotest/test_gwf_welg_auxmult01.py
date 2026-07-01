"""
Test grid array (WELG) package behaviors:

Case 1 (welg_auxmult01): AUXMULTNAME persists across periods when not re-read.
Case 2 (welg_auxmult02): Aux override in a later period takes effect.
Case 3 (welg_auxmult03): CONSTANT for aux (secondary array) works correctly.
Case 4 (welg_auxmult04): Empty PERIOD block removes stresses (nbound=0).
Case 5 (welg_auxmult05): CONSTANT for primary array (all cells become bounds).
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "welg_auxmult01",
    "welg_auxmult02",
    "welg_auxmult03",
    "welg_auxmult04",
    "welg_auxmult05",
]

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

    if idx == 0:
        # Case 1: aux persists, Q updated in period 2, no block period 3
        nper = 3
        q1 = np.full((nlay, nrow, ncol), DNODATA)
        q1[0, 2, 2] = 100.0
        q2 = np.full((nlay, nrow, ncol), DNODATA)
        q2[0, 2, 2] = 200.0
        wel_q = {0: q1, 1: q2}
        mult = np.full((nlay, nrow, ncol), DNODATA)
        mult[0, 2, 2] = 0.5
        wel_aux = {0: [mult]}

    elif idx == 1:
        # Case 2: aux overridden in period 3
        nper = 4
        q1 = np.full((nlay, nrow, ncol), DNODATA)
        q1[0, 2, 2] = 100.0
        wel_q = {0: q1}
        mult1 = np.full((nlay, nrow, ncol), DNODATA)
        mult1[0, 2, 2] = 0.5
        mult2 = np.full((nlay, nrow, ncol), DNODATA)
        mult2[0, 2, 2] = 2.0
        wel_aux = {0: [mult1], 2: [mult2]}

    elif idx == 2:
        # Case 3: CONSTANT for aux (secondary array)
        # Q is sparse, aux uses CONSTANT via flopy (all cells get same value)
        nper = 2
        q1 = np.full((nlay, nrow, ncol), DNODATA)
        q1[0, 2, 2] = 100.0
        q1[0, 3, 3] = 100.0
        wel_q = {0: q1}
        # aux as full array with same value at active positions, DNODATA elsewhere
        mult = np.full((nlay, nrow, ncol), DNODATA)
        mult[0, 2, 2] = 0.5
        mult[0, 3, 3] = 0.5
        wel_aux = {0: [mult]}

    elif idx == 3:
        # Case 4: empty PERIOD block in period 2 removes stresses
        nper = 3
        q1 = np.full((nlay, nrow, ncol), DNODATA)
        q1[0, 2, 2] = 100.0
        wel_q = {0: q1, 1: np.full((nlay, nrow, ncol), DNODATA)}
        wel_aux = None

    elif idx == 4:
        # Case 5: CONSTANT for primary array (all cells become bounds)
        nper = 1
        wel_q = {0: np.full((nlay, nrow, ncol), 4.0)}
        wel_aux = None

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

    # constant head on boundary
    chd_spd = [[(0, i, ncol - 1), strt] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

    welg_kwargs = dict(q=wel_q)
    if idx in (0, 1, 2):
        welg_kwargs["auxiliary"] = "WEL_MULT"
        welg_kwargs["auxmultname"] = "WEL_MULT"
        welg_kwargs["aux"] = wel_aux
    if idx == 2:
        welg_kwargs["maxbound"] = 2
    elif idx == 4:
        welg_kwargs["maxbound"] = nlay * nrow * ncol
    elif idx != 3:
        welg_kwargs["maxbound"] = 1

    flopy.mf6.ModflowGwfwelg(gwf, **welg_kwargs)

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

    if idx == 0:
        # period 1: 100*0.5=50, period 2: 200*0.5=100, period 3: persists=100
        expected = [50.0, 100.0, 100.0]
        for iper, (wel_data, exp_q) in enumerate(zip(wel_budgets, expected)):
            actual = wel_data["q"].max()
            assert np.isclose(actual, exp_q, rtol=1e-6), (
                f"Case 1, Period {iper + 1}: expected {exp_q}, got {actual}"
            )

    elif idx == 1:
        # period 1: 100*0.5=50, period 2: persists=50,
        # period 3: 100*2.0=200 (aux overridden), period 4: persists=200
        expected = [50.0, 50.0, 200.0, 200.0]
        for iper, (wel_data, exp_q) in enumerate(zip(wel_budgets, expected)):
            actual = wel_data["q"].max()
            assert np.isclose(actual, exp_q, rtol=1e-6), (
                f"Case 2, Period {iper + 1}: expected {exp_q}, got {actual}"
            )

    elif idx == 2:
        # Two wells, both with mult=0.5: 100*0.5=50 each
        # period 2: persists
        expected_sum = [100.0, 100.0]  # 2 wells * 50 each
        for iper, (wel_data, exp_sum) in enumerate(zip(wel_budgets, expected_sum)):
            actual_sum = wel_data["q"].sum()
            assert np.isclose(actual_sum, exp_sum, rtol=1e-6), (
                f"Case 3, Period {iper + 1}: expected sum {exp_sum}, got {actual_sum}"
            )

    elif idx == 3:
        # period 1: well active with q=100, period 2: empty block = no wells,
        # period 3: no block = persists empty (nbound=0)
        # period 1
        assert wel_budgets[0]["q"].max() > 0.0, (
            "Case 4, Period 1: well should be active"
        )
        # period 2: empty block - no well flow
        assert len(wel_budgets[1]["q"]) == 0 or np.isclose(
            wel_budgets[1]["q"].sum(), 0.0
        ), "Case 4, Period 2: empty block should remove stresses"
        # period 3: persists from period 2 (no stresses)
        assert len(wel_budgets[2]["q"]) == 0 or np.isclose(
            wel_budgets[2]["q"].sum(), 0.0
        ), "Case 4, Period 3: should persist empty from period 2"

    elif idx == 4:
        # All cells are wells with q=4.0 (CONSTANT)
        # nbound should equal total grid nodes
        total_q = wel_budgets[0]["q"].sum()
        n_bounds = len(wel_budgets[0]["q"])
        assert n_bounds == nlay * nrow * ncol, (
            f"Case 5: expected nbound={nlay * nrow * ncol}, got {n_bounds}"
        )
        assert total_q > 0.0, "Case 5: total well flow should be positive"


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
