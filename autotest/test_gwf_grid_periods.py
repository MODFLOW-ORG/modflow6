"""
Period input of RCHA and WELG: a missing period carries the previous input
forward, an empty grid period clears its stresses, and input carries forward
to NPER. The well rate is scaled by the second aux variable (AUXMULTNAME).
"""

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["grid-periods"]
nper = 5
nlay, nrow, ncol = 2, 1, 5
delr = delc = 10.0
# recharge area: the constant-head cell takes none
area = (ncol - 1) * delr * delc
# recharge rate and well aux multiplier given in these periods only
recharge = {0: 0.01, 2: 0.02}
auxmult = {0: 0.5, 3: 2.0}
qwell = -1.0
# expected rates, by period
expected_rch = [0.01, 0.01, 0.02, 0.02, 0.02]
expected_wel = [-0.5, 0.0, 0.0, -2.0, -2.0]


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=[-10.0, -20.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [((0, 0, 0), 0.0)]})
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    q = np.full((nlay, nrow, ncol), DNODATA)
    q[1, 0, ncol - 1] = qwell
    ones = np.ones((nlay, nrow, ncol))
    flopy.mf6.ModflowGwfwelg(
        gwf,
        auxiliary=["a1", "a2"],
        auxmultname="a2",
        # period 2 is empty: no well from period 2 until period 4
        q={0: q, 1: [], 3: q},
        aux={k: [ones, ones * m] for k, m in auxmult.items()},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    cbc = flopy.utils.CellBudgetFile(test.workspace / f"{name}.cbc", precision="double")
    for kper in range(nper):
        rch = sum(r["q"].sum() for r in cbc.get_data(text="RCH", kstpkper=(0, kper)))
        wel = sum(r["q"].sum() for r in cbc.get_data(text="WEL", kstpkper=(0, kper)))
        assert np.isclose(rch, expected_rch[kper] * area), (
            f"period {kper + 1}: RCH {rch}"
        )
        assert np.isclose(wel, expected_wel[kper]), f"period {kper + 1}: WEL {wel}"


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
