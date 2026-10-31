"""
Period input of RCHA and WELG: a missing period carries the previous input
forward, an empty grid period clears its stresses, and input carries forward
to NPER. The well rate is scaled by the second aux variable (AUXMULTNAME).

WELG's q and aux are each given in different periods (3 and 4): q restates
at period 3 without aux restating there, and aux restates at period 4
without q restating there. Each field fills forward on its own -- one
field restating a period does not clear or reset a sibling field left
silent at that same period.
"""

import flopy
import numpy as np
import pytest
from compare import Comparison
from framework import DNODATA, TestFramework

cases = ["grid-periods"]
nper = 5
nlay, nrow, ncol = 2, 1, 5
delr = delc = 10.0
# recharge area: the constant-head cell takes none
area = (ncol - 1) * delr * delc
# recharge rate and well aux multiplier given in these periods only
recharge = {0: 0.01, 2: 0.02}
auxmult = {0: 0.5, 4: 3.0}
qwell = -1.0
# expected rates, by period
expected_rch = [0.01, 0.01, 0.02, 0.02, 0.02]
# period 3: q restates (-1.0) while aux's mult (0.5, from period 0) stays
# silent and still applies. period 4: aux restates (3.0) while q (-1.0,
# from period 3) stays silent and still applies.
expected_wel = [-0.5, 0.0, 0.0, -0.5, -3.0]


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
    # aux must share q's node pattern (DNODATA elsewhere): when aux is the
    # only array restated in a period (period 4, below), it's the one
    # that establishes which cells are active, per MF6's IDM loader
    # (GridArrayLoadType%rp/param_load in Mf6FileGridArray.f90) -- a
    # dense array would wrongly mark every cell active instead of just
    # the well's.
    a1 = np.full((nlay, nrow, ncol), DNODATA)
    a1[1, 0, ncol - 1] = 1.0

    def _mult(m):
        a = np.full((nlay, nrow, ncol), DNODATA)
        a[1, 0, ncol - 1] = m
        return a

    flopy.mf6.ModflowGwfwelg(
        gwf,
        auxiliary=["a1", "a2"],
        auxmultname="a2",
        # period 2 is empty: no well from period 2 until period 4
        q={0: q, 1: [], 3: q},
        aux={k: [a1, _mult(m)] for k, m in auxmult.items()},
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
@pytest.mark.parametrize(
    "compare",
    [
        "auto",
        pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4),
        pytest.param(
            Comparison.FP4_STRUCTURED, marks=[pytest.mark.netcdf, pytest.mark.fp4]
        ),
        pytest.param(
            Comparison.FP4_LAYERED, marks=[pytest.mark.netcdf, pytest.mark.fp4]
        ),
    ],
)
def test_mf6model(idx, name, compare, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
        targets=targets,
    )
    test.run()
