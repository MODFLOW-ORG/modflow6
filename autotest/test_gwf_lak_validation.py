"""LAK PACKAGEDATA row-count and PERIOD outlet-setting bound validation.

- A duplicate PACKAGEDATA row and a missing one are both rejected, with
  both errors reported together.
- PERIOD outlet settings (RATE, INVERT, WIDTH, SLOPE, ROUGH) are bound-
  checked against NOUTLETS, not NLAKES, in both directions.
- CONNECTIONDATA's total row count (SUM(NLAKECONN)) is summed across
  lakes with different connection counts, not taken from just one.
"""

import flopy
import numpy as np
import pytest
from compare import Comparison
from framework import TestFramework

cases = [
    "lak-ifno-errors",
    "lak-noutlets-gt-nlakes",
    "lak-noutlets-lt-nlakes",
    "lak-nlakeconn-sum",
]


def build_ifno_errors(test):
    name = cases[0]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 2), 4.0]]
    )

    # nlakes=3: lake 0 given twice (duplicate), lake 1 given once, lake 2
    # never given (missing).
    packagedata = [
        (0, 5.0, 1),
        (0, 5.0, 1),
        (1, 5.0, 1),
    ]
    connectiondata = [
        (0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 5.0, 10.0, 10.0),
        (1, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 5.0, 10.0, 10.0),
    ]
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [(0, "status", "inactive")]},
        pname="lak-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def build_noutlets_gt(test):
    # NLAKES=1, NOUTLETS=2: outlet 2 only exists because of NOUTLETS, not
    # NLAKES -- its RATE setting must be accepted.
    name = cases[1]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="lak-gt", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 2), 100.0]]
    )
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        noutlets=2,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, 100.0, 1)],
        connectiondata=[(0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)],
        outlets=[
            (0, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0),
            (1, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0),
        ],
        perioddata={0: [(0, "status", "active"), (1, "rate", -0.001)]},
        pname="lak-1",
        observations={"lak_outlet.csv": [("out2", "outlet", 2)]},
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def build_noutlets_lt(test):
    # NLAKES=2, NOUTLETS=1: outlet 2 doesn't exist even though lake 2 does
    # -- its RATE setting must still be rejected.
    name = cases[2]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="lak-lt", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 4), 100.0]]
    )
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=2,
        noutlets=1,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, 100.0, 1), (1, 100.0, 1)],
        connectiondata=[
            (0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0),
            (1, 0, (0, 0, 3), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0),
        ],
        outlets=[(0, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0)],
        perioddata={
            0: [
                (0, "status", "active"),
                (1, "status", "active"),
                (1, "rate", -0.001),
            ]
        },
        pname="lak-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def build_nlakeconn_sum(test):
    # NLAKES=2, NLAKECONN 3 and 5 -- a wrong SUM would mis-size
    # CONNECTIONDATA's 8 total rows against the two lakes' actual counts.
    name = cases[3]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="lak-sum", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=9, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 8), 100.0]]
    )
    lake0 = [
        (0, j, (0, 0, j), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0) for j in range(3)
    ]
    lake1 = [
        (1, j, (0, 0, j + 3), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)
        for j in range(5)
    ]
    # interleaved, not grouped by lake -- row order must not matter
    connectiondata = [
        lake1[0],
        lake0[0],
        lake1[1],
        lake0[1],
        lake1[2],
        lake0[2],
        lake1[3],
        lake1[4],
    ]
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=2,
        save_flows=True,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, 100.0, 3), (1, 100.0, 5)],
        connectiondata=connectiondata,
        perioddata={0: [(0, "status", "active"), (1, "status", "active")]},
        pname="lak-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim


def build_models(idx, test):
    return [
        build_ifno_errors,
        build_noutlets_gt,
        build_noutlets_lt,
        build_nlakeconn_sum,
    ][idx](test)


def check_ifno_errors(test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert "data for lake 1 specified 2 times" in lst, (
        "expected the duplicate PACKAGEDATA row to be flagged; got:\n" + lst
    )
    assert "no data specified for lake 3" in lst, (
        "expected the missing lake to be flagged; got:\n" + lst
    )


def check_noutlets_gt(test):
    obs = np.genfromtxt(test.workspace / "lak_outlet.csv", delimiter=",", names=True)
    assert np.isclose(float(obs["OUT2"]), -0.001), (
        "outlet 2 (valid under NOUTLETS=2, invalid under NLAKES=1) should "
        "have applied its specified RATE"
    )


def check_noutlets_lt(test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert "NUMBER" in lst and "less than or equal to 1" in lst, (
        "outlet 2 (invalid under NOUTLETS=1, valid under NLAKES=2) should "
        "have been rejected against NOUTLETS; got:\n" + lst
    )


def check_nlakeconn_sum(test):
    name = cases[3]
    bobj = flopy.utils.CellBudgetFile(test.workspace / f"{name}.lak.cbc")
    records = bobj.get_data(text="GWF")[0]
    assert len(records) == 8, f"expected 8 GWF connections (3+5), got {len(records)}"
    assert (records["node"][:3] == 1).all() and (records["node"][3:] == 2).all(), (
        "connections were not grouped under the correct lake"
    )


def check_output(idx, test):
    [check_ifno_errors, check_noutlets_gt, check_noutlets_lt, check_nlakeconn_sum][idx](
        test
    )


@pytest.mark.parametrize(
    "compare",
    [
        "auto",
        pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4),
    ],
)
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, compare, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=(idx not in (1, 3)),
        compare=compare,
    )
    test.run()
