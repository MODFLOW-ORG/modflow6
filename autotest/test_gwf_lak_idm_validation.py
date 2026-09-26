"""LAK's PACKAGEDATA row-count and PERIOD-setting bound validation.

- A duplicate PACKAGEDATA row and a lake with no PACKAGEDATA row are both
  rejected, with both errors reported together.
- PERIOD outlet-domain settings (RATE, INVERT, WIDTH, SLOPE, ROUGH) must be
  sized and bound-checked against NOUTLETS, not NLAKES. Verifies both
  directions of a NLAKES/NOUTLETS mismatch: an outlet number beyond NLAKES
  must be accepted when NOUTLETS allows it, and an outlet number beyond
  NOUTLETS must still be rejected even when NLAKES would allow it. A model
  with NLAKES == NOUTLETS cannot distinguish either direction, so both cases
  here deliberately use mismatched counts.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-ifno-errors", "lak-noutlets-gt-nlakes", "lak-noutlets-lt-nlakes"]


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


def build_models(idx, test):
    return [build_ifno_errors, build_noutlets_gt, build_noutlets_lt][idx](test)


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


def check_output(idx, test):
    [check_ifno_errors, check_noutlets_gt, check_noutlets_lt][idx](test)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=(idx != 1),
    )
    test.run()
