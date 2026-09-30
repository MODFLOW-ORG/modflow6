"""UZF PERIOD-block input errors, rejected as clean input errors with the
erroring file reported:

  per_dup    (uzf-per-dup): a PERIOD block with more rows than NUZFCELLS
      (here, the same cell given twice in one period).
  ifno_range (uzf-ifno-range): an out-of-range IFNO, with a distinct
      message per bad row (value and row number).
"""

import flopy
import pytest
from framework import DNODATA, TestFramework

cases = ["uzf-per-dup", "uzf-ifno-range"]


def _get_per_dup(ws, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)

    packagedata = [(0, (0, 0, 0), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5)]
    perioddata = {
        0: [
            (0, 0.01, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
            (0, 0.02, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
        ],
    }
    flopy.mf6.ModflowGwfuzf(
        gwf,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _check_per_dup(test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert "line count exceeds input dimension" in lst, (
        "expected the excess PERIOD row to be flagged; got:\n" + lst
    )


def _get_ifno_range(ws, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
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

    # nuzfcells=2: both cells properly defined in PACKAGEDATA. PERIOD
    # reissues cell 0 (valid) then cell 5 -- out of range, doesn't exist.
    packagedata = [
        (0, (0, 0, 0), 1, -1, 0.5, 0.1, 0.2, 0.3, 0.25, 3.5),
        (1, (0, 0, 1), 1, -1, 0.5, 0.1, 0.2, 0.3, 0.25, 3.5),
    ]
    flopy.mf6.ModflowGwfuzf(
        gwf,
        nuzfcells=2,
        packagedata=packagedata,
        perioddata={
            0: [
                (0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                (5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ]
        },
        pname="uzf-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _check_ifno_range(test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert (
        "IFNO 6 on row 2 must be greater than 0 and less than or equal to 2." in lst
    ), "expected a distinct, row-and-value-specific IFNO range message; got:\n" + lst
    assert "ERROR OCCURRED WHILE READING FILE" in lst, (
        "expected the erroring input file to be reported; got:\n" + lst
    )


_builders = [_get_per_dup, _get_ifno_range]
_checkers = [_check_per_dup, _check_ifno_range]


def build_models(idx, test):
    name = cases[idx]
    return _builders[idx](test.workspace, name)


def check_output(idx, test):
    _checkers[idx](test)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=True,
    )
    test.run()
