"""A duplicate PACKAGEDATA row and a well with no PACKAGEDATA row are both
rejected, with both errors reported together.
"""

import flopy
import pytest
from framework import TestFramework

cases = ["maw-ifno-errors"]


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[-100.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 2), 4.0]]
    )

    # nmawwells=3: well 0 given twice (duplicate), well 1 given once, well 2
    # never given (missing).
    packagedata = [
        (0, 0.15, -100.0, 5.0, "thiem", 1),
        (0, 0.15, -100.0, 5.0, "thiem", 1),
        (1, 0.15, -100.0, 5.0, "thiem", 1),
    ]
    connectiondata = [
        (0, 0, (0, 0, 1), 10.0, -100.0, 1.0, 0.25),
        (1, 0, (0, 0, 1), 10.0, -100.0, 1.0, 0.25),
    ]
    flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [(0, "status", "inactive")]},
        pname="maw-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def check_output(idx, test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert "data for well 1 specified 2 times" in lst, (
        "expected the duplicate PACKAGEDATA row to be flagged; got:\n" + lst
    )
    assert "no data specified for well 3" in lst, (
        "expected the missing well to be flagged; got:\n" + lst
    )


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
