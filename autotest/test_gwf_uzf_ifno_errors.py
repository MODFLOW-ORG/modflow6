"""A duplicate PACKAGEDATA row and a cell with no PACKAGEDATA row are both
rejected, with both errors reported together.
"""

import flopy
import pytest
from framework import TestFramework

cases = ["uzf-ifno-errors"]


def build_models(idx, test):
    name = cases[idx]
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

    # nuzfcells=3: cell 0 given twice (duplicate), cell 1 given once, cell 2
    # never given (missing).
    packagedata = [
        (0, (0, 0, 0), 1, -1, 0.5, 0.1, 0.2, 0.3, 0.25, 3.5),
        (0, (0, 0, 0), 1, -1, 0.5, 0.1, 0.2, 0.3, 0.25, 3.5),
        (1, (0, 0, 1), 1, -1, 0.5, 0.1, 0.2, 0.3, 0.25, 3.5),
    ]
    flopy.mf6.ModflowGwfuzf(
        gwf,
        nuzfcells=3,
        packagedata=packagedata,
        perioddata={0: [(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]},
        pname="uzf-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def check_output(idx, test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    assert "data for uzf cell 1 specified 2 times" in lst, (
        "expected the duplicate PACKAGEDATA row to be flagged; got:\n" + lst
    )
    assert "no data specified for uzf cell 3" in lst, (
        "expected the missing cell to be flagged; got:\n" + lst
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
