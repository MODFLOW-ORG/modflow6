# test use of binary stress period data with auxiliary data
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "chd02",
]


def build_models(idx, test):
    name = cases[idx]
    nlay, nrow, ncol = 1, 1, 10
    sim = flopy.mf6.MFSimulation(sim_ws=test.workspace, sim_name=name)
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, complexity="simple")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, print_input=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=10.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
    )
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=10.0,
    )
    chd_data = [
        (0, 0, 0, 10.0, 1.0, 100.0),
        (0, 0, ncol - 1, 5.0, 0.0, 100.0),
    ]
    chd_data = {
        0: {
            "filename": "chd.bin",
            "binary": True,
            "data": chd_data,
        },
    }
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary=["conc", "something"],
        stress_period_data=chd_data,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    name = test.name

    fpth = os.path.join(test.workspace, f"{name}.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    head = hobj.get_data().flatten()

    # This is the answer to this problem.
    hres = np.array(
        [
            10.00,
            9.57481441,
            9.1298034,
            8.66189866,
            8.16714512,
            7.64029169,
            7.07409994,
            6.45808724,
            5.77600498,
            5.000,
        ]
    )
    assert np.allclose(hres, head), "simulated head does not match with known solution."


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
