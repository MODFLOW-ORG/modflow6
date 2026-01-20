"""
Test to make sure that NPF anisotropy ratio options are read and processed
correctly.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["npf05a"]


def build_models(idx, test):
    nlay, nrow, ncol = 2, 1, 5
    chdheads = [100.0]
    nper = len(chdheads)
    perlen = nper * [1.0]
    nstp = nper * [1]
    tsmult = nper * [1.0]

    delr = delc = 1.0
    strt = 100.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "npf"

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=100.0,
        botm=[50.0, 0.0],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    hk = 5.0
    k22 = 0.5
    k33 = 0.05
    k22overk = False
    k33overk = False
    aniso = True
    if aniso:
        k22 = k22 / hk
        k33 = k33 / hk
        k22overk = True
        k33overk = True
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        k22overk=k22overk,
        k33overk=k33overk,
        save_flows=True,
        icelltype=0,
        k=hk,
        k22=k22,
        k33=k33,
    )

    # chd files
    chdspd = {}
    for kper, chdval in enumerate(chdheads):
        chdspd[kper] = [[(0, 0, 0), 100], [(nlay - 1, 0, ncol - 1), 110.0]]
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

    rcha = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.01)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, "npf.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    heads = hobj.get_alldata()
    # answer was obtained from running problem without anisotropy
    answer = [
        100.0,
        100.00031999,
        100.00055998,
        100.00071997,
        100.00079997,
        109.99960002,
        109.99964002,
        109.99972002,
        109.99984001,
        110.0,
    ]
    answer = np.array(answer)
    assert np.allclose(heads.flatten(), answer)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
