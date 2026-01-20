"""
Test adaptive time step module
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_sto02a", "gwf_sto02b"]
ncols = [1, 2]
(nlay, nrow) = (1, 1)


def build_models(idx, test):
    perlen = [10]
    nper = len(perlen)
    nstp = [1]
    tsmult = nper * [1.0]
    delr = 100.0
    delc = 1.0
    top = 100.0
    botm = [0.0]
    strt = top
    hk = 1.0
    laytyp = 1
    ss = 0.0
    sy = 0.1

    tdis_rc = []
    for id in range(nper):
        tdis_rc.append((perlen[id], nstp[id], tsmult[id]))

    name = cases[idx]
    ncol = ncols[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions=None,
    )

    # create iterative model solution and register the gwf model with it
    nouter, ninner = 20, 10
    hclose, rclose, relax = 1e-6, 1e-6, 0.97
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        # under_relaxation="DBD",
        # under_relaxation_theta=0.7,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state={0: False},
        transient={0: True},
    )

    # wel files
    welspdict = {
        0: [[(0, 0, j), -10.0] for j in range(ncol)],
    }
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=welspdict,
        save_flows=False,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # write MODFLOW 6 files
    sim.write_simulation()

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwfname = "gwf_" + name

    # This will fail if budget numbers cannot be read
    fpth = os.path.join(test.workspace, f"{gwfname}.lst")
    mflist = flopy.utils.Mf6ListBudget(fpth)
    names = mflist.get_record_names()
    print(names)

    inc = mflist.get_incremental()
    print(inc)


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
