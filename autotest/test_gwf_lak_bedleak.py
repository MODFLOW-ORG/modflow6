import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["bedleak", "bedleak_fail", "bedleak_none"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 10, 10
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]

    lenx = 300.0
    delr = delc = lenx / float(nrow)
    strt = 100.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

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

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=90.0,
        botm=0.0,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=1.0, k33=0.01)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=0.1,
        steady_state={0: True},
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 100.0])
    chdlist0.append([(0, nrow - 1, ncol - 1), 95.0])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
    )

    # lak package
    if "fail" in name:
        bedleak = -100.0
    elif "none" in name:
        bedleak = "none"
    else:
        bedleak = DNODATA

    # <ifno> <strt> <nlakeconn> [<aux(naux)>] [<boundname>]
    packagedata = [
        [0, 100.0, 1, "lake1"],
        [1, 100.0, 1, "lake2"],
    ]
    # <ifno> <iconn> <cellid(ncelldim)> <claktype> <bedleak> <belev> <telev> ...
    #        <connlen> <connwidth>
    connectiondata = [
        [0, 0, (0, 1, 1), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [1, 0, (0, 2, 2), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
    ]
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        boundnames=True,
        surfdep=1.0,
        print_input=True,
        print_stage=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord=f"{name}.lak.bud",
        nlakes=len(packagedata),
        packagedata=packagedata,
        connectiondata=connectiondata,
    )
    # lak.remove()

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]

    # lak budget
    if "fail" not in name:
        fpth = os.path.join(test.workspace, f"{name}.lak.bud")
        bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        bobj.list_unique_records()
        records = bobj.get_data(text="GWF")
        for r in records:
            assert np.allclose(r["q"][0], -4.79616347e-12)
            assert np.allclose(r["q"][1], -6.19237994e-12)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail="fail" in str(function_tmpdir),
    )
    test.run()
