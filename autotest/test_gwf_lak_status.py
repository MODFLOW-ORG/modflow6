"""
Test to make sure lakes can be turned on and off using
the STATUS setting.

This is the configuration, where C is constant head
and L is the lake.
C 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 L L L 1 1 1 1
1 1 1 L L L 1 1 1 1
1 1 1 L L L 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 C

The lake is active for the first and last stress
period and inactive for the second stress period.

The test checks the lake observations, lake stage
file, lake budget file, and gwf head file.

"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf-lak-status"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 10, 10
    nper = 3
    perlen = nper * [1.0]
    nstp = nper * [1]
    tsmult = nper * [1.0]

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

    # <ifno> <strt> <nlakeconn> [<aux(naux)>] [<boundname>]
    packagedata = [
        [0, 100.0, 9, "lake1"],
    ]
    # <ifno> <iconn> <cellid(ncelldim)> <claktype> <bedleak> <belev> <telev> ...
    #        <connlen> <connwidth>
    bedleak = 1.0
    connectiondata = [
        [0, 0, (0, 3, 3), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 1, (0, 3, 4), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 2, (0, 3, 5), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 3, (0, 4, 3), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 4, (0, 4, 4), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 5, (0, 4, 5), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 6, (0, 5, 3), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 7, (0, 5, 4), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
        [0, 8, (0, 5, 5), "vertical", bedleak, 0.0, 0.0, 0.0, 0.0],
    ]
    lak_spd = {
        0: [[0, "rainfall", 0.1]],
        1: [[0, "status", "inactive"]],
        2: [[0, "status", "active"]],
    }

    # note: for specifying lake number, use fortran indexing!
    fname = f"{name}.lak.obs.csv"
    lak_obs = {
        fname: [
            ("lakestage", "stage", 1),
            ("lakevolume", "volume", 1),
            ("lak1", "lak", "lake1"),
        ],
        "digits": 10,
    }
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        boundnames=True,
        surfdep=1.0,
        print_input=True,
        print_stage=True,
        print_flows=True,
        save_flows=True,
        stage_filerecord=f"{name}.lak.stage",
        budget_filerecord=f"{name}.lak.bud",
        nlakes=len(packagedata),
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=lak_spd,
        observations=lak_obs,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_budget_file(idx, test):
    # lak budget
    name = cases[idx]
    fpth = test.workspace / f"{name}.lak.bud"
    print(f"Checking contents of {fpth.name}")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    times = bobj.get_times()
    print(bobj.list_unique_records())

    # check to make sure GWF is zero when lake is inactive
    # and non-zero otherwise
    for kper in range(3):
        print(f"  Checking binary budget for stress period {kper + 1}")
        record = bobj.get_data(text="GWF", totim=times[kper])[0]
        for r in record:
            if kper == 1:
                assert r["q"] == 0.0
            else:
                assert r["q"] != 0.0

    # check all the other terms to make sure they are zero for
    # the second period
    lake_terms = [
        "rainfall",
        "evaporation",
        "ext-inflow",
        "withdrawal",
        "ext-outflow",
        "storage",
        "constant",
    ]
    for lake_term in lake_terms:
        print(f"  Checking lake term {lake_term} for period 2.")
        record = bobj.get_data(text=lake_term, totim=times[1])[0]
        q = record["q"][0]
        msg = f"Lake term for stress period 2 is not zero ({lake_term}={q})"
        assert np.allclose(q, 0.0), msg


def check_stage_file(idx, test):
    # Check to make sure the stage is equal to dhnoflo
    # for stress period 2 and not dhnoflo otherwise
    dhnoflo = 1.0e30
    name = cases[idx]
    fpth = test.workspace / f"{name}.lak.stage"
    print(f"Checking contents of {fpth.name}")
    bobj = flopy.utils.HeadFile(fpth, text="stage", precision="double")
    times = bobj.get_times()
    for kper in range(3):
        print(f"  Checking binary stage for stress period {kper + 1}")
        stage = bobj.get_data(totim=times[kper]).flatten()
        print(stage)
        if kper == 1:
            assert stage[0] == dhnoflo
        else:
            assert stage[0] != dhnoflo


def check_head_file(idx, test):
    # Check to make sure the simulated head for the first and
    # last period is the same and different from the second period.
    name = cases[idx]
    fpth = test.workspace / f"{name}.hds"
    print(f"Checking contents of {fpth.name}")
    bobj = flopy.utils.HeadFile(fpth, text="head", precision="double")
    times = bobj.get_times()

    head0 = bobj.get_data(totim=times[0]).flatten()
    head1 = bobj.get_data(totim=times[1]).flatten()
    head2 = bobj.get_data(totim=times[2]).flatten()

    assert np.allclose(head0, head2), (
        "Simulated heads for period 1 and 3 should be the same"
    )
    assert np.all(np.less_equal(head1, 100.0)), (
        f"Simulated heads for period 2 should be less than or equal to 100. {head1}"
    )
    assert np.any(np.greater(head0, 100.0)), (
        f"Some simulated heads for period 1 should be greater than 100. {head0}"
    )


def check_lake_obs(idx, test):
    # Check to make sure the simulated head for the first and
    # last period is the same and different from the second period.
    dnodata = 3.0e30
    name = cases[idx]
    fpth = test.workspace / f"{name}.lak.obs.csv"
    print(f"Checking contents of {fpth.name}")
    obs = np.genfromtxt(fpth, names=True, delimiter=",")
    stage = obs["LAKESTAGE"].tolist()
    print(stage)
    assert stage[0] == stage[2], "Period 1 and period 3 stages should be equal."
    assert stage[1] == dnodata, "Period 2 stage should equal dnodata"


def check_output(idx, test):
    check_lake_obs(idx, test)
    check_head_file(idx, test)
    check_stage_file(idx, test)
    check_budget_file(idx, test)


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
