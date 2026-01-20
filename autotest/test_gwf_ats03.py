"""
Purpose of this test is to ensure that time series values are being
recalculated for the ATS failed retry option.  The first time step fails
multiple times until it can converge.  Each time it fails, it has to back
up in time, which requires the constant head time series to be reinterpolated
for the new totim value.  For this problem, the constant head starts at 100 at
time zero and drops to 50.0 at time 100.  So the constant head values, which
are observed and written to and obs output file must fall on a line between
(0, 100) and (100, 50), which is ensured by this test.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_ats03a"]
nlay, nrow, ncol = 1, 1, 10


def build_models(idx, test):
    perlen = [100.0]
    nper = len(perlen)
    nstp = [1]
    tsmult = nper * [1.0]
    delr = 100.0
    delc = 1.0
    top = 100.0
    botm = [0.0]
    strt = 100.0
    hk = 1.0
    laytyp = 1
    ss = 0.0
    sy = 0.1

    tdis_rc = []
    for id in range(nper):
        tdis_rc.append((perlen[id], nstp[id], tsmult[id]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    # set dt0, dtmin, dtmax, dtadj, dtfailadj
    dt0 = 100.0
    dtmin = 1.00e-5
    dtmax = 100.0
    dtadj = 2.0
    dtfailadj = 5.0
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )
    if True:
        ats_filerecord = name + ".ats"
        atsperiod = [(0, dt0, dtmin, dtmax, dtadj, dtfailadj)]
        tdis.ats.initialize(
            maxats=len(atsperiod),
            perioddata=atsperiod,
            filename=ats_filerecord,
        )

    # create gwf model
    gwfname = name
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        #            newtonoptions=newtonoptions,
    )

    # create iterative model solution and register the gwf model with it
    nouter, ninner = 15, 5
    hclose, rclose, relax = 1.5e-6, 1e-6, 0.97
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
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
        0: [[(0, 0, 0), -1.0]],
    }
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=welspdict,
        save_flows=False,
    )

    # chd files
    chdspdict = {
        0: [[(0, 0, ncol - 1), "tshead"]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=chdspdict,
        save_flows=False,
    )

    # chd ts package
    ts_recarray = [
        (0.0, 100.0),
        (100.0, 50),
    ]
    filename = name + ".chd.ts"
    time_series_namerecord = [("tshead",)]
    interpolation_methodrecord = [("linearend",)]
    chd.ts.initialize(
        filename=filename,
        timeseries=ts_recarray,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
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

    obs_lst = []
    obs_lst.append(["obs1", "head", (0, 0, 0)])
    obs_lst.append(["obs2", "head", (0, 0, ncol - 1)])
    obs_dict = {f"{gwfname}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    return sim, None


def check_output(idx, test):
    # ensure obs2 (a constant head time series) drops linearly from 100 to 50
    fpth = os.path.join(test.workspace, test.name + ".obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    x = np.array(tc["time"])
    answer = 100.0 - x * 0.5
    result = np.array(tc["OBS2"])
    msg = f"obs2 must drop linearly from 100 down to 50: {result}"
    assert np.allclose(answer, result), msg


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
