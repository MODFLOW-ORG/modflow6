"""
Test adaptive time step module
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_ats01a"]
nlay, nrow, ncol = 1, 1, 2

# set dt0, dtmin, dtmax, dtadj, dtfailadj
dt0 = 5
dtmin = 1.001e-5
dtmax = 10.0
dtadj = 2.0
dtfailadj = 5.0


def build_models(idx, test):
    perlen = [10]
    nper = len(perlen)
    nstp = [1]
    tsmult = nper * [1.0]
    delr = 100.0
    delc = 1.0
    top = 100.0
    botm = [0.0]
    strt = 50.0
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
    ats_filerecord = None
    tdis = flopy.mf6.ModflowTdis(
        sim,
        ats_filerecord=ats_filerecord,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )
    if True:
        ats_filerecord = name + ".ats"
        atsperiod = [
            (0, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (7, dt0, dtmin, dtmax, dtadj, dtfailadj),
        ]
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
    nouter, ninner = 10, 2
    hclose, rclose, relax = 1e-6, 1e-6, 0.97
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
        0: [[(0, 0, ncol - 1), -10.0]],
    }
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=welspdict,
        save_flows=False,
    )

    # ghb files
    ghbspdict = {
        0: [[(0, 0, 0), 50.0, 1.0]],
    }
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=ghbspdict,
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

    obs_lst = []
    obs_lst.append(["obs1", "head", (0, 0, 0)])
    obs_lst.append(["obs2", "head", (0, 0, 1)])
    obs_dict = {f"{gwfname}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    return sim, None


def check_output(idx, test):
    # This will fail if budget numbers cannot be read
    fpth = os.path.join(test.workspace, f"{test.name}.lst")
    mflist = flopy.utils.Mf6ListBudget(fpth)
    names = mflist.get_record_names()
    inc = mflist.get_incremental()
    msg = f"budget times not monotically increasing {inc['totim']}."
    assert np.all(np.diff(inc["totim"]) > dtmin), msg
    v = inc["totim"][-1]
    assert v == 10.0, f"Last time should be 10.  Found {v}"

    # ensure obs results changing monotonically
    fpth = os.path.join(test.workspace, test.name + ".obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    msg = f"obs times not monotically increasing {tc['time']}."
    assert np.all(np.diff(tc["time"]) > dtmin), msg
    for obsname in ["OBS1", "OBS2"]:
        v = tc[obsname]
        msg = f"{obsname} not monotically decreasing: {v}."
        assert np.all(np.diff(v) < 0), msg
    v = tc["time"][-1]
    assert v == 10.0, f"Last time should be 10.  Found {v}"


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
