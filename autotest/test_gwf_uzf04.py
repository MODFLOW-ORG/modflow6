"""
Test uzf mass balance.  One cell model with starting water table at -20
and GHB with stage of -25.  Uzf infiltration is applied, but water table
still falls.  This test looks at the simulated unsat zone storage and
unsat volume (stored as an auxiliary variable) and compares the results
to calculated values.  Although the Uzf unsat storage and unsat volume
should probably be for total water instead of just mobile water (theta -
thetar), this is not how Uzf was designed.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_uzf04a"]
nlay, nrow, ncol = 1, 1, 1
thts = 0.30  # saturated water content
thtr = 0.05  # residual water content
thti = 0.10  # initial water content
strt = -20.0


def build_models(idx, test):
    perlen = [1.0]
    nper = len(perlen)
    nstp = [1]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    delv = 30.0
    top = 0.0
    botm = [top - (k + 1) * delv for k in range(nlay)]
    laytyp = 1
    ss = 1.0e-5
    sy = 0.3

    # unsat props
    hk = 10.0
    infiltration_rate = 0.5 * hk
    evapotranspiration_rate = 0.0
    evt_extinction_depth = 2.0
    brooks_corey_epsilon = 3.5  # brooks corey exponent

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
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = name
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions=newtonoptions,
        save_flows=True,
    )

    # create iterative model solution and register the gwf model with it
    nouter, ninner = 100, 10
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
        linear_acceleration="BICGSTAB",
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

    # ghb
    ghbspdict = {
        0: [[(nlay - 1, 0, 0), -25.0, hk / (0.5 * delv)]],
    }
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=ghbspdict,
        save_flows=False,
    )

    # note: for specifying uzf number, use fortran indexing!
    uzf_obs = {
        name + ".uzf.obs.csv": [
            (f"wc{k + 1}", "water-content", 1, depth)
            for k, depth in enumerate(np.linspace(1, 20, 15))
        ]
    }

    surfdep = 1.0e-5
    uzf_pkdat = [
        [
            0,
            (0, 0, 0),
            1,
            -1,
            surfdep,
            hk,
            thtr,
            thts,
            thti,
            brooks_corey_epsilon,
            "uzf01",
        ]
    ]
    uzf_spd = {
        0: [
            [
                0,
                infiltration_rate,
                evapotranspiration_rate,
                evt_extinction_depth,
                thtr,
                0.0,
                0.0,
                0.0,
            ],
        ]
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=True,
        simulate_et=True,
        unsat_etwc=True,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{name}.uzf.bud",
        wc_filerecord=f"{name}.uzf.bin",
        observations=uzf_obs,
        filename=f"{name}.uzf",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    obs_lst = []
    obs_lst.append(["obs1", "head", (0, 0, 0)])
    obs_dict = {f"{gwfname}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    return sim, None


def check_output(idx, test):
    name = test.name
    ws = test.workspace

    fname = os.path.join(ws, f"{name}.uzf.bin")
    wobj = flopy.utils.HeadFile(fname, text="WATER-CONTENT")
    wc = wobj.get_alldata()

    fname = os.path.join(ws, f"{name}.hds")
    wobj = flopy.utils.HeadFile(fname)
    head = wobj.get_alldata()

    bpth = os.path.join(ws, name + ".uzf.bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    qstosimobj = bobj.get_data(text="STORAGE")[0]
    qstosim = qstosimobj["q"][0]
    volume_mobile_sim = qstosimobj["VOLUME"][0]
    print("volume mobile ", volume_mobile_sim)

    # calculate volume of mobile water in unsat zone
    wc = wc.flatten()[0]
    head = head.flatten()[0]
    print("Ending Simulated water content", wc)
    print("Ending Simulated head", head)
    top = 0
    vw = (top - head) * (wc - thtr)
    v0 = (top - strt) * (thti - thtr)
    qsto = -(vw - v0) / 1.0
    print("Starting volume of mobile water in unsat zone is ", v0)
    print("Ending volume of mobile water in unsat zone is ", vw)
    print("Storage change for mobile water in unsat zone should be ", qsto)
    print("Simulated storage is ", qstosim)
    assert np.allclose(qsto, qstosim), "Simulated storage not equal known storage"
    assert np.allclose(vw, volume_mobile_sim), (
        "Simulated mobile water volume in aux does not match known result"
    )


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
