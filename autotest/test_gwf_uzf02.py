"""
Test uzf for the vs2d comparison problem in the uzf documentation
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_uzf02a"]
nlay, nrow, ncol = 1, 1, 1


def build_models(idx, test):
    perlen = [17.7]
    nper = len(perlen)
    nstp = [177]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    delv = 30.0
    top = 0.0
    botm = [top - (k + 1) * delv for k in range(nlay)]
    strt = -22.0
    laytyp = 1
    ss = 0.0
    sy = 0.4

    # unsat props
    seconds_to_days = 60.0 * 60.0 * 24.0
    hk = 4.0e-6 * seconds_to_days  # saturated vertical conductivity
    thts = 0.4  # saturated water content
    thtr = 0.2  # residual water content
    thti = thtr  # initial water content
    infiltration_rate = 0.5 * hk
    evapotranspiration_rate = 5e-8 * seconds_to_days
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
        0: [[(nlay - 1, 0, 0), strt, hk / (0.5 * delv)]],
    }
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=ghbspdict,
        save_flows=False,
    )

    # note: for specifying lake number, use fortran indexing!
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
            1,
            surfdep,
            hk,
            thtr,
            thts,
            thti,
            brooks_corey_epsilon,
            "uzf01",
        ]
    ] + [
        [
            k,
            (k, 0, 0),
            0,
            k + 1,
            surfdep,
            hk,
            thtr,
            thts,
            thti,
            brooks_corey_epsilon,
            f"uzf0{k + 1}",
        ]
        for k in range(1, nlay)
    ]
    uzf_pkdat[-1][3] = -1
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


def make_plot(sim, obsvals):
    print("making plots...")

    # shows curves for times 2.5, 7.5, 12.6, 17.7
    # which are indices 24, 74, 125, and -1
    indices = [24, 74, 125, -1]

    obsvals = [list(row) for row in obsvals]
    obsvals = [obsvals[i] for i in indices]
    obsvals = np.array(obsvals)

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1)
    depth = np.arange(1, 20, 2.0)
    depth = np.linspace(1, 20, 15)
    for row in obsvals:
        label = f"time {row[0]}"
        ax.plot(row[1:], depth, label=label, marker="o")
    ax.set_ylim(0.0, 20.0)
    ax.set_xlim(0.15, 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("Water Content")
    ax.set_ylabel("Depth, in meters")
    plt.legend()

    fname = "fig-xsect.pdf"
    fname = os.path.join(sim.workspace, fname)
    plt.savefig(fname, bbox_inches="tight")


def check_output(idx, test):
    name = test.name
    ws = test.workspace

    # check binary grid file
    fname = os.path.join(ws, name + ".dis.grb")
    grbobj = flopy.mf6.utils.MfGrdFile(fname)
    ia = grbobj._datadict["IA"] - 1
    ja = grbobj._datadict["JA"] - 1

    bpth = os.path.join(ws, name + ".uzf.bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    gwf_recharge = bobj.get_data(text="GWF")

    bpth = os.path.join(ws, name + ".bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    flow_ja_face = bobj.get_data(text="FLOW-JA-FACE")
    uzf_recharge = bobj.get_data(text="UZF-GWRCH")
    errmsg = "uzf rch is not equal to negative gwf rch"
    for gwr, uzr in zip(gwf_recharge, uzf_recharge):
        assert np.allclose(gwr["q"], -uzr["q"]), errmsg

    # Check on residual, which is stored in diagonal position of
    # flow-ja-face.  Residual should be less than convergence tolerance,
    # or this means the residual term is not added correctly.
    for fjf in flow_ja_face:
        fjf = fjf.flatten()
        res = fjf[ia[:-1]]
        errmsg = f"min or max residual too large {res.min()} {res.max()}"
        assert np.allclose(res, 0.0, atol=1.0e-6), errmsg

    bpth = os.path.join(ws, name + ".uzf.bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    uzet = bobj.get_data(text="UZET")
    uz_answer = [-0.00432]
    for uz in uzet[20:]:  # start at 20 when et can be met
        assert np.allclose(uz["q"], uz_answer), "unsat ET is not correct"

    # Make plot of obs
    fpth = os.path.join(test.workspace, name + ".uzf.obs.csv")
    try:
        obsvals = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    # make_plot(sim, obsvals)


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
