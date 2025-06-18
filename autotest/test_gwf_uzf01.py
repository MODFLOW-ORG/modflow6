"""
Test the ability of a uzf to route waves through a simple 1d vertical
column.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_uzf01a"]
nlay, nrow, ncol = 100, 1, 1

crs = "EPSG:26916"


def build_models(idx, test):
    name = cases[idx]

    perlen = [500.0]
    nper = len(perlen)
    nstp = [10]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    delv = 1.0
    top = 100.0
    botm = [top - (k + 1) * delv for k in range(nlay)]
    strt = 0.5
    hk = 1.0
    laytyp = 1
    ss = 0.0
    sy = 0.1

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

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
    )

    # create gwf model
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions=newtonoptions,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        crs=crs,
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
        0: [[(nlay - 1, 0, 0), 1.5, 1.0]],
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
        f"{name}.uzf.obs.csv": [
            ("wc 02", "water-content", 2, 0.5),
            ("wc 50", "water-content", 50, 0.5),
            ("wcbn 02", "water-content", "uzf 002", 0.5),
            ("wcbn 50", "water-content", "UZF 050", 0.5),
            ("rch 02", "uzf-gwrch", "uzf 002"),
            ("rch 50", "uzf-gwrch", "uzf 050"),
        ]
    }

    sd = 0.1
    vks = hk
    thtr = 0.05
    thti = thtr
    thts = sy
    eps = 4
    uzf_pkdat = [[0, (0, 0, 0), 1, 1, sd, vks, thtr, thts, thti, eps, "uzf 001"]] + [
        [k, (k, 0, 0), 0, k + 1, sd, vks, thtr, thts, thti, eps, f"uzf {k + 1:03d}"]
        for k in range(1, nlay - 1)
    ]
    uzf_pkdat[-1][3] = -1
    infiltration = 2.01
    uzf_spd = {0: [[0, infiltration, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]}
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=True,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{name}.uzf.bud",
        budgetcsv_filerecord=f"{name}.uzf.bud.csv",
        observations=uzf_obs,
        filename=f"{name}.uzf",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.bud",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    obs_lst = []
    obs_lst.append(["obs1", "head", (0, 0, 0)])
    obs_lst.append(["obs2", "head", (1, 0, 0)])
    obs_dict = {f"{name}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    return sim, None


def check_output(idx, test):
    name = test.name
    ws = test.workspace

    # check binary grid file
    fname = os.path.join(ws, name + ".dis.grb")
    grbobj = flopy.mf6.utils.MfGrdFile(fname)
    ia = grbobj._datadict["IA"] - 1
    ja = grbobj._datadict["JA"] - 1
    grb_crs = grbobj._datadict["CRS"]

    # verify crs data string in grb version 2 file
    assert grb_crs == crs

    upth = os.path.join(ws, name + ".uzf.bud")
    uobj = flopy.utils.CellBudgetFile(upth, precision="double")
    gwf_recharge = uobj.get_data(text="GWF")

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

    # Open the uzf observation file
    fpth = os.path.join(ws, f"{name}.uzf.obs.csv")
    obs_obj = flopy.utils.Mf6Obs(fpth)
    names = obs_obj.get_obsnames()
    obs = {
        names[-2]: obs_obj.get_data(obsname=names[-2]),
        names[-1]: obs_obj.get_data(obsname=names[-1]),
    }
    cbc = uobj.get_ts(idx=[[0, 0, 1], [0, 0, 49]], text="GWF")
    for i, key in enumerate(obs.keys()):
        assert np.allclose(obs[key][key], -cbc[:, i + 1]), (
            f"observation data for {key} is not the same as "
            "data in the cell-by-cell file."
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
