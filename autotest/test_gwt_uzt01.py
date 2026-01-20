"""
Test uzt for one-d transport in a vertical column.  This problem is based
on test_gwf_uzf03.py.  Infiltration is assigned a concentration of 100.  The
uzet concentration is also assigned as 100. so that the calculated uzf cells
should have a concentration of 100.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["uzt01a"]
nlay, nrow, ncol = 15, 1, 1


def build_models(idx, test):
    perlen = [17.7]
    nper = len(perlen)
    nstp = [177]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    delv = 2.0
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
    gwfname = "gwf_" + name
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
        gwfname + ".uzf.obs.csv": [
            (f"wc{k + 1}", "water-content", k + 1, 0.5 * delv) for k in range(nlay)
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
        budget_filerecord=f"{gwfname}.uzf.bud",
        observations=uzf_obs,
        pname="UZF-1",
        filename=f"{gwfname}.uzf",
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
    obs_lst.append(["obs2", "head", (1, 0, 0)])
    obs_dict = {f"{gwfname}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(
        sim,
        save_flows=True,
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # storage
    porosity = sy
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")
    # sources
    sourcerecarray = [
        (),
        # ('WEL-1', 'AUX', 'CONCENTRATION'),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    uztpackagedata = [(iuz, 0.0, f"myuzt{iuz + 1}") for iuz in range(nlay)]
    uztperioddata = [
        (0, "INFILTRATION", 100.0),
        (0, "UZET", 100.0),
    ]

    ncv = nlay
    uzt_obs = {}
    for obstype in [
        "CONCENTRATION",
        "STORAGE",
        "CONSTANT",
        "FROM-MVR",
        "UZT",
        "INFILTRATION",
        "REJ-INF",
        "UZET",
        "REJ-INF-TO-MVR",
    ]:
        fname = f"{gwtname}.uzt.obs.{obstype.lower()}.csv"
        obs1 = [(f"uzt{i + 1}", obstype, i + 1) for i in range(ncv)]
        obs2 = [(f"buzt{i + 1}", obstype, f"myuzt{i + 1}") for i in range(ncv)]
        uzt_obs[fname] = obs1 + obs2

    obstype = "FLOW-JA-FACE"
    fname = f"{gwtname}.uzt.obs.{obstype.lower()}.csv"
    obs1 = []
    for id1 in range(ncv):
        id2list = []
        if id1 > 0:
            id2list.append(id1 - 1)
        if id1 < ncv - 1:
            id2list.append(id1 + 1)
        for id2 in id2list:
            obs1.append((f"uzt{id1 + 1}x{id2 + 1}", obstype, id1 + 1, id2 + 1))
    obs2 = [(f"buzt{i + 1}", obstype, f"myuzt{i + 1}") for i in range(ncv)]
    uzt_obs[fname] = obs1 + obs2

    # append additional obs attributes to obs dictionary
    uzt_obs["digits"] = 15
    uzt_obs["print_input"] = True
    uzt_obs["filename"] = gwtname + ".uzt.obs"

    uzt = flopy.mf6.modflow.ModflowGwtuzt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".uzt.bin",
        budget_filerecord=gwtname + ".uzt.bud",
        packagedata=uztpackagedata,
        uztperioddata=uztperioddata,
        observations=uzt_obs,
        flow_package_name="UZF-1",
        pname="UZT-1",
    )
    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.bud",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[
            ("CONCENTRATION", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("CONCENTRATION", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def make_plot(sim, obsvals):
    print("making plots...")

    name = sim.name
    ws = sim.workspace

    # shows curves for times 2.5, 7.5, 12.6, 17.7
    # which are indices 24, 74, 125, and -1
    indices = [24, 74, 125, -1]

    obsvals = [list(row) for row in obsvals]
    obsvals = [obsvals[i] for i in indices]
    obsvals = np.array(obsvals)

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1)
    depth = np.arange(1, 31, 2.0)
    for row in obsvals:
        label = f"time {row[0]}"
        ax.plot(row[1:], depth, label=label, marker="o")
    ax.set_ylim(0.0, 30.0)
    ax.set_xlim(0.0, 100.0)
    ax.invert_yaxis()
    ax.set_xlabel("Concentration")
    ax.set_ylabel("Depth, in meters")
    plt.legend()

    fname = "fig-xsect.pdf"
    fname = os.path.join(ws, fname)
    plt.savefig(fname, bbox_inches="tight")

    return


def check_obs(sim):
    print("checking obs...")
    name = sim.name
    ws = sim.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    ncv = nlay

    # extract uzt concentrations from binary output file
    conc_uzt = gwt.uzt.output.concentration().get_alldata()
    ntimes = conc_uzt.shape[0]
    conc_uzt = conc_uzt.reshape((ntimes, ncv))

    # ensure uzt obs are the same whether specified by
    # boundname or by control volume
    csvfiles = gwt.uzt.obs.output.obs_names
    for csvfile in csvfiles:
        if ".flow-ja-face.csv" in csvfile:
            continue
        print(f"Checking csv file: {csvfile}")
        conc_ra = gwt.uzt.obs.output.obs(f=csvfile).data
        success = True
        # check boundname observations with numeric ID observations
        for icv in range(ncv):
            # print(f"  Checking control volume {icv + 1}")

            if ".concentration.csv" in csvfile:
                is_same = np.allclose(conc_ra[f"BUZT{icv + 1}"], conc_uzt[:, icv])
                if not is_same:
                    success = False
                    print(
                        "Binary concentrations do not match with "
                        "observation concentrations for uzt1"
                    )
                    print(conc_ra["BUZT1"], conc_uzt)

            is_same = np.allclose(conc_ra[f"UZT{icv + 1}"], conc_ra[f"BUZT{icv + 1}"])
            if not is_same:
                success = False
                for t, x, y in zip(
                    conc_ra["totim"],
                    conc_ra[f"UZT{icv + 1}"],
                    conc_ra[f"BUZT{icv + 1}"],
                ):
                    print(t, x, y)

    # Sum individual iconn uzt rates and compare with total rate
    csvfile = f"{gwtname}.uzt.obs.flow-ja-face.csv"
    print(f"Checking csv file: {csvfile}")
    conc_ra = gwt.uzt.obs.output.obs(f=csvfile).data
    ntimes = conc_ra.shape[0]
    for iuzt in range(ncv):
        connection_sum = np.zeros(ntimes)
        for column_name in conc_ra.dtype.names:
            if f"UZT{icv + 1}X" in column_name:
                connection_sum += conc_ra[column_name]
        is_same = np.allclose(connection_sum, conc_ra[f"BUZT{icv + 1}"])
        if not is_same:
            success = False
            diff = connection_sum - conc_ra[f"BMWTUZT{icv + 1}"]
            print(
                f"Problem with UZT {icv + 1}; "
                f"mindiff {diff.min()} and maxdiff {diff.max()}"
            )

    assert success, "One or more UZT obs checks did not pass"


def check_output(idx, test):
    name = test.name
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    ws = test.workspace

    # check binary grid file
    fname = os.path.join(ws, gwfname + ".dis.grb")
    grbobj = flopy.mf6.utils.MfGrdFile(fname)
    ia = grbobj._datadict["IA"] - 1
    ja = grbobj._datadict["JA"] - 1

    cpth = os.path.join(ws, gwtname + ".uzt.bin")
    cobj = flopy.utils.HeadFile(cpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_alldata()
    for conc_this_time in conc:
        c = conc_this_time.flatten()
        errmsg = f"conc[0] must be 100 and conc[-1] must be 0: {c}"
        # assert np.allclose(c[0], 100.), errmsg
        assert np.allclose(c[-1], 0.0), errmsg

    bpth = os.path.join(ws, gwtname + ".uzt.bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    gwt_recharge = bobj.get_data(text="GWF")

    bpth = os.path.join(ws, gwtname + ".bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    flow_ja_face = bobj.get_data(text="FLOW-JA-FACE")
    uzt_recharge = bobj.get_data(text="UZT")
    errmsg = "uzt rch is not equal to negative gwt rch"
    for gwr, uzr in zip(gwt_recharge, uzt_recharge):
        assert np.allclose(gwr["q"], -uzr["q"]), errmsg

    # Check on residual, which is stored in diagonal position of
    # flow-ja-face.  Residual should be less than convergence tolerance,
    # or this means the residual term is not added correctly.
    for fjf in flow_ja_face:
        fjf = fjf.flatten()
        res = fjf[ia[:-1]]
        errmsg = f"min or max residual too large {res.min()} {res.max()}"
        assert np.allclose(res, 0.0, atol=1.0e-6), errmsg

    bpth = os.path.join(ws, gwtname + ".uzt.bud")
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    uzet = bobj.get_data(text="UZET")
    uz_answer = [-0.432] + 14 * [0.0]
    for uz in uzet[100:]:  # Need to look later in simulation when ET demand is met
        msg = f"unsat ET not correct.  Found {uz['q']}.  Should be {uz_answer}"
        assert np.allclose(uz["q"], uz_answer), msg

    uzinfil = bobj.get_data(text="INFILTRATION")
    uz_answer = [17.28] + 14 * [0.0]
    for uz in uzinfil[20:]:
        assert np.allclose(uz["q"], uz_answer), "unsat ET is not correct"

    # Check ending concentrations
    fname = os.path.join(ws, gwtname + ".uzt.bin")
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    c = cobj.get_data().flatten()
    canswer = 10 * [100.0] + 5 * [0.0]
    msg = f"Ending uzf concentrations {c} do not match known concentrations {canswer}"
    assert np.allclose(c, canswer)

    # check observations
    check_obs(test)

    # Make plot of obs
    fpth = os.path.join(test.workspace, gwtname + ".uzt.obs.concentration.csv")
    obsvals = np.genfromtxt(fpth, names=True, delimiter=",")

    # make_plot(sim, obsvals)


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
