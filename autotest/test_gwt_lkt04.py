"""
Simple one-layer model with a lak.  Purpose is to test a lake
with a variable stage and variable concentration.  The lake
starts at a concentration of 100. and slowly decreases as
fresh groundwater flows into it.  Concentrations in the aquifer
should remain at zero.
"""

import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["lkt_04"]


def build_models(idx, test):
    lx = 5.0
    lz = 1.0
    nlay = 1
    nrow = 1
    ncol = 5
    nper = 1
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = [0.0, 0.0, -0.90, 0.0, 0.0]
    botm = list(top - np.arange(delz, nlay * delz + delz, delz))
    botm[2] = -1.0

    perlen = [0.1]
    nstp = [10]
    kstp = perlen[0] / nstp[0]
    tsmult = [1.0]

    Kh = 20.0
    Kv = 20.0

    steady = [True]
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    single_matrix = False
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name

    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )

    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )

    idomain = np.full((nlay, nrow, ncol), 1)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=Kh,
        k33=Kv,
    )

    # chd files
    chdlist1 = [
        [(0, 0, 0), -0.5, 0.0],
        [(0, 0, ncol - 1), -0.5, 0.0],
    ]
    chd1 = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.chd",
    )

    # note: this is the number of connections for a lake,
    # not total number of connections
    nlakeconn = 3
    # pak_data = [ifno, strt, nlakeconn, CONC, dense, boundname]
    pak_data = [(0, -0.4, nlakeconn, 0.0, 1025.0)]

    connlen = connwidth = delr / 2.0
    con_data = []
    # con_data=(ifno,iconn,(cellid),claktype,bedleak,belev,telev,connlen,connwidth )
    con_data.append(
        (0, 0, (0, 0, 1), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth)
    )
    con_data.append(
        (0, 1, (0, 0, 3), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth)
    )
    con_data.append((0, 2, (0, 0, 2), "VERTICAL", DNODATA, 10, 10, connlen, connwidth))
    p_data = [
        (0, "STATUS", "ACTIVE"),
        (0, "STAGE", -0.4),
        (0, "RAINFALL", 0.1),
        (0, "EVAPORATION", 0.2),
        (0, "RUNOFF", 0.1 * delr * delc),
        (0, "WITHDRAWAL", 0.1),
    ]
    # <outletno> <lakein> <lakeout> <couttype> <invert> <width> <rough> <slope>
    outlets = [(0, 0, -1, "SPECIFIED", 999.0, 999.0, 999.0, 999.0)]
    outletperioddata = [(0, "RATE", -0.1)]

    # note: for specifying lake number, use fortran indexing!
    lak_obs = {
        ("lak_obs.csv"): [
            ("lakestage", "stage", 1),
            ("lakevolume", "volume", 1),
            ("lak1", "lak", 1, 1),
            ("lak2", "lak", 1, 2),
            ("lak3", "lak", 1, 3),
        ]
    }

    lak = flopy.mf6.modflow.ModflowGwflak(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        stage_filerecord="stage",
        budget_filerecord="lakebud",
        budgetcsv_filerecord=f"{gwfname}.lak.bud.csv",
        nlakes=1,
        ntables=0,
        noutlets=1,
        packagedata=pak_data,
        outlets=outlets,
        pname="LAK-1",
        connectiondata=con_data,
        perioddata=p_data + outletperioddata,
        observations=lak_obs,
        auxiliary=["CONCENTRATION", "DENSITY"],
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

    if not single_matrix:
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
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # storage
    porosity = 0.30
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")
    # sources
    sourcerecarray = [
        ("CHD-1", "AUX", "CONCENTRATION"),
        # ('WEL-1', 'AUX', 'CONCENTRATION'),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    lktpackagedata = [
        (0, 100.0, 99.0, 999.0, "mylake"),
    ]
    lktperioddata = [
        (0, "STATUS", "ACTIVE"),
        (0, "RAINFALL", 25.0),
        (0, "EVAPORATION", 0.0),
        (0, "RUNOFF", 25.0),
    ]

    lkt_obs = {
        (gwtname + ".lkt.obs.csv",): [
            ("lkt-1-conc", "CONCENTRATION", 1),
            ("lkt-1-extinflow", "EXT-INFLOW", 1),
            ("lkt-1-rain", "RAINFALL", 1),
            ("lkt-1-roff", "RUNOFF", 1),
            ("lkt-1-evap", "EVAPORATION", 1),
            ("lkt-1-wdrl", "WITHDRAWAL", 1),
            ("lkt-1-stor", "STORAGE", 1),
            ("lkt-1-const", "CONSTANT", 1),
            ("lkt-1-gwt2", "LKT", 1, 1),
            ("lkt-1-gwt4", "LKT", 1, 3),
            ("lkt-1-gwt3", "LKT", 1, 2),
            ("lkt-1-mylake", "LKT", "MYLAKE"),
        ],
    }
    # append additional obs attributes to obs dictionary
    lkt_obs["digits"] = 7
    lkt_obs["print_input"] = True
    lkt_obs["filename"] = gwtname + ".lkt.obs"

    lkt = flopy.mf6.modflow.ModflowGwtlkt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".lkt.bin",
        budget_filerecord="gwtlak1.bud",
        budgetcsv_filerecord=f"{gwtname}.lkt.bud.csv",
        packagedata=lktpackagedata,
        lakeperioddata=lktperioddata,
        observations=lkt_obs,
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="CONCENTRATION",
        pname="LKT-1",
        auxiliary=["aux1", "aux2"],
    )
    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
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


def get_mfsim(testsim):
    ws = testsim.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    return sim


def eval_csv_information(testsim):
    sim = get_mfsim(testsim)
    name = testsim.name
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    success = True
    atol = 0.002

    # Lake budget checks
    lak_budget = gwf.lak.output.budgetcsv().data
    result = lak_budget["PERCENT_DIFFERENCE"]
    for pd in result:
        if abs(pd) > atol:
            success = False
            print(f"Lake package balance error ({pd}) > tolerance ({atol})")

    # Lake transport budget checks
    lkt_budget = gwt.lkt.output.budgetcsv().data
    result = lkt_budget["PERCENT_DIFFERENCE"]
    for pd in result:
        if abs(pd) > atol:
            success = False
            print(f"Lake transport package balance error ({pd}) > tolerance ({atol})")

    assert success, "One or more errors encountered in budget checks"


def check_output(idx, test):
    # eval csv files
    eval_csv_information(test)

    # ensure lake concentrations were saved
    name = test.name
    gwtname = "gwt_" + name
    fname = gwtname + ".lkt.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    # load the lake concentrations and make sure all values are 100.
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    clak = cobj.get_alldata().flatten()
    clak_answer = np.array(
        [
            99.6180852,
            99.23811519,
            98.86008004,
            98.48396992,
            98.10977501,
            97.73748558,
            97.36709191,
            96.99858435,
            96.63195329,
            96.26718919,
        ]
    )
    assert np.allclose(clak, clak_answer), f"{clak} {clak_answer}"

    # load the aquifer concentrations and make sure all values are correct
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    caq = cobj.get_alldata()
    answer = np.zeros(5)
    assert np.allclose(caq[-1].flatten(), answer), f"{caq[-1].flatten()} {answer}"

    # lkt observation results
    fpth = os.path.join(test.workspace, gwtname + ".lkt.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    res = tc["LKT1CONC"]
    answer = clak_answer
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1EXTINFLOW"]
    answer = np.ones(10) * 0.0
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1RAIN"]
    answer = np.ones(10) * 2.5
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1ROFF"]
    answer = np.ones(10) * 2.5
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1EVAP"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1WDRL"]
    answer = clak_answer * -0.1
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1STOR"]
    answer = np.array(
        [
            14.92362,
            14.84762,
            14.77202,
            14.69679,
            14.62196,
            14.5475,
            14.47342,
            14.39972,
            14.32639,
            14.25344,
        ]
    )
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1CONST"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1GWT2"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1GWT4"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1GWT3"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1MYLAKE"]
    answer = np.zeros(10)
    assert np.allclose(res, answer), f"{res} {answer}"


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
