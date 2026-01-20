"""
Simple one-layer model with a lak.  Purpose is to test outlets that
move solute from one lake to another.
"""

import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["lkt_02"]


def build_models(idx, test):
    lx = 7.0
    lz = 1.0
    nlay = 1
    nrow = 1
    ncol = 7
    nper = 1
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = [0.0, 0.0, -0.90, -0.90, -0.90, 0.0, 0.0]
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
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
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
        [(0, 0, 0), -0.4, 100.0],
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

    # pak_data = [ifno, strt, nlakeconn, CONC, dense, boundname]
    pak_data = [
        (0, -0.4, 2, 0.0, 1025.0),
        (1, -0.4, 1, 0.0, 1025.0),
        (2, -0.4, 2, 0.0, 1025.0),
    ]

    connlen = connwidth = delr / 2.0
    con_data = []
    # con_data=(ifno,iconn,(cellid),claktype,bedleak,belev,telev,connlen,connwidth )
    # lake 1
    con_data.append(
        (0, 0, (0, 0, 1), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth)
    )
    con_data.append((0, 1, (0, 0, 2), "VERTICAL", DNODATA, 10, 10, connlen, connwidth))
    # lake 2
    con_data.append((1, 0, (0, 0, 3), "VERTICAL", DNODATA, 10, 10, connlen, connwidth))
    # lake 3
    con_data.append((2, 0, (0, 0, 4), "VERTICAL", DNODATA, 10, 10, connlen, connwidth))
    con_data.append(
        (2, 1, (0, 0, 5), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth)
    )

    p_data = [
        (0, "RAINFALL", 0.0),
        (1, "RAINFALL", 0.0),
        (2, "RAINFALL", 0.0),
    ]
    # <outletno> <lakein> <lakeout> <couttype> <invert> <width> <rough> <slope>
    outlets = [
        (0, 0, 1, "SPECIFIED", 999.0, 999.0, 999.0, 999.0),
        (1, 1, 2, "SPECIFIED", 999.0, 999.0, 999.0, 999.0),
    ]
    outletperioddata = [
        (0, "RATE", -0.1),
        (1, "RATE", -0.1),
    ]

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
        stage_filerecord=gwfname + ".lak.stg",
        budget_filerecord=gwfname + ".lak.bud",
        nlakes=3,
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
    ic = flopy.mf6.ModflowGwtic(
        gwt,
        strt=[100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        filename=f"{gwtname}.ic",
    )

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
        (0, 0.0, 99.0, 999.0, "mylake1"),
        (1, 0.0, 99.0, 999.0, "mylake2"),
        (2, 0.0, 99.0, 999.0, "mylake3"),
    ]
    lktperioddata = [
        (0, "STATUS", "ACTIVE"),
        (1, "STATUS", "ACTIVE"),
        (2, "STATUS", "ACTIVE"),
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
            ("lkt-1-gwt1", "LKT", 1, 1),
            ("lkt-1-gwt2", "LKT", 1, 2),
            ("lkt-2-gwt1", "LKT", 2, 1),
            ("lkt-1-mylake1", "LKT", "MYLAKE1"),
            ("lkt-1-fjf", "FLOW-JA-FACE", 1, 2),
            ("lkt-2-fjf", "FLOW-JA-FACE", 2, 1),
            ("lkt-3-fjf", "FLOW-JA-FACE", 2, 3),
            ("lkt-4-fjf", "FLOW-JA-FACE", 3, 2),
            ("lkt-5-fjf", "FLOW-JA-FACE", "MYLAKE1"),
            ("lkt-6-fjf", "FLOW-JA-FACE", "MYLAKE2"),
            ("lkt-7-fjf", "FLOW-JA-FACE", "MYLAKE3"),
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
        budget_filerecord=gwtname + ".lkt.bud",
        packagedata=lktpackagedata,
        lakeperioddata=lktperioddata,
        observations=lkt_obs,
        pname="LAK-1",
        auxiliary=["aux1", "aux2"],
    )
    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
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


def check_output(idx, test):
    # ensure lake concentrations were saved
    name = test.name
    gwtname = "gwt_" + name
    fname = gwtname + ".lkt.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    # load the lake concentrations and make sure all values are correct
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    clak = cobj.get_data()
    answer = np.array([2.20917292e-01, 2.06602236e-03, 1.64115202e-05])
    assert np.allclose(clak, answer), f"{clak} {answer}"

    # load the aquifer concentrations and make sure all values are correct
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    caq = cobj.get_data()
    answer = np.array(
        [
            1.00000000e02,
            8.50692758e00,
            5.71603342e-01,
            1.30064549e-02,
            2.38402430e-04,
            3.30714374e-06,
            7.33450747e-08,
        ]
    )
    assert np.allclose(caq, answer), f"{caq.flatten()} {answer}"

    # lkt observation results
    fpth = os.path.join(test.workspace, gwtname + ".lkt.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    res = tc["LKT1CONC"]
    answer = [
        0.00418356,
        0.01249386,
        0.02487472,
        0.04127051,
        0.06162619,
        0.08588733,
        0.114,
        0.145911,
        0.1815675,
        0.2209173,
    ]
    answer = np.array(answer)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1STOR"]
    answer = [
        -0.198852,
        -0.3950041,
        -0.5884846,
        -0.7793216,
        -0.9675429,
        -1.153176,
        -1.336248,
        -1.516786,
        -1.694817,
        -1.870368,
    ]
    answer = np.array(answer)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1MYLAKE1"]
    answer = [
        0.1992704,
        0.3962535,
        0.5909721,
        0.7834487,
        0.9737055,
        1.161765,
        1.347648,
        1.531377,
        1.712974,
        1.89246,
    ]
    answer = np.array(answer)
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT1FJF"]
    answer = -tc["LKT2FJF"]
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT3FJF"]
    answer = -tc["LKT4FJF"]
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT5FJF"]
    answer = tc["LKT1FJF"]
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT6FJF"]
    answer = tc["LKT2FJF"] + tc["LKT3FJF"]
    assert np.allclose(res, answer), f"{res} {answer}"
    res = tc["LKT7FJF"]
    answer = tc["LKT4FJF"]
    assert np.allclose(res, answer), f"{res} {answer}"

    # load the lake budget file
    fname = gwtname + ".lkt.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=False)
    # check the flow-ja-face terms
    res = bobj.get_data(text="flow-ja-face")[-1]
    answer = [
        (1, 2, -0.02209173),
        (2, 1, 0.02209173),
        (2, 3, -0.0002066),
        (3, 2, 0.0002066),
    ]
    dt = [("node", "<i4"), ("node2", "<i4"), ("q", "<f8")]
    answer = np.array(answer, dtype=dt)
    for dtname, dttype in dt:
        assert np.allclose(res[dtname], answer[dtname]), f"{res} {answer}"
    # check the storage terms, which include the total mass in the lake
    # as an aux variable
    res = bobj.get_data(text="storage")[-1]
    answer = [
        (1, 1, -1.87036784e00, 1.05006009e-01),
        (2, 2, -2.18851269e-02, 8.85969084e-04),
        (3, 3, -2.10991316e-04, 6.88879732e-06),
    ]
    dt = [("node", "<i4"), ("node2", "<i4"), ("q", "<f8"), ("MASS", "<f8")]
    answer = np.array(answer, dtype=dt)
    for dtname, dttype in dt:
        assert np.allclose(res[dtname], answer[dtname]), f"{res} {answer}"


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
