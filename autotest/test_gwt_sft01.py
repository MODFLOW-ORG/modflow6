"""
Simple one-layer model with sfr on top.  Purpose is to test transport in a
one-d sfr network.  Flows in the sfr network were specified to exactly
match flows in the aquifer.  A constant concentration boundary is
specified upgradient for both the stream and aquifer, so concentrations
in the stream should exactly equal the concentrations in the aquifer.
There is no flow between the stream and the aquifer.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sft_01"]


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
    top = 0.0
    botm = [top - (k + 1) * delz for k in range(nlay)]

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
        [(0, 0, ncol - 1), 0.0, 0.0],
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

    # wel files
    wellist1 = [
        [(0, 0, 0), 1.0, 0.0],
    ]
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wellist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.wel",
    )

    # pak_data = [<rno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #             <man> <ncon> <ustrf> <ndv> [<aux(naux)>] [<boundname>]]
    rlen = delr
    rwid = delc
    rgrd = 1.0
    rtp = 0.0
    rbth = 0.1
    rhk = 0.0  # 0.01
    rman = 1.0
    ncon = 2
    ustrf = 1.0
    ndv = 0
    pak_data = []
    for irno in range(ncol):
        ncon = 2
        if irno in [0, ncol - 1]:
            ncon = 1
        cellid = (0, 0, irno)
        t = (irno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, rman, ncon, ustrf, ndv)
        pak_data.append(t)

    con_data = []
    for irno in range(ncol):
        if irno == 0:
            t = (irno, -(irno + 1))
        elif irno == ncol - 1:
            t = (irno, irno - 1)
        else:
            t = (irno, irno - 1, -(irno + 1))
        con_data.append(t)

    p_data = [
        (0, "INFLOW", 1.0),
    ]

    # note: for specifying sfr number, use fortran indexing!
    fname = gwfname + ".sfr.obs.csv"
    sfr_obs = {fname: [("sfrstage1", "stage", 1), ("sfrgwf1", "sfr", 1)]}

    sfr = flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        stage_filerecord=gwfname + ".sfr.stg",
        budget_filerecord=gwfname + ".sfr.bud",
        nreaches=ncol,
        packagedata=pak_data,
        pname="SFR-1",
        connectiondata=con_data,
        perioddata=p_data,
        observations=sfr_obs,
        # auxiliary=['CONCENTRATION',
        #           'DENSITY'],
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
        strt=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        filename=f"{gwtname}.ic",
    )

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # storage
    porosity = 1.0
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")
    # sources
    sourcerecarray = [
        ("CHD-1", "AUX", "CONCENTRATION"),
        ("WEL-1", "AUX", "CONCENTRATION"),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # cnc files
    cnclist1 = [
        [(0, 0, 0), 100.0],
    ]
    cnc1 = flopy.mf6.ModflowGwtcnc(
        gwt,
        stress_period_data=cnclist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
    )

    sftpackagedata = []
    for irno in range(ncol):
        t = (irno, 0.0, 99.0, 999.0, f"myreach{irno + 1}")
        sftpackagedata.append(t)

    sftperioddata = [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 100.0)]

    sft_obs = {
        (gwtname + ".sft.obs.csv",): [
            (f"sft-{i + 1}-conc", "CONCENTRATION", i + 1) for i in range(7)
        ]
        + [
            ("sft-1-extinflow", "EXT-INFLOW", 1),
            ("sft-1-rain", "RAINFALL", 1),
            ("sft-1-roff", "RUNOFF", 1),
            ("sft-1-evap", "EVAPORATION", 1),
            ("sft-1-const", "CONSTANT", 1),
            ("sft-1-gwt1", "SFT", 1, 1),
            ("sft-1-gwt2", "SFT", 2, 1),
            ("sft-1-gwt3", "SFT", 3, 1),
            ("sft-1-myreach1", "SFT", "MYREACH1"),
        ],
    }
    # append additional obs attributes to obs dictionary
    sft_obs["digits"] = 7
    sft_obs["print_input"] = True
    sft_obs["filename"] = gwtname + ".sft.obs"

    sft = flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sftpackagedata,
        reachperioddata=sftperioddata,
        observations=sft_obs,
        pname="SFR-1",
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
    fname = gwtname + ".sft.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    # load the lake concentrations and make sure all values are correct
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    csft = cobj.get_data().flatten()

    # load the aquifer concentrations and make sure all values are correct
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    caq = cobj.get_data().flatten()

    assert np.allclose(csft, caq), f"{csft} {caq}"

    # sft observation results
    fpth = os.path.join(test.workspace, gwtname + ".sft.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    # compare observation concs with binary file concs
    for i in range(7):
        oname = f"SFT{i + 1}CONC"
        assert np.allclose(tc[oname][-1], csft[i]), f"{tc[oname][-1]} {csft[i]}"

    # load the sft budget file
    fname = gwtname + ".sft.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=False)

    # check the flow-ja-face terms
    res = bobj.get_data(text="flow-ja-face")[-1]
    # print(res)

    # check the storage terms, which include the total mass in the reach
    # as an aux variable
    res = bobj.get_data(text="storage")[-1]
    # print(res)

    # check the constant term
    res = bobj.get_data(text="constant")[-1]
    qs = res["q"]
    qa = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    msg = f"{qs} /= {qa}"
    assert np.allclose(qs, qa), msg


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
