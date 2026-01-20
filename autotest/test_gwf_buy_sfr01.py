"""
Simple one-layer model with sfr on top.  Purpose is to test buy package in a
one-d sfr network.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["buy_sfr_01"]


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
    top = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
    gwtname = "gwt_" + name

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

    buy_on = True
    if buy_on:
        pd = [(0, 0.7, 0.0, gwtname, "CONCENTRATION")]
        buy = flopy.mf6.ModflowGwfbuy(gwf, denseref=1000.0, packagedata=pd)

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
    rhk = 0.01
    rman = 1.0
    ncon = 2
    ustrf = 1.0
    ndv = 0
    pak_data = []
    denseaux = 1000.0
    concaux = -99.0
    for irno in range(ncol):
        ncon = 2
        if irno in [0, ncol - 1]:
            ncon = 1
        cellid = (0, 0, irno)
        t = (
            irno,
            cellid,
            rlen,
            rwid,
            rgrd,
            rtp,
            rbth,
            rhk,
            rman,
            ncon,
            ustrf,
            ndv,
            denseaux,
            concaux,
        )
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
        auxiliary=["DTEMP", "CONCENTRATION"],
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
        flow_package_auxiliary_name="CONCENTRATION",
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
    # assign names
    gwtname = "gwt_" + test.name
    gwfname = "gwf_" + test.name

    # load the sft concentrations and make sure all values are correct
    fname = gwtname + ".sft.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    csftall = cobj.get_alldata()
    csft = csftall[-2].flatten()  # because it's lagged, get two time steps back

    # load the aquifer concentrations
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    cgwfall = cobj.get_alldata()
    cgwf = cgwfall[-2].flatten()

    # load the aquifer heads
    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    hobj = flopy.utils.HeadFile(fname, text="HEAD")
    headall = hobj.get_alldata()
    head = headall[-1].flatten()

    # load the sfr budget file and get sfr/gwf flows
    fname = gwfname + ".sfr.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=False)
    qsfrgwfsimall = bobj.get_data(text="GWF")
    qsfrgwfsim = qsfrgwfsimall[-1]
    qsfrgwfsim = qsfrgwfsim["q"]

    # load the sfr budget and check to make sure that concentrations are set
    # correctly from sft concentrations
    fname = gwfname + ".sfr.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=False)
    b = bobj.get_data(text="AUXILIARY")
    b = b[-1]
    b = b["CONCENTRATION"]
    errmsg = f"SFR aux conc not equal to SFT sim conc\n{b}\n{csft}"
    assert np.allclose(b, csft), errmsg

    # load the sfr stage file
    # load the aquifer concentrations and make sure all values are correct
    fname = gwfname + ".sfr.stg"
    fname = os.path.join(test.workspace, fname)
    stgobj = flopy.utils.HeadFile(fname, text="STAGE")
    stageall = stgobj.get_alldata()
    stage = stageall[-1]
    stage = stage.flatten()

    # calculate sfr/gwf flow rate and compare with modeled value
    cond = 0.1
    for n in range(6):
        hsfr = stage[n]
        hgwf = head[n]
        rhogwf = 1000.0 + 0.7 * cgwf[n]
        rhosfr = 1000.0 + 0.7 * csft[n]
        avgdense = 0.5 * rhogwf + 0.5 * rhosfr
        elevgwf = -0.5
        elevsfr = 0.0
        avgelev = 0.5 * elevgwf + 0.5 * elevsfr
        avghead = 0.5 * hgwf + 0.5 * hsfr
        qcalc = cond * (hsfr - hgwf)
        qcalc += cond * (avgdense / 1000.0 - 1.0) * (hsfr - hgwf)
        qcalc += cond * (avghead - avgelev) * (rhosfr - rhogwf) / 1000.0
        qsim = -qsfrgwfsim[n]
        # print(n, hsfr, hgwf, rhosfr, rhogwf, qcalc, qsim)
        # if not np.allclose(qcalc, qsim):
        #    print('reach {} flow {} not equal {}'.format(n, qcalc, qsim))
        assert np.allclose(qcalc, qsim), f"reach {n} flow {qcalc} not equal {qsim}"


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
