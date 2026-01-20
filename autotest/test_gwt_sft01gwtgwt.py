"""
Based on sft01, but split into two gwf models and two gwt models
in order to test gwf-gwf and gwt-gwt.  There are sfr and sft models
for flow and transport.  The sfr flows and the sft concentrations
should match exactly with the gwf flows and the gwf concentrations.

      flow1                        flow2
 sfr  1 2 3 4 5 6 7  gwfgwf-mvr => 1 2 3 4 5 6 7
      -------------                -------------   (sfr leakance is zero so no flow between sfr and gwf)
 gwf  1 2 3 4 5 6 7  gwfgwf     => 1 2 3 4 5 6 7
          |                            |
 gwfgwt (flow1-transport1)    gwfgwt (flow2-transport2)
          |                            |
      transport1                   transport2
 sft  1 2 3 4 5 6 7  gwtgwt-mvt => 1 2 3 4 5 6 7
      -------------                -------------
 gwt  1 2 3 4 5 6 7  gwtgwt     => 1 2 3 4 5 6 7
"""  # noqa

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sft01gwtgwt"]

# properties for each model combination
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
Kh = 20.0
Kv = 20.0

# names for the models
gwfnames = ["flow1", "flow2"]
gwtnames = ["transport1", "transport2"]

# set functionality for flow
sto_on = False
rch_on = True
sfr_on = True
lake_on = True
within_model_mvr_on = True
across_model_mvr_on = True

# set functionality for transport
transport_on = True
lkt_on = True and lake_on
sft_on = True and sfr_on
within_model_mvt_on = True and within_model_mvr_on
across_model_mvt_on = True and across_model_mvr_on


def build_models(idx, test):
    name = "mf6sim"
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
        continue_=False,
    )

    perlen = [0.1]
    nstp = [10]
    tsmult = [1.0]
    tdis_rc = []
    nper = len(perlen)
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
        filename="mfsim.tdis",
    )

    # Flow solver
    hclose = 1.0e-8
    rclose = 1.0e-6
    nouter = 700
    ninner = 300
    relax = 0.97
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="flow.ims",
    )

    # Transport solver
    imsgwt = None
    if transport_on:
        hclose = 0.001
        rclose = 0.1
        nouter = 50
        ninner = 20
        relax = 0.97
        imsgwt = flopy.mf6.ModflowIms(
            sim,
            print_option="ALL",
            outer_dvclose=hclose,
            outer_maximum=nouter,
            inner_maximum=ninner,
            under_relaxation="DBD",
            under_relaxation_theta=0.7,
            inner_dvclose=hclose,
            rcloserecord=rclose,
            linear_acceleration="BICGSTAB",
            scaling_method="NONE",
            reordering_method="NONE",
            relaxation_factor=relax,
            filename="transport.ims",
        )

    # build models and add to simulation
    for icombo in [1, 2]:
        gwfname = gwfnames[icombo - 1]
        gwtname = gwtnames[icombo - 1]
        build_gwfgwt_combo(sim, gwfname, gwtname, icombo)

    # assign models to solutions
    sim.register_ims_package(imsgwf, gwfnames)
    sim.register_ims_package(imsgwt, gwtnames)

    # add a gwf-gwf exchange
    gwfgwf_data = [
        ((0, 0, ncol - 1), (0, 0, 0), 1, delr / 2.0, delr / 2.0, delc, 0.0, delr)
    ]

    # GWF-GWF
    mvr_filerecord = None
    if across_model_mvr_on:
        mvr_filerecord = f"{gwfnames[0]}-{gwfnames[1]}.exg.mvr"
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwfnames[0],  # north
        exgmnameb=gwfnames[1],  # south
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvr_filerecord=mvr_filerecord,
        dev_interfacemodel_on=False,
        filename=f"{gwfnames[0]}-{gwfnames[1]}.exg",
    )

    # simulation GWF-GWF Mover
    if across_model_mvr_on:
        maxmvr, maxpackages = 1, 2
        mvrpack_sim = [["flow1", "sfr-1"], ["flow2", "sfr-1"]]
        mvrspd = [["flow1", "sfr-1", ncol - 1, "flow2", "sfr-1", 0, "FACTOR", 1.00]]
        gwfgwf.mvr.initialize(
            modelnames=True,
            maxmvr=maxmvr,
            print_flows=True,
            maxpackages=maxpackages,
            packages=mvrpack_sim,
            perioddata=mvrspd,
            filename=mvr_filerecord,
        )

    # GWT-GWT
    mvt_filerecord = None
    if across_model_mvt_on:
        mvt_filerecord = f"{gwtnames[0]}-{gwtnames[1]}.exg.mvr"
    gwtgwt = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwtnames[0],  # north
        exgmnameb=gwtnames[1],  # south
        gwfmodelname1=gwfnames[0],
        gwfmodelname2=gwfnames[1],
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvt_filerecord=mvt_filerecord,
        # dev_interfacemodel_on=False,
        filename=f"{gwtnames[0]}-{gwtnames[1]}.exg",
    )

    # simulation GWT-GWT Mover
    if across_model_mvt_on:
        gwtgwt.mvt.initialize(filename=mvt_filerecord)

    regression = None
    return sim, regression


def build_gwfgwt_combo(sim, gwfname, gwtname, icombo):
    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
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

    if icombo == 2:
        # add chd to right edge
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
    if icombo == 1:
        # inject water into left edge
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

    p_data = None
    if icombo == 1:
        p_data = [
            (0, "INFLOW", 1.0),
        ]

    # note: for specifying sfr number
    fname = gwfname + ".sfr.obs.csv"
    sfr_obs = {fname: [("sfrstage1", "stage", (0,)), ("sfrgwf1", "sfr", (0,))]}

    sfr = flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        mover=True,
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
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=ncol * [0.0])

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # storage
    porosity = 1.0
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = []
    if icombo == 1:
        sourcerecarray.append(("WEL-1", "AUX", "CONCENTRATION"))
    if icombo == 2:
        sourcerecarray.append(("CHD-1", "AUX", "CONCENTRATION"))
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # cnc files
    if icombo == 1:
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

    sftperioddata = None
    if icombo == 1:
        sftperioddata = [
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", 100.0),
        ]

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
        filename=f"{gwfname}-{gwtname}.exg",
    )

    return sim, None


def check_output(idx, test):
    # load the simulations
    ws = test.workspace
    test = flopy.mf6.MFSimulation.load(sim_ws=ws)

    # construct head and conc for combined models
    gwf1 = test.gwf[0]
    gwf2 = test.gwf[1]
    gwt1 = test.gwt[0]
    gwt2 = test.gwt[1]
    head = list(gwf1.output.head().get_data().flatten()) + list(
        gwf2.output.head().get_data().flatten()
    )
    conc = list(gwt1.output.concentration().get_data().flatten()) + list(
        gwt2.output.concentration().get_data().flatten()
    )

    # construct concentration in sfr for combined models
    sfrconc = []
    sfrconc += list(gwt1.sft.output.concentration().get_data().flatten())
    sfrconc += list(gwt2.sft.output.concentration().get_data().flatten())

    print("cell/sfr  head  conc  sfrconc")
    for i in range(ncol * 2):
        print(f"{i + 1}  {head[i]}  {conc[i]}  {sfrconc[i]}")
    assert np.allclose(conc, sfrconc), (
        "aquifer concentration does not equal sfr concentration"
    )


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
