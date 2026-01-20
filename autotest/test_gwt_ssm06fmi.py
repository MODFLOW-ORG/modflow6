"""
Derived from test_gwt_ssm01fmi.py, but drops RIV and adds SFR.
In that autotest, flow and transport run separately, which is the case here
as well. However, by adding SFR we can now invoke MVR, in this case
WEL -> SFR at the same time as invoking auxiliary variables.  A
companion autotest runs this same simulation (uses SFR in place of RIV for
receiving MVR water), but runs both flow and transport simultaneously which,
at one time, wrongly threw an input error.  However, running the models
separately never threw the error.
"""

import os

import flopy
import numpy as np

testgroup = "ssm06fmi"

nlay = 1
nrow = 10
ncol = 10
delr = 10.0
delc = 10.0
top = 100.0
botm = 0.0

# Add SFR for serving as a MVR receiver (something's up when multiple packages
# appear in SSM and MVR is active. When MVR is inactive, all seem to work well.
# However, things break as soon as MVR is activated.

conns = [(0, -1), (1, 0, -2), (2, 1, -3), (3, 2, -4), (4, 3)]

sfrcells = [(0, 4, 5), (0, 4, 6), (0, 4, 7), (0, 4, 8), (0, 4, 9)]

rlen = [100.0, 100.0, 100.0, 100.0, 100.0]

rbt = [99.409676, 99.320812, 99.221775, 99.146317, 99.074970]

rgrd = 0.12e-03
rwid = 20
rbth = 1.5
rbhk = 0.1
man = 0.04
ustrf = 1.0
ndv = 0


def run_flow_model(dir, exe):
    global idomain
    name = "flow"
    gwfname = name
    wsf = os.path.join(dir, testgroup, name)
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=wsf, exe_name=exe)
    tdis_rc = [(100.0, 1, 1.0), (100.0, 1, 1.0)]
    nper = len(tdis_rc)
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    # ims
    hclose = 1.0e-6
    rclose = 1.0e-6
    nouter = 1000
    ninner = 100
    relax = 0.99
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
        filename=f"{gwfname}.ims",
    )

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

    ic = flopy.mf6.ModflowGwfic(gwf, strt=100.0)

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=[1],
        k=10.0,
    )

    sto_on = False
    if sto_on:
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            save_flows=True,
            iconvert=[1],
            ss=1.0e-5,
            sy=0.3,
            steady_state={0: True},
            transient={0: False},
        )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    rch_on = False
    if rch_on:
        rch = flopy.mf6.ModflowGwfrcha(gwf, recharge={0: 4.79e-3}, pname="RCH-1")

    # wel
    wellist = []
    for i in np.arange(3, 7, 2):
        wellist.append(((0, i, 2), -100.0, 0.0))
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        mover=True,
        stress_period_data=wellist,
        auxiliary=["concentration"],
        pname="WEL-1",
    )

    # ghb
    rows = [0, 1, 2, 3]
    for ipak, i in enumerate(rows):
        blist = []
        blist.append(((0, i, ncol - 1), 50.0, 1000.0, 100.0))
        fname = f"flow.{ipak + 1}.ghb"
        ghb = flopy.mf6.ModflowGwfghb(
            gwf,
            stress_period_data=blist,
            auxiliary=["concentration"],
            filename=fname,
            pname=f"GHB-{ipak + 1}",
        )

    # drn
    rows = [7, 8, 9]
    for ipak, i in enumerate(rows):
        blist = []
        blist.append(((0, i, ncol - 1), 50.0, 1000.0, 100.0))
        fname = f"flow.{ipak + 1}.drn"
        drn = flopy.mf6.ModflowGwfdrn(
            gwf,
            stress_period_data=blist,
            auxiliary=["concentration"],
            filename=fname,
            pname=f"DRN-{ipak + 1}",
        )

    # sfr - stream starts in the middle of domain and goes due east
    sfr_pkdat = []
    for i in np.arange(len(rlen)):
        ncon = len(conns[i]) - 1
        sfr_pkdat.append(
            (
                i,
                sfrcells[i],
                rlen[i],
                rwid,
                rgrd,
                rbt[i],
                rbth,
                rbhk,
                man,
                ncon,
                ustrf,
                ndv,
            )
        )

    sfrspd = {0: [[0, "INFLOW", 86400.0]]}
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        mover=True,
        stage_filerecord=gwfname + ".sfr.stg",
        budget_filerecord=gwfname + ".sfr.bud",
        save_flows=True,
        pname="SFR-1",
        length_conversion=3.28084,
        time_conversion=86400.0,
        boundnames=False,
        nreaches=len(conns),
        packagedata=sfr_pkdat,
        connectiondata=conns,
        perioddata=sfrspd,
        filename=f"{gwfname}.sfr",
    )

    # mvr
    mvrpack = [["SFR-1"], ["WEL-1"]]
    static_mvrperioddata = []
    wel_idx = 0
    for wl in np.arange(2):  # There are only a maximum of 2 wells
        static_mvrperioddata.append(("WEL-1", wel_idx, "SFR-1", 0, "FACTOR", 1.0))
        wel_idx += 1

    mvrspd = {0: static_mvrperioddata}
    maxmvr = len(static_mvrperioddata)
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=maxmvr,
        print_flows=False,
        maxpackages=len(mvrpack),
        packages=mvrpack,
        perioddata=mvrspd,
        budget_filerecord=gwfname + ".mvr.bud",
    )

    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)
    errmsg = f"flow model did not terminate successfully\n{buff}"
    assert success, errmsg


def run_transport_model(dir, exe):
    name = "transport"
    gwtname = name
    wst = os.path.join(dir, testgroup, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=exe,
        sim_ws=wst,
        continue_=False,
    )

    tdis_rc = [(100.0, 10, 1.0), (100.0, 10, 1.0)]
    nper = len(tdis_rc)
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    # ims
    hclose = 0.001
    rclose = 0.001
    nouter = 50
    ninner = 20
    relax = 0.97
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
    )
    ic = flopy.mf6.ModflowGwtic(gwt, strt=100.0)
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=20.0, ath1=2, atv=0.2)

    # Create the ssm sources block information
    sourcerecarray = []
    sourcerecarray += [("WEL-1", "AUX", "CONCENTRATION")]
    sourcerecarray += [(f"GHB-{i + 1}", "AUX", "CONCENTRATION") for i in [0, 1, 2, 3]]
    sourcerecarray += [(f"DRN-{i + 1}", "AUX", "CONCENTRATION") for i in [0, 1, 2]]

    ssm = flopy.mf6.ModflowGwtssm(gwt, print_flows=True, sources=sourcerecarray)

    pd = [
        ("GWFHEAD", "../flow/flow.hds", None),
        ("GWFBUDGET", "../flow/flow.bud", None),
        ("GWFMOVER", "../flow/flow.mvr.bud", None),
        ("SFR-1", "../flow/flow.sfr.bud", None),
    ]
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd, flow_imbalance_correction=True)

    sftpkdat = []
    for irno in range(len(sfrcells)):
        t = (irno, 1.0)
        sftpkdat.append(t)

    sftspd = {0: [[0, "INFLOW", 1.0]]}
    sft = flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=False,
        flow_package_name="SFR-1",
        print_concentration=True,
        save_flows=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sftpkdat,
        reachperioddata=sftspd,
        pname="SFT-1",
    )

    mvt = flopy.mf6.modflow.ModflowGwtmvt(gwt, pname="MVT-1")

    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        budgetcsv_filerecord=f"{gwtname}.cbc.csv",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)
    errmsg = f"transport model did not terminate successfully\n{buff}"
    assert success, errmsg

    # ensure budget table can be parsed
    fname = gwtname + ".lst"
    fname = os.path.join(wst, fname)
    budl = flopy.utils.Mf6ListBudget(fname, budgetkey="MASS BUDGET FOR ENTIRE MODEL")
    d0 = budl.get_budget()[0]

    # Load the csv representation of the budget for confirming that the model ran
    fname = f"{gwtname}.cbc.csv"
    fname = os.path.join(wst, fname)
    d0 = np.genfromtxt(fname, names=True, delimiter=",", deletechars="")
    print(d0.dtype.names)


def test_ssm06fmi(function_tmpdir, targets):
    run_flow_model(str(function_tmpdir), targets["mf6"])
    run_transport_model(str(function_tmpdir), targets["mf6"])
