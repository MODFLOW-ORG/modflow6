"""
Test to compare MODFLOW 6 groundwater transport simulation results to MT3DMS
results.  This test was first documented in Zheng and Wang (1999) (MT3DMS:
A Modular Three-Dimensional Multispecies Transport Model for Simulation of
Advection, Dispersion, and Chemical Reactions of Contaminants in Groundwater
Systems; Documentation and User's Guide) on page 130.  This is a 1D set of 4
test problems that apply incrementally varying combinations of advection,
dispersion, and reaction (sorption and decay):

* Case 1a: Advection only
* Case 1b: Advection and dispersion
* Case 1c: Advection, dispersion, and sorption
* Case 1d: Advection, dispersion, sorption, and decay
* Case 1e: Advection, dispersion, sorption, decay, immobile domain
* Case 1f: Advection, dispersion, sorption, decay, immobile domain (do not
           specify decay_sorbed in mst input file so that mf6 assumes that
           decay_sorbed = decay_aqueous.  Results should be same as Case 1e.
* Case 1g: Advection and zero-order growth
"""

import os
from pathlib import Path

import flopy
import numpy as np
from conftest import try_get_target

testgroup = "mt3dms_p01"


def p01mt3d(
    model_ws,
    al,
    retardation,
    rc1,
    mixelm,
    zeta=None,
    prsity2=None,
    rc2=None,
    zero_order_decay=False,
    mf2005s="mf2005s",
    mt3dms="mt3dms",
):
    nlay = 1
    nrow = 1
    ncol = 101
    delr = 10.0
    delc = 1.0
    delv = 1.0
    top = 0.0
    botm = [top - delv]
    Lx = (ncol - 1) * delr
    v = 0.24
    prsity = 0.25
    q = v * prsity

    perlen = 2000.0
    dt0 = perlen / 10.0
    hk = 1.0
    laytyp = 1
    rhob = 0.25
    kd = (retardation - 1.0) * prsity / rhob

    modelname_mf = "p01_mf"
    mf = flopy.modflow.Modflow(
        modelname=modelname_mf, model_ws=model_ws, exe_name=mf2005s
    )
    dis = flopy.modflow.ModflowDis(
        mf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        perlen=perlen,
    )
    ibound = np.ones((nlay, nrow, ncol), dtype=int)
    ibound[0, 0, 0] = -1
    ibound[0, 0, -1] = -1
    strt = np.zeros((nlay, nrow, ncol), dtype=float)
    h1 = q * Lx
    strt[0, 0, 0] = h1
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk, laytyp=laytyp)
    pcg = flopy.modflow.ModflowPcg(mf)
    lmt = flopy.modflow.ModflowLmt(mf)
    mf.write_input()
    mf.run_model(silent=True)

    modelname_mt = "p01_mt"
    mt = flopy.mt3d.Mt3dms(
        modelname=modelname_mt,
        model_ws=model_ws,
        exe_name=mt3dms,
        modflowmodel=mf,
    )
    c0 = 1.0
    icbund = np.ones((nlay, nrow, ncol), dtype=int)
    icbund[0, 0, 0] = -1
    sconc = np.zeros((nlay, nrow, ncol), dtype=float)
    sconc[0, 0, 0] = c0
    btn = flopy.mt3d.Mt3dBtn(
        mt,
        laycon=laytyp,
        icbund=icbund,
        prsity=prsity,
        sconc=sconc,
        dt0=dt0,
        ifmtcn=1,
    )
    dceps = 1.0e-5
    nplane = 1
    npl = 0
    nph = 4
    npmin = 0
    npmax = 8
    nlsink = nplane
    npsink = nph
    adv = flopy.mt3d.Mt3dAdv(
        mt,
        mixelm=mixelm,
        dceps=dceps,
        nplane=nplane,
        npl=npl,
        nph=nph,
        npmin=npmin,
        npmax=npmax,
        nlsink=nlsink,
        npsink=npsink,
        percel=0.5,
    )
    dsp = flopy.mt3d.Mt3dDsp(mt, al=al)
    isothm = 1
    if zeta is not None:
        isothm = 6
    ireact = 1
    if zero_order_decay:
        ireact = 100
    if rc2 is None:
        rc2 = rc1
    rct = flopy.mt3d.Mt3dRct(
        mt,
        isothm=isothm,
        ireact=ireact,
        igetsc=0,
        rhob=rhob,
        sp1=kd,
        sp2=zeta,
        prsity2=prsity2,
        rc1=rc1,
        rc2=rc2,
    )
    ssm = flopy.mt3d.Mt3dSsm(mt)
    gcg = flopy.mt3d.Mt3dGcg(mt, mxiter=10)
    mt.write_input()
    fname = os.path.join(model_ws, "MT3D001.UCN")
    if os.path.isfile(fname):
        os.remove(fname)
    mt.run_model(silent=True)

    fname = os.path.join(model_ws, "MT3D001.UCN")
    ucnobj = flopy.utils.UcnFile(fname)
    times = ucnobj.get_times()
    conc = ucnobj.get_alldata()

    fname = os.path.join(model_ws, "MT3D001.OBS")
    if os.path.isfile(fname):
        cvt = mt.load_obs(fname)
    else:
        cvt = None

    fname = os.path.join(model_ws, "MT3D001.MAS")
    mvt = mt.load_mas(fname)

    return mf, mt, conc, cvt, mvt


def p01mf6(
    model_ws,
    al,
    retardation,
    decay_rate,
    mixelm,
    zeta=None,
    prsity2=None,
    onelambda=False,
    zero_order_decay=False,
    exe="mf6",
):
    name = "p01"
    nlay, nrow, ncol = 1, 1, 101
    nper = 1
    perlen = [2000.0]
    nstp = [10]
    tsmult = [1.0]
    steady = [True]
    delr = 10.0
    delc = 1.0
    delv = 1.0
    top = 0.0
    botm = [top - delv]
    strt = 1.0
    hk = 1.0
    laytyp = 1

    Lx = (ncol - 1) * delr
    v = 0.24
    prsity = 0.25
    q = v * prsity

    rhob = 0.25
    kd = (retardation - 1.0) * prsity / rhob

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    ws = model_ws
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)
    from flopy.mf6.mfbase import VerbosityLevel

    sim.simulation_data.verbosity_level = VerbosityLevel.quiet
    sim.name_file.memory_print_option = "all"

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
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
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    strt = np.zeros((nlay, nrow, ncol), dtype=float)
    h1 = q * Lx
    strt[0, 0, 0] = h1
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=laytyp,
        k=hk,
        k33=hk,
        save_specific_discharge=True,
    )

    # chd files
    chdspd = [[(0, 0, 0), h1], [(0, 0, ncol - 1), 0.0]]
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(
        sim,
        modelname=gwtname,
        save_flows=True,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
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
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{gwtname}.adv")

    # dispersion
    dsp = flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, alh=al, ath1=0.1)

    # mass storage and transfer
    if onelambda:
        # assign sorbed decay to decay rate
        decay_rate_sorbed = decay_rate
    else:
        decay_rate_sorbed = decay_rate

    porosity_mobile = prsity
    porosity_immobile = None
    if prsity2 is not None:
        # immobile domain is active
        volfrac_immobile = 0.5
        volfrac_mobile = 1.0 - volfrac_immobile
        theta_immobile = prsity2
        porosity_immobile = theta_immobile / volfrac_immobile
        porosity_mobile = prsity / volfrac_mobile

    first_order_decay = True
    if zero_order_decay:
        first_order_decay = False
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity_mobile,
        first_order_decay=first_order_decay,
        zero_order_decay=zero_order_decay,
        decay=decay_rate,
        decay_sorbed=decay_rate_sorbed,
        sorption="linear",
        bulk_density=rhob,
        distcoef=kd,
    )

    # constant concentration
    c0 = 1.0
    cncspd = [[(0, 0, 0), c0]]
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        maxbound=len(cncspd),
        stress_period_data=cncspd,
        save_flows=False,
        pname="CNC-1",
    )

    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[[]], filename=f"{gwtname}.ssm")

    if zeta is not None:
        ist = flopy.mf6.ModflowGwtist(
            gwt,
            sorption="LINEAR",
            first_order_decay=first_order_decay,
            zero_order_decay=zero_order_decay,
            bulk_density=rhob,
            distcoef=kd,
            decay=decay_rate,
            decay_sorbed=decay_rate_sorbed,
            zetaim=zeta,
            porosity=porosity_immobile,
            volfrac=volfrac_immobile,
            filename=f"{gwtname}.ist",
            pname="IST-1",
        )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.bud",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    sim.write_simulation()
    fname = os.path.join(model_ws, gwtname + ".ucn")
    if os.path.isfile(fname):
        os.remove(fname)
    success, buff = sim.run_simulation(silent=True, report=True)
    if not success:
        print(buff)

    # load concentrations
    fname = os.path.join(model_ws, gwtname + ".ucn")
    ucnobj = flopy.utils.HeadFile(fname, precision="double", text="CONCENTRATION")
    times = ucnobj.get_times()
    conc = ucnobj.get_alldata()

    return sim, conc


def get_binaries(targets) -> tuple[Path, Path, Path]:
    return (
        targets["mf6"],
        try_get_target(targets, "mf2005s"),
        try_get_target(targets, "mt3dms"),
    )


def test_mt3dmsp01a(function_tmpdir, targets):
    longitudinal_dispersivity = 0.0
    retardation = 1.0
    decay_rate = 0.00
    mixelm = 0
    zeta = None
    prsity2 = None

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}a"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-4), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )

    # load transport budget
    # budget text:
    #     STORAGE-AQUEOUS, DECAY-AQUEOUS, STORAGE-SORBED,
    #     DECAY-SORBED, FLOW-JA-FACE, SOURCE-SINK MIX, CONSTANT CONC
    gwtname = "gwt_p01"
    fname = os.path.join(mf6_ws, f"{gwtname}.bud")
    try:
        bobj = flopy.utils.CellBudgetFile(fname, precision="double")
        budra = bobj.get_data(kstpkper=(9, 0), text="DECAY-AQUEOUS")[0]
    except:
        assert False, f'could not load data from "{fname}"'

    # ensure decay aqueous is zero
    decay_aqueous = bobj.get_data(kstpkper=(9, 0), text="DECAY-AQUEOUS")[0]
    assert np.allclose(0.0, decay_aqueous)

    # ensure decay sorbed is zero
    decay_sorbed = bobj.get_data(kstpkper=(9, 0), text="DECAY-SORBED")[0]
    assert np.allclose(0.0, decay_sorbed)

    # ensure storage sorbed is zero
    storage_sorbed = bobj.get_data(kstpkper=(9, 0), text="STORAGE-SORBED")[0]
    bobj.file.close()
    assert np.allclose(0.0, storage_sorbed), f"{storage_sorbed}"


def test_mt3dmsp01b(function_tmpdir, targets):
    longitudinal_dispersivity = 10.0
    retardation = 1.0
    decay_rate = 0.00
    mixelm = 0
    zeta = None
    prsity2 = None

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}b"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-4), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )


def test_mt3dmsp01c(function_tmpdir, targets):
    longitudinal_dispersivity = 10.0
    retardation = 1.5
    decay_rate = 0.00
    mixelm = 0
    zeta = None
    prsity2 = None

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}c"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-4), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )


def test_mt3dmsp01d(function_tmpdir, targets):
    longitudinal_dispersivity = 10.0
    retardation = 1.5
    decay_rate = 0.002
    mixelm = 0
    zeta = None
    prsity2 = None

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}d"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-4), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )


def test_mt3dmsp01e(function_tmpdir, targets):
    longitudinal_dispersivity = 10.0
    retardation = 1.5
    decay_rate = 0.002
    mixelm = 0
    zeta = 0.1
    prsity2 = 0.05

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}e"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-1), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )


def test_mt3dmsp01f(function_tmpdir, targets):
    longitudinal_dispersivity = 10.0
    retardation = 1.5
    decay_rate = 0.002
    mixelm = 0
    zeta = 0.1
    prsity2 = 0.05

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}f"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        onelambda=True,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1e-1), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )


def test_mt3dmsp01g(function_tmpdir, targets):
    longitudinal_dispersivity = 0.0
    retardation = 1.0
    decay_rate = -1.0
    mixelm = 0
    zeta = None
    prsity2 = None

    mf6, mf2005, mt3dms = get_binaries(targets)
    mf6_ws = function_tmpdir / f"{testgroup}g"
    mt3d_ws = mf6_ws / "mt3d"

    sim, conc_mf6 = p01mf6(
        mf6_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        zero_order_decay=True,
        exe=mf6,
    )

    mf, mt, conc_mt3d, cvt, mvt = p01mt3d(
        mt3d_ws,
        longitudinal_dispersivity,
        retardation,
        decay_rate,
        mixelm,
        zeta,
        prsity2,
        rc2=0.0,
        zero_order_decay=True,
        mf2005s=mf2005,
        mt3dms=mt3dms,
    )

    assert np.allclose(conc_mt3d, conc_mf6, atol=1.0e-4), (
        f"concentrations not equal {conc_mt3d} {conc_mf6}"
    )
