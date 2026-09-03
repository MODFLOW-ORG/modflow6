"""
LKT/LKE time-series, PERIOD AUXILIARY, and cross-period persistence
integration tests.

Case lkt_01ts — Exercises two GWT advanced package TS re-sync paths:

  PERIOD TS fields: lkt_source_period and lkt_ad_ts per-timestep re-sync.
            RAINFALL is linked to a STEPWISE TS that changes value
            mid-period (at t=0.05).  The lkt-1-rain observation verifies that
            lkt_ad_ts picks up the new value at each timestep.

  PACKAGEDATA AUX TS: bndext_ad per-timestep re-sync of pkg_auxvar
            into featureauxvar.  LKT AUX1 column carries TS name "aux1_ts"
            which also steps at t=0.05.  The LKT binary budget AUXILIARY term
            verifies that featureauxvar reflects the TS-advanced value each step.

Case lkt_01_peraux — Exercises GWT PERIOD AUXILIARY without a TS file:

  PERIOD AUXILIARY literal values: apt_rp applies PERIOD_AUXILIARY rows to
            featureauxvar even when ts_active=False, so the PACKAGEDATA
            baseline (99.0) is overridden by the PERIOD literal (55.0).

Cases lke_01ts and lke_01paux are the GWE-LKE counterparts of lkt_01ts and
lkt_01paux respectively.  They reuse the same GWF/LAK geometry but attach a
GWE model with LKE instead of GWT with LKT.

Cases lktraincont/lkeraincont — a PERIOD RAINFALL value linked to a time
series keeps tracking that time series in a later period whose own PERIOD
block reappears (for STATUS/CONCENTRATION) without repeating RAINFALL.
RAINFALL is a package-specific (non-AUX) field, wired in gwt-lkt.f90/
gwe-lke.f90 -- confirms the generalized apply_period_settings mechanism
beyond AUXILIARY.

All cases share the same GWF/LAK geometry from test_gwt_lkt01.
"""

import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = [
    "lkt_01ts",
    "lkt_01paux",
    "lke_01ts",
    "lke_01paux",
    "lktraincont",
    "lkeraincont",
]

# TS step-change time: halfway through the 0.1-day period
TS_STEP = 0.05

# RAINFALL concentration: 25.0 for steps 1-5, 0.0 for steps 6-10
RAIN_CONC_EARLY = 25.0
RAIN_CONC_LATE = 0.0

# AUX1 value for TS case: 99.0 for steps 1-5, 0.0 for steps 6-10
AUX1_EARLY = 99.0
AUX1_LATE = 0.0

# LAK volumetric rainfall rate (from GWF LAK setup)
Q_RAIN = 0.1

# Expected rainfall budget obs = Q_RAIN * RAIN_CONC
RAIN_OBS_EARLY = Q_RAIN * RAIN_CONC_EARLY  # 2.5
RAIN_OBS_LATE = Q_RAIN * RAIN_CONC_LATE  # 0.0

# PERIOD AUXILIARY case constants
PKGDATA_AUX1 = 99.0  # PACKAGEDATA literal baseline
PERIOD_AUX1 = 55.0  # PERIOD AUXILIARY literal override (should win)

# cross-period persistence: distinct value per period, 3 periods x 5 steps
TS_VAL_P1 = 10.0
TS_VAL_P2 = 20.0
TS_VAL_P3 = 30.0
QRAIN = 0.1  # constant GWF LAK rainfall flow rate; obs = QRAIN * concrain


def _build_gwf(
    sim,
    name,
    delr,
    delc,
    nlay,
    nrow,
    ncol,
    top,
    botm,
    idomain,
    hclose,
    rclose,
    relax,
    nouter,
    ninner,
):
    """Build shared GWF model (same for both cases)."""
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    flopy.mf6.ModflowIms(
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
    flopy.mf6.ModflowGwfdis(
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=20.0,
        k33=20.0,
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), -0.5, 0.0],
            [(0, 0, ncol - 1), -0.5, 0.0],
        ],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.chd",
    )
    nlakeconn = 3
    connlen = connwidth = delr / 2.0
    flopy.mf6.modflow.ModflowGwflak(
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
        packagedata=[(0, -0.4, nlakeconn, 0.0, 1025.0)],
        connectiondata=[
            (0, 0, (0, 0, 1), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", DNODATA, 10, 10, connlen, connwidth),
        ],
        outlets=[(0, 0, -1, "SPECIFIED", 999.0, 999.0, 999.0, 999.0)],
        perioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "STAGE", -0.4),
            (0, "RAINFALL", 0.1),
            (0, "EVAPORATION", 0.2),
            (0, "RUNOFF", 0.1 * delr * delc),
            (0, "WITHDRAWAL", 0.1),
            (0, "RATE", -0.1),
        ],
        pname="LAK-1",
        auxiliary=["CONCENTRATION", "DENSITY"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def _build_gwt_base(
    sim,
    name,
    delr,
    delc,
    nlay,
    nrow,
    ncol,
    top,
    botm,
    idomain,
    hclose,
    rclose,
    relax,
    nouter,
    ninner,
):
    """Build GWT model without LKT package; returns (sim, gwt, gwtname)."""
    gwtname = "gwt_" + name
    gwf_name = "gwf_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    flopy.mf6.ModflowIms(
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
    sim.register_ims_package(sim.get_package(f"{gwtname}.ims"), [gwt.name])
    flopy.mf6.ModflowGwtdis(
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
    flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30, filename=f"{gwtname}.sto")
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[("CHD-1", "AUX", "CONCENTRATION")],
        filename=f"{gwtname}.ssm",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_name,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return gwt, gwtname


def _grid_params():
    lx = 5.0
    lz = 1.0
    nlay, nrow, ncol = 1, 1, 5
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = [0.0, 0.0, -0.90, 0.0, 0.0]
    botm = list(top - np.arange(delz, nlay * delz + delz, delz))
    botm[2] = -1.0
    idomain = np.full((nlay, nrow, ncol), 1)
    return nlay, nrow, ncol, delr, delc, top, botm, idomain


def _build_gwf_default(
    sim, name, aux_name, delr, delc, nlay, nrow, ncol, top, botm, idomain
):
    """Build GWF+LAK model with flopy's default IMS settings (no custom
    solver tuning), used by the cross-period persistence cases."""
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim, model_type="gwf6", modelname=gwfname, model_nam_file=f"{gwfname}.nam"
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwfname}.ims",
    )
    flopy.mf6.ModflowGwfdis(
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=20.0, k33=20.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), -0.5, 0.0], [(0, 0, ncol - 1), -0.5, 0.0]],
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    nlakeconn = 3
    connlen = connwidth = delr / 2.0
    flopy.mf6.modflow.ModflowGwflak(
        gwf,
        save_flows=True,
        budget_filerecord="lakebud",
        nlakes=1,
        noutlets=1,
        packagedata=[(0, -0.4, nlakeconn, 0.0)],
        connectiondata=[
            (0, 0, (0, 0, 1), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", DNODATA, 10, 10, connlen, connwidth),
        ],
        outlets=[(0, 0, -1, "SPECIFIED", 999.0, 999.0, 999.0, 999.0)],
        perioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "STAGE", -0.4),
            (0, "RAINFALL", QRAIN),
        ],
        pname="LAK-1",
        auxiliary=[aux_name],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def _build_gwt_base_default(
    sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
):
    """Build GWT model without LKT, flopy IMS defaults; returns (gwt, gwtname)."""
    gwtname = "gwt_" + name
    gwfname = "gwf_" + name
    gwt = flopy.mf6.MFModel(
        sim, model_type="gwt6", modelname=gwtname, model_nam_file=f"{gwtname}.nam"
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(ims, [gwt.name])
    flopy.mf6.ModflowGwtdis(
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
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return gwt, gwtname


def _build_gwe_base_default(
    sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
):
    """Build GWE model without LKE, flopy IMS defaults; returns (gwe, gwename)."""
    gwename = "gwe_" + name
    gwfname = "gwf_" + name
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(ims, [gwename])
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=CPW,
        density_water=RHOW,
        latent_heat_vaporization=LHV,
        heat_capacity_solid=CPS,
        density_solid=RHOS,
    )
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.2700)
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "TEMPERATURE")])
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    return gwe, gwename


def _rain_perioddata(depvarkey):
    return {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, depvarkey, 0.0),
            (0, "RAINFALL", "rain_ts"),
        ],
        1: [
            (0, "STATUS", "CONSTANT"),
            (0, depvarkey, 0.0),
            # no RAINFALL setting -- period 0's TS-linked value should
            # persist and keep tracking the TS
        ],
        2: [(0, "STATUS", "CONSTANT"), (0, depvarkey, 0.0)],
    }


def build_models_lktraincont(name, ws):
    """Case lktraincont: PERIOD RAINFALL linked to a TS keeps tracking it
    in a later period whose own PERIOD block doesn't mention it, GWT-LKT."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim, name, "CONCENTRATION", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwt, gwtname = _build_gwt_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    lkt_obs = {(gwtname + ".lkt.obs.csv",): [("lkt1-rain", "rainfall", 1)]}
    lkt = flopy.mf6.ModflowGwtlkt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".lkt.bud",
        packagedata=[(0, 0.0, 999.0, "mylkt1")],
        lakeperioddata=_rain_perioddata("CONCENTRATION"),
        observations=lkt_obs,
        pname="LAK-1",
        auxiliary=["aux1"],
    )
    lkt.ts.initialize(
        filename="rain.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="rain_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_lkeraincont(name, ws):
    """Case lkeraincont: PERIOD RAINFALL continues across periods, GWE-LKE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim, name, "TEMPERATURE", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwe, gwename = _build_gwe_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    lke_obs = {(gwename + ".lke.obs.csv",): [("lke1-rain", "rainfall", 1)]}
    lke = flopy.mf6.ModflowGwelke(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".lke.bud",
        packagedata=[(0, 0.0, KTF, RBTHCND, 999.0, "mylke1")],
        lakeperioddata=_rain_perioddata("TEMPERATURE"),
        observations=lke_obs,
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="TEMPERATURE",
        pname="LAK-1",
        auxiliary=["aux1"],
    )
    lke.ts.initialize(
        filename="rain.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="rain_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_ts(name, ws):
    """Case lkt_01ts: PERIOD TS + PACKAGEDATA AUX TS."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(0.1, 10, 1.0)],
    )

    _build_gwf(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )
    gwt, gwtname = _build_gwt_base(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )

    # PACKAGEDATA: AUX1 uses TS name "aux1_ts"
    lktpackagedata = [(0, 35.0, "aux1_ts", 999.0, "mylake")]

    # PERIOD: RAINFALL uses TS name "rain_conc"
    lktperioddata = [
        (0, "STATUS", "CONSTANT"),
        (0, "CONCENTRATION", 100.0),
        (0, "RAINFALL", "rain_conc"),
        (0, "EVAPORATION", 25.0),
        (0, "RUNOFF", 25.0),
    ]

    lkt_obs = {
        (gwtname + ".lkt.obs.csv",): [
            ("lkt-1-conc", "CONCENTRATION", 1),
            ("lkt-1-rain", "RAINFALL", 1),
            ("lkt-1-evap", "EVAPORATION", 1),
            ("lkt-1-roff", "RUNOFF", 1),
        ],
    }
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
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="CONCENTRATION",
        pname="LKT-1",
        auxiliary=["aux1", "aux2"],
    )

    # STEPWISE series: rain_conc and aux1_ts both step at t=0.05
    lkt.ts.initialize(
        filename=f"{gwtname}.lkt.ts",
        timeseries=[
            (0.0, RAIN_CONC_EARLY, AUX1_EARLY),
            (TS_STEP, RAIN_CONC_LATE, AUX1_LATE),
            (0.1, RAIN_CONC_LATE, AUX1_LATE),
        ],
        time_series_namerecord=[("rain_conc", "aux1_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_peraux(name, ws):
    """Case lkt_01_peraux: PERIOD AUXILIARY literal without TS.

    PACKAGEDATA aux1=99.0; PERIOD AUXILIARY overrides to 55.0.
    Expected: AUXILIARY budget = 55.0 every step.
    """
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(0.1, 10, 1.0)],
    )

    _build_gwf(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )
    gwt, gwtname = _build_gwt_base(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )

    # PACKAGEDATA: AUX1 = 99.0 literal baseline (no TS)
    lktpackagedata = [(0, 35.0, PKGDATA_AUX1, 999.0, "mylake")]

    # PERIOD: includes literal AUXILIARY override for aux1
    # Tuple format: (ifno, "AUXILIARY", auxname, value)
    lktperioddata = [
        (0, "STATUS", "CONSTANT"),
        (0, "CONCENTRATION", 100.0),
        (0, "RAINFALL", 1.0),
        (0, "EVAPORATION", 25.0),
        (0, "RUNOFF", 25.0),
        (0, "AUXILIARY", "aux1", PERIOD_AUX1),
    ]

    lkt_obs = {
        (gwtname + ".lkt.obs.csv",): [
            ("lkt-1-conc", "CONCENTRATION", 1),
            ("lkt-1-rain", "RAINFALL", 1),
        ],
    }
    lkt_obs["digits"] = 7
    lkt_obs["print_input"] = True
    lkt_obs["filename"] = gwtname + ".lkt.obs"

    flopy.mf6.modflow.ModflowGwtlkt(
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
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="CONCENTRATION",
        pname="LKT-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


# GWE EST heat-capacity parameters (used for lke_01ts and lke_01paux)
CPW = 4183.0
RHOW = 999.728
LHV = 2500.0
CPS = 800.0
RHOS = 2650.0
KTF = 1.0
RBTHCND = 0.1


def _build_gwf_lke(
    sim,
    name,
    delr,
    delc,
    nlay,
    nrow,
    ncol,
    top,
    botm,
    idomain,
    hclose,
    rclose,
    relax,
    nouter,
    ninner,
):
    """GWF model variant for GWE coupling: CHD/LAK auxiliary = TEMPERATURE."""
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    flopy.mf6.ModflowIms(
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
    flopy.mf6.ModflowGwfdis(
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=20.0,
        k33=20.0,
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), -0.5, 0.0],
            [(0, 0, ncol - 1), -0.5, 0.0],
        ],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary="TEMPERATURE",
        filename=f"{gwfname}.chd",
    )
    nlakeconn = 3
    connlen = connwidth = delr / 2.0
    flopy.mf6.modflow.ModflowGwflak(
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
        packagedata=[(0, -0.4, nlakeconn, 0.0)],
        connectiondata=[
            (0, 0, (0, 0, 1), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", DNODATA, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", DNODATA, 10, 10, connlen, connwidth),
        ],
        outlets=[(0, 0, -1, "SPECIFIED", 999.0, 999.0, 999.0, 999.0)],
        perioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "STAGE", -0.4),
            (0, "RAINFALL", 0.1),
            (0, "EVAPORATION", 0.2),
            (0, "RUNOFF", 0.1 * delr * delc),
            (0, "WITHDRAWAL", 0.1),
            (0, "RATE", -0.1),
        ],
        pname="LAK-1",
        auxiliary=["TEMPERATURE"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def _build_gwe_base(
    sim,
    name,
    delr,
    delc,
    nlay,
    nrow,
    ncol,
    top,
    botm,
    idomain,
    hclose,
    rclose,
    relax,
    nouter,
    ninner,
):
    """Build GWE model without LKE package; returns (gwe, gwename)."""
    gwename = "gwe_" + name
    gwf_name = "gwf_" + name
    gwe = flopy.mf6.ModflowGwe(
        sim,
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    imsgwe = flopy.mf6.ModflowIms(
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
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwename])
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0, filename=f"{gwename}.ic")
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM", filename=f"{gwename}.adv")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=CPW,
        density_water=RHOW,
        latent_heat_vaporization=LHV,
        heat_capacity_solid=CPS,
        density_solid=RHOS,
        filename=f"{gwename}.est",
    )
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        ktw=0.5918,
        kts=0.2700,
        filename=f"{gwename}.cnd",
    )
    flopy.mf6.ModflowGwessm(
        gwe,
        sources=[("CHD-1", "AUX", "TEMPERATURE")],
        filename=f"{gwename}.ssm",
    )
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwf_name,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )
    return gwe, gwename


def build_models_lke_ts(name, ws):
    """Case lke_01ts: PERIOD TS + PACKAGEDATA AUX TS for GWE-LKE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(0.1, 10, 1.0)],
    )

    _build_gwf_lke(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )
    gwe, gwename = _build_gwe_base(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )

    # PACKAGEDATA: ktf, rbthcnd, AUX1 uses TS name "aux1_ts"
    lkepackagedata = [(0, 35.0, KTF, RBTHCND, "aux1_ts", 999.0, "mylake")]

    # PERIOD: RAINFALL uses TS name "rain_temp"
    lkeperioddata = [
        (0, "STATUS", "CONSTANT"),
        (0, "TEMPERATURE", 100.0),
        (0, "RAINFALL", "rain_temp"),
        (0, "EVAPORATION", 25.0),
        (0, "RUNOFF", 25.0),
    ]

    lke_obs = {
        (gwename + ".lke.obs.csv",): [
            ("lke-1-temp", "TEMPERATURE", 1),
            ("lke-1-rain", "RAINFALL", 1),
        ],
    }
    lke_obs["digits"] = 7
    lke_obs["print_input"] = True
    lke_obs["filename"] = gwename + ".lke.obs"

    lke = flopy.mf6.ModflowGwelke(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".lke.bin",
        budget_filerecord=gwename + ".lke.bud",
        packagedata=lkepackagedata,
        lakeperioddata=lkeperioddata,
        observations=lke_obs,
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="TEMPERATURE",
        pname="LKE-1",
        auxiliary=["aux1", "aux2"],
    )

    # STEPWISE series: rain_temp and aux1_ts both step at t=0.05
    lke.ts.initialize(
        filename=f"{gwename}.lke.ts",
        timeseries=[
            (0.0, RAIN_CONC_EARLY, AUX1_EARLY),
            (TS_STEP, RAIN_CONC_LATE, AUX1_LATE),
            (0.1, RAIN_CONC_LATE, AUX1_LATE),
        ],
        time_series_namerecord=[("rain_temp", "aux1_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_lke_peraux(name, ws):
    """Case lke_01paux: PERIOD AUXILIARY literal without TS for GWE-LKE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(0.1, 10, 1.0)],
    )

    _build_gwf_lke(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )
    gwe, gwename = _build_gwe_base(
        sim,
        name,
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        hclose,
        rclose,
        relax,
        nouter,
        ninner,
    )

    # PACKAGEDATA: ktf, rbthcnd, AUX1 = 99.0 literal baseline (no TS)
    lkepackagedata = [(0, 35.0, KTF, RBTHCND, PKGDATA_AUX1, 999.0, "mylake")]

    # PERIOD: includes literal AUXILIARY override for aux1
    lkeperioddata = [
        (0, "STATUS", "CONSTANT"),
        (0, "TEMPERATURE", 100.0),
        (0, "RAINFALL", 1.0),
        (0, "EVAPORATION", 25.0),
        (0, "RUNOFF", 25.0),
        (0, "AUXILIARY", "aux1", PERIOD_AUX1),
    ]

    lke_obs = {
        (gwename + ".lke.obs.csv",): [
            ("lke-1-temp", "TEMPERATURE", 1),
            ("lke-1-rain", "RAINFALL", 1),
        ],
    }
    lke_obs["digits"] = 7
    lke_obs["print_input"] = True
    lke_obs["filename"] = gwename + ".lke.obs"

    flopy.mf6.ModflowGwelke(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".lke.bin",
        budget_filerecord=gwename + ".lke.bud",
        packagedata=lkepackagedata,
        lakeperioddata=lkeperioddata,
        observations=lke_obs,
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="TEMPERATURE",
        pname="LKE-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


_BUILD_FNS = {
    "lkt_01ts": build_models_ts,
    "lkt_01paux": build_models_peraux,
    "lke_01ts": build_models_lke_ts,
    "lke_01paux": build_models_lke_peraux,
    "lktraincont": build_models_lktraincont,
    "lkeraincont": build_models_lkeraincont,
}


def build_models(idx, test):
    name = cases[idx]
    return _BUILD_FNS[name](name, test.workspace)


def check_output_ts(test):
    """Verify PERIOD TS and PACKAGEDATA AUX TS re-sync for lkt_01ts."""
    name = test.name
    gwtname = "gwt_" + name
    ws = test.workspace

    # ------------------------------------------------------------------ T1 --
    fpth = os.path.join(ws, gwtname + ".lkt.obs.csv")
    assert os.path.isfile(fpth), f"LKT obs file not found: {fpth}"
    tc = np.genfromtxt(fpth, names=True, delimiter=",")

    rain_obs = tc["LKT1RAIN"]
    assert len(rain_obs) == 10, f"Expected 10 obs, got {len(rain_obs)}"

    assert np.allclose(rain_obs[:5], RAIN_OBS_EARLY, rtol=1e-5), (
        f"Steps 1-5 rain obs expected {RAIN_OBS_EARLY}, got {rain_obs[:5]}"
    )
    assert np.allclose(rain_obs[5:], RAIN_OBS_LATE, rtol=1e-5), (
        f"Steps 6-10 rain obs expected {RAIN_OBS_LATE}, got {rain_obs[5:]}"
    )

    conc_obs = tc["LKT1CONC"]
    assert np.allclose(conc_obs, 100.0, rtol=1e-5), (
        f"Lake concentration expected 100.0 throughout, got {conc_obs}"
    )

    # ------------------------------------------------------------------ T2 --
    bud_fname = os.path.join(ws, gwtname + ".lkt.bud")
    assert os.path.isfile(bud_fname), f"LKT budget file not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")

    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records, got {len(aux_records)}"
    )

    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, (
        f"Expected AUX1 field in AUXILIARY budget record; found: {fnames}"
    )

    aux1_early = np.array([rec["AUX1"][0] for rec in aux_records[:5]])
    aux1_late = np.array([rec["AUX1"][0] for rec in aux_records[5:]])

    assert np.allclose(aux1_early, AUX1_EARLY, rtol=1e-5), (
        f"Steps 1-5 AUX1 expected {AUX1_EARLY}, got {aux1_early}"
    )
    assert np.allclose(aux1_late, AUX1_LATE, atol=1e-10), (
        f"Steps 6-10 AUX1 expected {AUX1_LATE}, got {aux1_late}"
    )


def check_output_peraux(test):
    """Verify PERIOD AUXILIARY literal (no TS) for lkt_01_peraux.

    PERIOD AUXILIARY aux1=55.0 must override the PACKAGEDATA baseline (99.0).
    The AUXILIARY budget term should show 55.0 for all 10 steps.
    """
    name = test.name
    gwtname = "gwt_" + name
    ws = test.workspace

    bud_fname = os.path.join(ws, gwtname + ".lkt.bud")
    assert os.path.isfile(bud_fname), f"LKT budget file not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")

    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records, got {len(aux_records)}"
    )

    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, (
        f"Expected AUX1 field in AUXILIARY budget record; found: {fnames}"
    )

    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    assert np.allclose(aux1_vals, PERIOD_AUX1, rtol=1e-5), (
        f"All steps AUX1 expected {PERIOD_AUX1} (PERIOD AUXILIARY), "
        f"got {aux1_vals} — PACKAGEDATA baseline {PKGDATA_AUX1} leaked through"
    )


def check_output_lke_ts(test):
    """Verify PERIOD TS and PACKAGEDATA AUX TS re-sync for lke_01ts."""
    name = test.name
    gwename = "gwe_" + name
    ws = test.workspace

    # LKE RAINFALL obs = Q * T * Cpw * rhow (energy flux, J/day)
    rain_energy_early = Q_RAIN * RAIN_CONC_EARLY * CPW * RHOW
    rain_energy_late = Q_RAIN * RAIN_CONC_LATE * CPW * RHOW  # = 0.0

    # ------------------------------------------------------------------ T1 --
    fpth = os.path.join(ws, gwename + ".lke.obs.csv")
    assert os.path.isfile(fpth), f"LKE obs file not found: {fpth}"
    tc = np.genfromtxt(fpth, names=True, delimiter=",")

    rain_obs = tc["LKE1RAIN"]
    assert len(rain_obs) == 10, f"Expected 10 obs, got {len(rain_obs)}"

    assert np.allclose(rain_obs[:5], rain_energy_early, rtol=1e-4), (
        f"Steps 1-5 rain obs expected {rain_energy_early}, got {rain_obs[:5]}"
    )
    assert np.allclose(rain_obs[5:], rain_energy_late, atol=1e-3), (
        f"Steps 6-10 rain obs expected {rain_energy_late}, got {rain_obs[5:]}"
    )

    temp_obs = tc["LKE1TEMP"]
    assert np.allclose(temp_obs, 100.0, rtol=1e-5), (
        f"Lake temperature expected 100.0 throughout, got {temp_obs}"
    )

    # ------------------------------------------------------------------ T2 --
    bud_fname = os.path.join(ws, gwename + ".lke.bud")
    assert os.path.isfile(bud_fname), f"LKE budget file not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")

    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records, got {len(aux_records)}"
    )

    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, (
        f"Expected AUX1 field in AUXILIARY budget record; found: {fnames}"
    )

    aux1_early = np.array([rec["AUX1"][0] for rec in aux_records[:5]])
    aux1_late = np.array([rec["AUX1"][0] for rec in aux_records[5:]])

    assert np.allclose(aux1_early, AUX1_EARLY, rtol=1e-5), (
        f"Steps 1-5 AUX1 expected {AUX1_EARLY}, got {aux1_early}"
    )
    assert np.allclose(aux1_late, AUX1_LATE, atol=1e-10), (
        f"Steps 6-10 AUX1 expected {AUX1_LATE}, got {aux1_late}"
    )


def check_output_lke_peraux(test):
    """Verify PERIOD AUXILIARY literal (no TS) for lke_01paux."""
    name = test.name
    gwename = "gwe_" + name
    ws = test.workspace

    bud_fname = os.path.join(ws, gwename + ".lke.bud")
    assert os.path.isfile(bud_fname), f"LKE budget file not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")

    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records, got {len(aux_records)}"
    )

    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, (
        f"Expected AUX1 field in AUXILIARY budget record; found: {fnames}"
    )

    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    assert np.allclose(aux1_vals, PERIOD_AUX1, rtol=1e-5), (
        f"All steps AUX1 expected {PERIOD_AUX1} (PERIOD AUXILIARY), "
        f"got {aux1_vals} — PACKAGEDATA baseline {PKGDATA_AUX1} leaked through"
    )


def _check_raincont(test, gwname, ext, budname, thermal):
    fname = test.workspace / f"{gwname}.{ext}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")
    rain_obs = tc["LKT1RAIN"] if "LKT1RAIN" in tc.dtype.names else tc["LKE1RAIN"]
    assert len(rain_obs) == 15, (
        f"Expected 15 obs records (5/period, 3 periods), got {len(rain_obs)}"
    )
    period1, period2, period3 = rain_obs[:5], rain_obs[5:10], rain_obs[10:]
    scale = CPW * RHOW if thermal else 1.0
    exp1 = QRAIN * TS_VAL_P1 * scale
    exp2 = QRAIN * TS_VAL_P2 * scale
    exp3 = QRAIN * TS_VAL_P3 * scale
    assert np.allclose(period1, exp1), (
        f"{budname} period 1 RAINFALL obs expected {exp1}, got {period1}"
    )
    assert np.allclose(period2, exp2), (
        f"{budname} period 2 RAINFALL obs expected {exp2} (TS should still "
        f"be tracked even though period 1's PERIOD block doesn't mention "
        f"RAINFALL), got {period2}"
    )
    assert np.allclose(period3, exp3), (
        f"{budname} period 3 RAINFALL obs expected {exp3}, got {period3}"
    )


def check_output_lktraincont(test):
    _check_raincont(test, "gwt_" + test.name, "lkt", "LKT", thermal=False)


def check_output_lkeraincont(test):
    _check_raincont(test, "gwe_" + test.name, "lke", "LKE", thermal=True)


_CHECK_FNS = {
    "lkt_01ts": check_output_ts,
    "lkt_01paux": check_output_peraux,
    "lke_01ts": check_output_lke_ts,
    "lke_01paux": check_output_lke_peraux,
    "lktraincont": check_output_lktraincont,
    "lkeraincont": check_output_lkeraincont,
}


def check_output(idx, test):
    _CHECK_FNS[cases[idx]](test)


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
