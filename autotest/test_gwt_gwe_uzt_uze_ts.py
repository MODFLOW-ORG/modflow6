"""
UZT/UZE time-series, PERIOD AUXILIARY, and cross-period persistence
integration tests.

Covers GWT-UZT and GWE-UZE in paired cases that share the same GWF+UZF model.

  PERIOD TS: CONCENTRATION/TEMPERATURE (CONSTANT status) and INFILTRATION
             (uzt_ad_ts/uze_ad_ts's own re-sync) linked to a STEPWISE TS
             that steps at t=0.05, applied to feature 2 (the only cell not
             CHD-adjacent, so the only one with measurable INFILTRATION).

  PACKAGEDATA AUX TS: AUX1 carries TS name "aux1_ts"; the binary budget
             AUXILIARY term verifies per-timestep re-sync.

  PERIOD AUXILIARY literal: PERIOD AUXILIARY overrides PACKAGEDATA baseline
             (99.0 -> 55.0). Binary budget AUXILIARY must show 55.0.

  Cross-period persistence: a PERIOD INFILTRATION value linked to a TS
             keeps tracking that TS in a later period whose own PERIOD
             block reappears (for STATUS/CONCENTRATION) without repeating
             INFILTRATION. INFILTRATION is a package-specific (non-AUX)
             field, wired in gwt-uzt.f90/gwe-uze.f90 -- confirms the
             generalized apply_period_settings mechanism generalizes
             without per-package changes. The UZF infiltration flow rate
             isn't a simple prescribed constant (variably saturated flow
             accounting can adjust actual infiltration vs. FINF), so this
             case derives its expected scale factor from period 1's own
             observed value rather than assuming an exact analytic
             constant, then checks periods 2/3 track the TS ratios
             (2x, 3x) relative to that.

Cases:
  uzt_01ts       — CONCENTRATION TS + PACKAGEDATA AUX TS for GWT-UZT
  uzt_01paux     — PERIOD AUXILIARY literal for GWT-UZT
  uze_01ts       — TEMPERATURE TS + PACKAGEDATA AUX TS for GWE-UZE
  uze_01paux     — PERIOD AUXILIARY literal for GWE-UZE
  uztinflcont    — INFILTRATION TS continues across periods, GWT-UZT
  uzeinflcont    — INFILTRATION TS continues across periods, GWE-UZE
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "uzt_01ts",
    "uzt_01paux",
    "uze_01ts",
    "uze_01paux",
    "uztinflcont",
    "uzeinflcont",
]

TS_STEP = 0.05  # TS step-change time within the 0.1-day period
CONC_TEMP_EARLY = 25.0  # UZF cell conc/temp (CONSTANT status) for steps 1-5
CONC_TEMP_LATE = 0.0  # UZF cell conc/temp for steps 6-10
AUX1_EARLY = 99.0
AUX1_LATE = 0.0
PKGDATA_AUX1 = 99.0
PERIOD_AUX1 = 55.0
INFL_EARLY = 25.0  # INFILTRATION conc/temp, steps 1-5
INFL_LATE = 0.0  # steps 6-10
# GWF UZF cell 2 (0-based iuzno=1, mf6 ifno=2) is the only non-CHD cell;
# INFILTRATION obs = finf * area * conc = 1e-5 * 10000 * conc
Q_INFL = 1e-5 * 10000.0

# GWE EST heat-capacity parameters
CPW = 4183.0
RHOW = 999.728
LHV = 2500.0
CPS = 800.0
RHOS = 2650.0

# cross-period persistence: distinct value per period, 3 periods x 10 steps
TS_VAL_P1 = 10.0
TS_VAL_P2 = 20.0
TS_VAL_P3 = 30.0


def _grid_params():
    nlay, nrow, ncol = 1, 1, 3
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = [-10.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    return nlay, nrow, ncol, delr, delc, top, botm, idomain


def _build_gwf(sim, name, aux_name, delr, delc, nlay, nrow, ncol, top, botm, idomain):
    """Build GWF+UZF model; aux_name is 'CONCENTRATION' or 'TEMPERATURE'."""
    gwfname = "gwf_" + name
    nouter, ninner = 600, 100
    hclose, rclose, relax = 1e-6, 0.1, 1.0
    kv = 1e-4

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions="NEWTON",
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
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
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kv, k33=kv)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-6,
        sy=0.35,
        steady_state={0: False},
        transient={0: True},
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), 0.0, 0.0],
            [(0, 0, ncol - 1), 0.0, 0.0],
        ],
        print_input=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    # UZF: one surface cell per column, no vertical connection, no ET
    nuzfcells = ncol
    uzf_pkdat = [
        (iuzno, (0, 0, iuzno), 1, -1, 1.0, kv, 0.2, 0.4, 0.3, 3.5)
        for iuzno in range(nuzfcells)
    ]
    # iuzno, finf, pet, extdp, extwc, ha, hroot, rootact
    # finf=1e-5: cells 0 and ncol-1 are CHD-saturated (rejected regardless);
    # middle cell(s) actually infiltrate, used for the INFILTRATION TS check
    uzf_spd = [
        [iuzno, 1e-5, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0] for iuzno in range(nuzfcells)
    ]
    flopy.mf6.ModflowGwfuzf(
        gwf,
        print_input=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord=f"{gwfname}.uzf.bud",
        nuzfcells=nuzfcells,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        pname="UZF-1",
        filename=f"{gwfname}.uzf",
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
    """Build GWT model without UZT; returns (gwt, gwtname)."""
    gwtname = "gwt_" + name
    gwfname = "gwf_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
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
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.35, filename=f"{gwtname}.sto")
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
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return gwt, gwtname


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
    """Build GWE model without UZE; returns (gwe, gwename)."""
    gwename = "gwe_" + name
    gwfname = "gwf_" + name
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
        porosity=0.35,
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
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )
    return gwe, gwename


def _build_gwf_default(
    sim, name, aux_name, delr, delc, nlay, nrow, ncol, top, botm, idomain
):
    """Build GWF+UZF model with flopy's default IMS settings (no custom
    solver tuning), used by the cross-period persistence cases."""
    gwfname = "gwf_" + name
    kv = 1e-4
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions="NEWTON",
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
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
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kv, k33=kv)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-6,
        sy=0.35,
        steady_state={0: False},
        transient={0: True},
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 0.0, 0.0], [(0, 0, ncol - 1), 0.0, 0.0]],
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    nuzfcells = ncol
    uzf_pkdat = [
        (iuzno, (0, 0, iuzno), 1, -1, 1.0, kv, 0.2, 0.4, 0.3, 3.5)
        for iuzno in range(nuzfcells)
    ]
    uzf_spd = [
        [iuzno, 1e-5, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0] for iuzno in range(nuzfcells)
    ]
    flopy.mf6.ModflowGwfuzf(
        gwf,
        save_flows=True,
        budget_filerecord=f"{gwfname}.uzf.bud",
        nuzfcells=nuzfcells,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        pname="UZF-1",
        auxiliary=[aux_name],
        filename=f"{gwfname}.uzf",
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
    """Build GWT model without UZT, flopy IMS defaults; returns (gwt, gwtname)."""
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
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
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
    """Build GWE model without UZE, flopy IMS defaults; returns (gwe, gwename)."""
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
        porosity=0.35,
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


def _infil_perioddata(depvarkey):
    return {
        0: [
            (0, "STATUS", "ACTIVE"),
            (0, depvarkey, 0.0),
            (0, "INFILTRATION", "infil_ts"),
        ],
        1: [
            (0, "STATUS", "ACTIVE"),
            (0, depvarkey, 0.0),
            # no INFILTRATION setting -- period 0's TS-linked value should
            # persist and keep tracking the TS
        ],
        2: [(0, "STATUS", "ACTIVE"), (0, depvarkey, 0.0)],
    }


def build_models_uztinflcont(name, ws):
    """Case uztinflcont: PERIOD INFILTRATION linked to a TS keeps tracking
    it in a later period whose own PERIOD block doesn't mention it,
    GWT-UZT."""
    nlay, nrow, ncol = 1, 1, 3
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = [-10.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    nuzfcells = ncol

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(0.1, 10, 1.0)] * 3
    )
    _build_gwf_default(
        sim, name, "CONCENTRATION", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwt, gwtname = _build_gwt_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    uztpackagedata = [(iuz, 0.0, 999.0, f"myuzt{iuz + 1}") for iuz in range(nuzfcells)]
    uzt_obs = {(gwtname + ".uzt.obs.csv",): [("uzt1-infl", "infiltration", "myuzt1")]}
    uzt = flopy.mf6.modflow.ModflowGwtuzt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".uzt.bud",
        packagedata=uztpackagedata,
        uztperioddata=_infil_perioddata("CONCENTRATION"),
        observations=uzt_obs,
        pname="UZF-1",
        auxiliary=["aux1"],
    )
    uzt.ts.initialize(
        filename="infil.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (0.1, TS_VAL_P2),
            (0.2, TS_VAL_P3),
            (0.3, TS_VAL_P3),
        ],
        time_series_namerecord="infil_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_uzeinflcont(name, ws):
    """Case uzeinflcont: PERIOD INFILTRATION continues across periods, GWE-UZE."""
    nlay, nrow, ncol = 1, 1, 3
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = [-10.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    nuzfcells = ncol

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(0.1, 10, 1.0)] * 3
    )
    _build_gwf_default(
        sim, name, "TEMPERATURE", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwe, gwename = _build_gwe_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    uzepackagedata = [(iuz, 0.0, 999.0, f"myuze{iuz + 1}") for iuz in range(nuzfcells)]
    uze_obs = {(gwename + ".uze.obs.csv",): [("uze1-infl", "infiltration", "myuze1")]}
    uze = flopy.mf6.ModflowGweuze(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".uze.bud",
        packagedata=uzepackagedata,
        uzeperioddata=_infil_perioddata("TEMPERATURE"),
        observations=uze_obs,
        flow_package_name="UZF-1",
        flow_package_auxiliary_name="TEMPERATURE",
        pname="UZF-1",
        auxiliary=["aux1"],
    )
    uze.ts.initialize(
        filename="infil.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (0.1, TS_VAL_P2),
            (0.2, TS_VAL_P3),
            (0.3, TS_VAL_P3),
        ],
        time_series_namerecord="infil_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_uzt_ts(name, ws):
    """Case uzt_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWT-UZT."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97
    nuzfcells = ncol

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
        sim, name, "CONCENTRATION", delr, delc, nlay, nrow, ncol, top, botm, idomain
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

    uzt_obs = {
        (gwtname + ".uzt.obs.csv",): [
            ("uzt2-conc", "CONCENTRATION", 2),
            ("uzt2-infl", "INFILTRATION", 2),
        ],
    }
    uzt_obs["digits"] = 7
    uzt_obs["print_input"] = True
    uzt_obs["filename"] = gwtname + ".uzt.obs"

    # PACKAGEDATA: (ifno, strt, aux1_ts, aux2, boundname)
    uztpackagedata = [
        (iuz, 0.0, "aux1_ts", 999.0, f"myuzt{iuz + 1}") for iuz in range(nuzfcells)
    ]

    # STATUS=CONSTANT; CONCENTRATION and INFILTRATION (uzt_ad_ts's own
    # re-sync) link to TS. Feature 2 (0-based ifno=1) is the only cell not
    # CHD-adjacent, so it's the one with a measurable INFILTRATION flow.
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
        uztperioddata=[
            (1, "STATUS", "CONSTANT"),
            (1, "CONCENTRATION", "conc_ts"),
            (1, "INFILTRATION", "infl_ts"),
        ],
        observations=uzt_obs,
        pname="UZF-1",
        auxiliary=["aux1", "aux2"],
    )

    uzt.ts.initialize(
        filename=f"{gwtname}.uzt.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, INFL_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, INFL_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, INFL_LATE),
        ],
        time_series_namerecord=[("conc_ts", "aux1_ts", "infl_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_uzt_peraux(name, ws):
    """Case uzt_01paux: PERIOD AUXILIARY literal (T3) for GWT-UZT."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97
    nuzfcells = ncol

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
        sim, name, "CONCENTRATION", delr, delc, nlay, nrow, ncol, top, botm, idomain
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

    uztpackagedata = [
        (iuz, 0.0, PKGDATA_AUX1, 999.0, f"myuzt{iuz + 1}") for iuz in range(nuzfcells)
    ]

    flopy.mf6.modflow.ModflowGwtuzt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".uzt.bin",
        budget_filerecord=gwtname + ".uzt.bud",
        packagedata=uztpackagedata,
        uztperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "CONCENTRATION", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="UZF-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


def build_models_uze_ts(name, ws):
    """Case uze_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWE-UZE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97
    nuzfcells = ncol

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
        sim, name, "TEMPERATURE", delr, delc, nlay, nrow, ncol, top, botm, idomain
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

    uze_obs = {
        (gwename + ".uze.obs.csv",): [
            ("uze2-temp", "TEMPERATURE", 2),
            ("uze2-infl", "INFILTRATION", 2),
        ],
    }
    uze_obs["digits"] = 7
    uze_obs["print_input"] = True
    uze_obs["filename"] = gwename + ".uze.obs"

    uzepackagedata = [
        (iuz, 0.0, "aux1_ts", 999.0, f"myuze{iuz + 1}") for iuz in range(nuzfcells)
    ]

    # feature 2 (0-based ifno=1) is the only cell not CHD-adjacent
    uze = flopy.mf6.modflow.ModflowGweuze(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".uze.bin",
        budget_filerecord=gwename + ".uze.bud",
        packagedata=uzepackagedata,
        uzeperioddata=[
            (1, "STATUS", "CONSTANT"),
            (1, "TEMPERATURE", "temp_ts"),
            (1, "INFILTRATION", "infl_ts"),
        ],
        observations=uze_obs,
        pname="UZF-1",
        auxiliary=["aux1", "aux2"],
    )

    uze.ts.initialize(
        filename=f"{gwename}.uze.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, INFL_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, INFL_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, INFL_LATE),
        ],
        time_series_namerecord=[("temp_ts", "aux1_ts", "infl_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_uze_peraux(name, ws):
    """Case uze_01paux: PERIOD AUXILIARY literal (T3) for GWE-UZE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    nper = 1
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97
    nuzfcells = ncol

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
        sim, name, "TEMPERATURE", delr, delc, nlay, nrow, ncol, top, botm, idomain
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

    uzepackagedata = [
        (iuz, 0.0, PKGDATA_AUX1, 999.0, f"myuze{iuz + 1}") for iuz in range(nuzfcells)
    ]

    flopy.mf6.modflow.ModflowGweuze(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".uze.bin",
        budget_filerecord=gwename + ".uze.bud",
        packagedata=uzepackagedata,
        uzeperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "TEMPERATURE", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="UZF-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


_BUILD_FNS = {
    "uzt_01ts": build_models_uzt_ts,
    "uzt_01paux": build_models_uzt_peraux,
    "uze_01ts": build_models_uze_ts,
    "uze_01paux": build_models_uze_peraux,
    "uztinflcont": build_models_uztinflcont,
    "uzeinflcont": build_models_uzeinflcont,
}


def build_models(idx, test):
    name = cases[idx]
    return _BUILD_FNS[name](name, test.workspace)


def _check_ts(bud_fname, obs_fname, pkg_label, conc_col, infl_col, thermal=False):
    """Shared T1+T2 check for UZT and UZE. infl_col checks uzt_ad_ts/
    uze_ad_ts's own INFILTRATION TS re-sync, not just the base-class one."""
    # ------------------------------------------------------------------ T1 --
    assert os.path.isfile(obs_fname), f"{pkg_label} obs file not found: {obs_fname}"
    tc = np.genfromtxt(obs_fname, names=True, delimiter=",")

    conc_obs = tc[conc_col]
    assert len(conc_obs) == 10, f"Expected 10 obs, got {len(conc_obs)}"
    assert np.allclose(conc_obs[:5], CONC_TEMP_EARLY, rtol=1e-5), (
        f"{pkg_label} steps 1-5 conc/temp expected {CONC_TEMP_EARLY}, "
        f"got {conc_obs[:5]}"
    )
    assert np.allclose(conc_obs[5:], CONC_TEMP_LATE, atol=1e-10), (
        f"{pkg_label} steps 6-10 conc/temp expected {CONC_TEMP_LATE}, "
        f"got {conc_obs[5:]}"
    )

    # INFILTRATION obs is a mass/energy flow rate (Q_INFL * conc, scaled by
    # RHOW*CPW for the thermal/GWE case), not the conc/temp itself
    scale = RHOW * CPW if thermal else 1.0
    infl_early = Q_INFL * INFL_EARLY * scale
    infl_obs = tc[infl_col]
    assert np.allclose(infl_obs[:5], infl_early, rtol=1e-5), (
        f"{pkg_label} steps 1-5 infiltration expected {infl_early}, got {infl_obs[:5]}"
    )
    assert np.allclose(infl_obs[5:], 0.0, atol=1e-10 * max(scale, 1.0)), (
        f"{pkg_label} steps 6-10 infiltration expected 0.0, got {infl_obs[5:]}"
    )

    # ------------------------------------------------------------------ T2 --
    assert os.path.isfile(bud_fname), f"{pkg_label} budget not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY records, got {len(aux_records)}"
    )
    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, f"AUX1 not in AUXILIARY record fields: {fnames}"

    aux1_early = np.array([rec["AUX1"][0] for rec in aux_records[:5]])
    aux1_late = np.array([rec["AUX1"][0] for rec in aux_records[5:]])
    assert np.allclose(aux1_early, AUX1_EARLY, rtol=1e-5), (
        f"{pkg_label} steps 1-5 AUX1 expected {AUX1_EARLY}, got {aux1_early}"
    )
    assert np.allclose(aux1_late, AUX1_LATE, atol=1e-10), (
        f"{pkg_label} steps 6-10 AUX1 expected {AUX1_LATE}, got {aux1_late}"
    )


def _check_peraux(bud_fname, pkg_label):
    """Shared T3 check for UZT and UZE."""
    assert os.path.isfile(bud_fname), f"{pkg_label} budget not found: {bud_fname}"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY records, got {len(aux_records)}"
    )
    fnames = aux_records[0].dtype.names
    assert "AUX1" in fnames, f"AUX1 not in AUXILIARY record fields: {fnames}"

    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    assert np.allclose(aux1_vals, PERIOD_AUX1, rtol=1e-5), (
        f"{pkg_label} AUX1 expected {PERIOD_AUX1}, got {aux1_vals} "
        f"— PACKAGEDATA baseline {PKGDATA_AUX1} leaked through"
    )


def check_output_uzt_ts(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_ts(
        os.path.join(ws, f"{gwtname}.uzt.bud"),
        os.path.join(ws, f"{gwtname}.uzt.obs.csv"),
        "UZT",
        "UZT2CONC",
        "UZT2INFL",
    )


def check_output_uzt_peraux(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_peraux(os.path.join(ws, f"{gwtname}.uzt.bud"), "UZT")


def check_output_uze_ts(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_ts(
        os.path.join(ws, f"{gwename}.uze.bud"),
        os.path.join(ws, f"{gwename}.uze.obs.csv"),
        "UZE",
        "UZE2TEMP",
        "UZE2INFL",
        thermal=True,
    )


def check_output_uze_peraux(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_peraux(os.path.join(ws, f"{gwename}.uze.bud"), "UZE")


def _check_inflcont(test, gwname, ext, budname):
    fname = test.workspace / f"{gwname}.{ext}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")
    infl_obs = tc["UZT1INFL"] if "UZT1INFL" in tc.dtype.names else tc["UZE1INFL"]
    assert len(infl_obs) == 30, (
        f"Expected 30 obs records (10/period, 3 periods), got {len(infl_obs)}"
    )
    period1, period2, period3 = infl_obs[:10], infl_obs[10:20], infl_obs[20:]

    assert np.allclose(period1, period1[0]), (
        f"{budname} period 1 INFILTRATION obs not constant: {period1}"
    )
    scale = period1[0] / TS_VAL_P1

    exp2 = scale * TS_VAL_P2
    exp3 = scale * TS_VAL_P3
    assert np.allclose(period2, exp2), (
        f"{budname} period 2 INFILTRATION obs expected {exp2} (TS should "
        f"still be tracked even though period 1's PERIOD block doesn't "
        f"mention INFILTRATION), got {period2}"
    )
    assert np.allclose(period3, exp3), (
        f"{budname} period 3 INFILTRATION obs expected {exp3}, got {period3}"
    )


def check_output_uztinflcont(test):
    _check_inflcont(test, "gwt_" + test.name, "uzt", "UZT")


def check_output_uzeinflcont(test):
    _check_inflcont(test, "gwe_" + test.name, "uze", "UZE")


_CHECK_FNS = {
    "uzt_01ts": check_output_uzt_ts,
    "uzt_01paux": check_output_uzt_peraux,
    "uze_01ts": check_output_uze_ts,
    "uze_01paux": check_output_uze_peraux,
    "uztinflcont": check_output_uztinflcont,
    "uzeinflcont": check_output_uzeinflcont,
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
