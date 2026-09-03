"""
MWT/MWE time-series, PERIOD AUXILIARY, and cross-period persistence
integration tests.

Covers GWT-MWT and GWE-MWE in paired cases that share the same GWF+MAW model.

  PERIOD TS: CONCENTRATION/TEMPERATURE and RATE (mwt_ad_ts/mwe_ad_ts's own
             re-sync) linked to a STEPWISE TS that changes mid-period. The
             well's flow is injection (RATE > 0) so the RATE obs reflects it.

  PACKAGEDATA AUX TS: AUX1 carries TS name "aux1_ts"; the binary budget
             AUXILIARY term verifies per-timestep re-sync.

  PERIOD AUXILIARY literal: PERIOD AUXILIARY overrides PACKAGEDATA baseline
             (99.0 -> 55.0). Binary budget AUXILIARY must show 55.0.

  Cross-period persistence: a PERIOD setting issued once must persist into
             a later period whose own PERIOD block doesn't repeat it -- for
             AUXILIARY (generic, apply_period_auxiliary), RATE (package-
             specific mem_setptr wiring in gwt-mwt.f90/gwe-mwe.f90), and
             CONCENTRATION/TEMPERATURE (shared tsp-apt.f90 base class,
             gated by STATUS).

  AUX/PACKAGEDATA-TS collision: a literal PERIOD AUXILIARY override must
             not be silently overwritten by a still-active PACKAGEDATA AUX
             time series on the next timestep. Depends on
             StructArray.f90's ts_update AUX branch addressing the TS link
             by its position in the AUX array, not the row-schema column,
             so a superseding PERIOD override can find and remove it.

Cases:
  mwt_01ts     — RATE TS + PACKAGEDATA AUX TS, GWT-MWT
  mwt_01paux   — PERIOD AUXILIARY literal, GWT-MWT
  mwe_01ts     — RATE TS + PACKAGEDATA AUX TS, GWE-MWE
  mwe_01paux   — PERIOD AUXILIARY literal, GWE-MWE
  mwtauxper    — AUXILIARY literal persists across periods, GWT-MWT
  mweauxper    — AUXILIARY literal persists across periods, GWE-MWE
  mwtauxovr    — literal AUXILIARY override outlives an active PACKAGEDATA
                 AUX TS, GWT-MWT
  mwtconccont  — CONCENTRATION TS continues across periods, GWT-MWT
  mweconccont  — TEMPERATURE TS continues across periods, GWE-MWE
  mwtratecont  — RATE TS continues across periods, GWT-MWT
  mweratecont  — RATE TS continues across periods, GWE-MWE
  mwtcontinue  — AUXILIARY TS continues across periods, GWT-MWT
  mwecontinue  — AUXILIARY TS continues across periods, GWE-MWE
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "mwt_01ts",
    "mwt_01paux",
    "mwe_01ts",
    "mwe_01paux",
    "mwtauxper",
    "mweauxper",
    "mwtauxovr",
    "mwtconccont",
    "mweconccont",
    "mwtratecont",
    "mweratecont",
    "mwtcontinue",
    "mwecontinue",
]

TS_STEP = 0.05  # TS step-change time within the 0.1-day period
CONC_TEMP_EARLY = 25.0  # well conc/temp (CONSTANT status) for steps 1-5
CONC_TEMP_LATE = 0.0  # well conc/temp for steps 6-10
AUX1_EARLY = 99.0
AUX1_LATE = 0.0
PKGDATA_AUX1 = 99.0
PERIOD_AUX1 = 55.0
RATE_EARLY = 25.0  # RATE conc/temp, steps 1-5
RATE_LATE = 0.0  # steps 6-10
Q_INJ = 0.5  # GWF MAW well 0 injection rate; RATE obs = Q_INJ * conc

# cross-period persistence: distinct value per period, 3 periods x 5 steps
TS_VAL_P1 = 10.0
TS_VAL_P2 = 20.0
TS_VAL_P3 = 30.0

# aux_pkgdata_override: PACKAGEDATA AUX TS values, distinct from the
# literal PERIOD override so a leak back to the TS is unambiguous
PKGD_AUX_VALS = [99.0, 150.0, 200.0]

# GWE EST heat-capacity parameters
CPW = 4183.0
RHOW = 999.728
LHV = 2500.0
CPS = 800.0
RHOS = 2650.0
KTF = 1.0
FTHK = 0.1


def _grid_params():
    nlay, nrow, ncol = 3, 1, 5
    lx, lz = 5.0, 3.0
    delr = lx / ncol
    delc = 1.0
    delz = lz / nlay
    top = [0.0] * ncol
    botm = list(0 - np.arange(delz, nlay * delz + delz, delz))
    idomain = np.full((nlay, nrow, ncol), 1)
    return nlay, nrow, ncol, delr, delc, top, botm, idomain


def _build_gwf(
    sim,
    name,
    aux_name,
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
    """Build GWF+MAW model; aux_name is 'CONCENTRATION' or 'TEMPERATURE'."""
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
            [(0, 0, 0), 0.0, 100.0],
            [(0, 0, ncol - 1), 0.0, 0.0],
        ],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    wellbottom = -3.0
    flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{gwfname}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        packagedata=[[0, 0.1, wellbottom, 0.0, "THIEM", 3]],
        connectiondata=[
            [0, 0, (0, 0, 2), 0.0, -1, 1.0, 0.1],
            [0, 1, (1, 0, 2), -1.0, -2, 1.0, 0.1],
            [0, 2, (2, 0, 2), -2.0, -3, 1.0, 0.1],
        ],
        perioddata=[[0, "rate", 0.5]],  # injection: RATE conc/temp only
        # applies (via concrate) when qbnd > 0; see mwt_rate_term
        pname="MAW-1",
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
    """Build GWT model without MWT; returns (gwt, gwtname)."""
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
    """Build GWE model without MWE; returns (gwe, gwename)."""
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
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )
    return gwe, gwename


def _build_gwf_default(
    sim, name, aux_name, delr, delc, nlay, nrow, ncol, top, botm, idomain, well_rate=0.5
):
    """Build GWF+MAW model with flopy's default IMS settings (no custom
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
        stress_period_data=[[(0, 0, 0), 0.0, 100.0], [(0, 0, ncol - 1), 0.0, 0.0]],
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    wellbottom = -3.0
    flopy.mf6.ModflowGwfmaw(
        gwf,
        packagedata=[[0, 0.1, wellbottom, 0.0, "THIEM", 3]],
        connectiondata=[
            [0, 0, (0, 0, 2), 0.0, -1, 1.0, 0.1],
            [0, 1, (1, 0, 2), -1.0, -2, 1.0, 0.1],
            [0, 2, (2, 0, 2), -2.0, -3, 1.0, 0.1],
        ],
        perioddata=[[0, "rate", well_rate]],
        pname="MAW-1",
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
    """Build GWT model without MWT, flopy IMS defaults; returns (gwt, gwtname)."""
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
    """Build GWE model without MWE, flopy IMS defaults; returns (gwe, gwename)."""
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


def build_models_mwt_ts(name, ws):
    """Case mwt_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWT-MWT."""
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
        "CONCENTRATION",
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

    mwt_obs = {
        (gwtname + ".mwt.obs.csv",): [
            ("mwt1-conc", "CONCENTRATION", 1),
            ("mwt1-rate", "RATE", 1),
        ],
    }
    mwt_obs["digits"] = 7
    mwt_obs["print_input"] = True
    mwt_obs["filename"] = gwtname + ".mwt.obs"

    # STATUS=CONSTANT; CONCENTRATION and RATE (mwt_ad_ts's own re-sync) link
    # to TS; PACKAGEDATA aux1 links to TS
    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".mwt.bin",
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, "aux1_ts", 999.0, "mymwt1")],
        mwtperioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", "conc_ts"),
            (0, "RATE", "rate_ts"),
        ],
        observations=mwt_obs,
        pname="MAW-1",
        auxiliary=["aux1", "aux2"],
    )

    mwt.ts.initialize(
        filename=f"{gwtname}.mwt.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, RATE_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, RATE_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, RATE_LATE),
        ],
        time_series_namerecord=[("conc_ts", "aux1_ts", "rate_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_mwt_peraux(name, ws):
    """Case mwt_01paux: PERIOD AUXILIARY literal (T3) for GWT-MWT."""
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
        "CONCENTRATION",
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

    flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".mwt.bin",
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, PKGDATA_AUX1, 999.0, "mymwt1")],
        mwtperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "CONCENTRATION", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="MAW-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


def build_models_mwe_ts(name, ws):
    """Case mwe_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWE-MWE."""
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
        "TEMPERATURE",
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

    mwe_obs = {
        (gwename + ".mwe.obs.csv",): [
            ("mwe1-temp", "TEMPERATURE", 1),
            ("mwe1-rate", "RATE", 1),
        ],
    }
    mwe_obs["digits"] = 7
    mwe_obs["print_input"] = True
    mwe_obs["filename"] = gwename + ".mwe.obs"

    # STATUS=CONSTANT; TEMPERATURE and RATE (mwe_ad_ts's own re-sync) link
    # to TS; PACKAGEDATA aux1 links to TS
    mwe = flopy.mf6.modflow.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".mwe.bin",
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, "aux1_ts", 999.0, "mymwe1")],
        mweperioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "TEMPERATURE", "temp_ts"),
            (0, "RATE", "rate_ts"),
        ],
        observations=mwe_obs,
        pname="MAW-1",
        auxiliary=["aux1", "aux2"],
    )

    mwe.ts.initialize(
        filename=f"{gwename}.mwe.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, RATE_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, RATE_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, RATE_LATE),
        ],
        time_series_namerecord=[("temp_ts", "aux1_ts", "rate_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_mwe_peraux(name, ws):
    """Case mwe_01paux: PERIOD AUXILIARY literal (T3) for GWE-MWE."""
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
        "TEMPERATURE",
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

    flopy.mf6.modflow.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".mwe.bin",
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, PKGDATA_AUX1, 999.0, "mymwe1")],
        mweperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "TEMPERATURE", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="MAW-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


def build_models_mwt_auxper(name, ws):
    """Case mwtauxper: PERIOD AUXILIARY literal persists into a later period
    whose own PERIOD block doesn't repeat it, GWT-MWT."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=2, perioddata=[(1.0, 5, 1.0)] * 2
    )
    _build_gwf_default(
        sim, name, "CONCENTRATION", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwt, gwtname = _build_gwt_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mwtperioddata = {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", 100.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        1: [
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", 100.0),
        ],
    }
    flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, PKGDATA_AUX1, "mymwt1")],
        mwtperioddata=mwtperioddata,
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    return sim, None


def build_models_mwe_auxper(name, ws):
    """Case mweauxper: PERIOD AUXILIARY literal persists into a later period
    whose own PERIOD block doesn't repeat it, GWE-MWE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=2, perioddata=[(1.0, 5, 1.0)] * 2
    )
    _build_gwf_default(
        sim, name, "TEMPERATURE", delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    gwe, gwename = _build_gwe_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mweperioddata = {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, "TEMPERATURE", 100.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        1: [
            (0, "STATUS", "CONSTANT"),
            (0, "TEMPERATURE", 100.0),
        ],
    }
    flopy.mf6.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, PKGDATA_AUX1, "mymwe1")],
        mweperioddata=mweperioddata,
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    return sim, None


def build_models_mwt_auxovr(name, ws):
    """Case mwtauxovr: a literal PERIOD AUXILIARY override, issued once,
    must persist across later periods that don't reissue it -- not
    silently revert to an active PACKAGEDATA AUX time series, GWT-MWT.

    Exercises StructArray.f90's ts_update: apply_period_auxiliary() finds
    the well's existing PACKAGEDATA-level AUX TS link by its position in
    the AUX array (1..naux), matching how the link itself is registered,
    so it can remove the stale link once the PERIOD override supersedes
    it -- otherwise TimeSeriesManager's per-timestep resync would silently
    overwrite the literal override back to the PACKAGEDATA TS value on
    the very next timestep.

    mwtauxper (above) doesn't exercise this path: its PACKAGEDATA AUX
    value is a plain literal, so no TS file is attached (ts_active is
    false) and the resync loop never runs.
    """
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    ts_times = [0.0, 1.0, 2.0, 3.0]
    ts_pkgd_aux = [PKGD_AUX_VALS[0]] + PKGD_AUX_VALS

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

    mwtperioddata = {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", 100.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        # periods 1, 2 reissue STATUS/CONCENTRATION only -- AUXILIARY is
        # never mentioned again and must keep the period-0 override, not
        # revert to the still-active PACKAGEDATA AUX time series
        1: [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 100.0)],
        2: [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 100.0)],
    }
    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, "pkgd_aux1", "mymwt1")],
        mwtperioddata=mwtperioddata,
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    mwt.ts.initialize(
        filename=f"{gwtname}.mwt.ts",
        timeseries=list(zip(ts_times, ts_pkgd_aux)),
        time_series_namerecord=["pkgd_aux1"],
        interpolation_methodrecord=["linearend"],
    )
    return sim, None


def build_models_mwt_conccont(name, ws):
    """Case mwtconccont: PERIOD CONCENTRATION linked to a TS keeps tracking
    it in a later period whose own PERIOD block reappears (for RATE)
    without repeating CONCENTRATION, GWT-MWT.

    CONCENTRATION/TEMPERATURE is handled directly in the shared
    tsp-apt.f90 base class, not per-package code, and is gated by STATUS.
    """
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim,
        name,
        "CONCENTRATION",
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        well_rate=-1.0,
    )
    gwt, gwtname = _build_gwt_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mwt_obs = {(gwtname + ".mwt.obs.csv",): [("mwt1-conc", "CONCENTRATION", 1)]}
    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, "mymwt1")],
        mwtperioddata=_conc_perioddata("CONCENTRATION"),
        observations=mwt_obs,
        pname="MAW-1",
    )
    mwt.ts.initialize(
        filename="conc.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="conc_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_mwe_conccont(name, ws):
    """Case mweconccont: PERIOD TEMPERATURE continues across periods, GWE-MWE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim,
        name,
        "TEMPERATURE",
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        well_rate=-1.0,
    )
    gwe, gwename = _build_gwe_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mwe_obs = {(gwename + ".mwe.obs.csv",): [("mwe1-conc", "TEMPERATURE", 1)]}
    mwe = flopy.mf6.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, "mymwe1")],
        mweperioddata=_conc_perioddata("TEMPERATURE"),
        observations=mwe_obs,
        pname="MAW-1",
    )
    mwe.ts.initialize(
        filename="conc.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="conc_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def _conc_perioddata(depvarkey):
    return {
        0: [(0, "STATUS", "CONSTANT"), (0, depvarkey, "conc_ts")],
        1: [
            (0, "STATUS", "CONSTANT"),
            (0, "RATE", 0.0),
            # no CONCENTRATION/TEMPERATURE setting -- period 0's TS-linked
            # value should persist and keep tracking the TS
        ],
        2: [(0, "STATUS", "CONSTANT"), (0, "RATE", 0.0)],
    }


def build_models_mwt_ratecont(name, ws):
    """Case mwtratecont: PERIOD RATE linked to a TS keeps tracking it in a
    later period whose own PERIOD block doesn't mention it, GWT-MWT. RATE
    is a package-specific (non-AUX) field wired in gwt-mwt.f90."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim,
        name,
        "CONCENTRATION",
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        well_rate=Q_INJ,
    )
    gwt, gwtname = _build_gwt_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mwt_obs = {(gwtname + ".mwt.obs.csv",): [("mwt1-rate", "RATE", 1)]}
    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, 999.0, "mymwt1")],
        mwtperioddata=_rate_perioddata("CONCENTRATION"),
        observations=mwt_obs,
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    mwt.ts.initialize(
        filename="rate.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="rate_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_mwe_ratecont(name, ws):
    """Case mweratecont: PERIOD RATE continues across periods, GWE-MWE."""
    nlay, nrow, ncol, delr, delc, top, botm, idomain = _grid_params()
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    _build_gwf_default(
        sim,
        name,
        "TEMPERATURE",
        delr,
        delc,
        nlay,
        nrow,
        ncol,
        top,
        botm,
        idomain,
        well_rate=Q_INJ,
    )
    gwe, gwename = _build_gwe_base_default(
        sim, name, delr, delc, nlay, nrow, ncol, top, botm, idomain
    )
    mwe_obs = {(gwename + ".mwe.obs.csv",): [("mwe1-rate", "RATE", 1)]}
    mwe = flopy.mf6.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, 999.0, "mymwe1")],
        mweperioddata=_rate_perioddata("TEMPERATURE"),
        observations=mwe_obs,
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    mwe.ts.initialize(
        filename="rate.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="rate_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def _rate_perioddata(depvarkey):
    return {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, depvarkey, 0.0),
            (0, "RATE", "rate_ts"),
        ],
        1: [
            (0, "STATUS", "CONSTANT"),
            (0, depvarkey, 0.0),
            # no RATE setting -- period 0's TS-linked value should persist
            # and keep tracking the TS
        ],
        2: [(0, "STATUS", "CONSTANT"), (0, depvarkey, 0.0)],
    }


def build_models_mwt_continue(name, ws):
    """Case mwtcontinue: PERIOD AUXILIARY linked to a TS keeps tracking it
    in a later period whose own PERIOD block doesn't mention it, GWT-MWT."""
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
    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=[(0, 0.0, 999.0, "mymwt1")],
        mwtperioddata=_auxts_perioddata("CONCENTRATION"),
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    mwt.ts.initialize(
        filename="aux1.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="aux1_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_mwe_continue(name, ws):
    """Case mwecontinue: PERIOD AUXILIARY continues across periods, GWE-MWE."""
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
    mwe = flopy.mf6.ModflowGwemwe(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".mwe.bud",
        packagedata=[(0, 0.0, KTF, FTHK, 999.0, "mymwe1")],
        mweperioddata=_auxts_perioddata("TEMPERATURE"),
        pname="MAW-1",
        auxiliary=["aux1"],
    )
    mwe.ts.initialize(
        filename="aux1.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="aux1_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def _auxts_perioddata(depvarkey):
    return {
        0: [
            (0, "STATUS", "CONSTANT"),
            (0, depvarkey, 100.0),
            (0, "AUXILIARY", "aux1", "aux1_ts"),
        ],
        1: [(0, "STATUS", "CONSTANT"), (0, depvarkey, 100.0)],
        2: [(0, "STATUS", "CONSTANT"), (0, depvarkey, 100.0)],
    }


_BUILD_FNS = {
    "mwt_01ts": build_models_mwt_ts,
    "mwt_01paux": build_models_mwt_peraux,
    "mwe_01ts": build_models_mwe_ts,
    "mwe_01paux": build_models_mwe_peraux,
    "mwtauxper": build_models_mwt_auxper,
    "mweauxper": build_models_mwe_auxper,
    "mwtauxovr": build_models_mwt_auxovr,
    "mwtconccont": build_models_mwt_conccont,
    "mweconccont": build_models_mwe_conccont,
    "mwtratecont": build_models_mwt_ratecont,
    "mweratecont": build_models_mwe_ratecont,
    "mwtcontinue": build_models_mwt_continue,
    "mwecontinue": build_models_mwe_continue,
}


def build_models(idx, test):
    name = cases[idx]
    return _BUILD_FNS[name](name, test.workspace)


def _check_ts(bud_fname, obs_fname, pkg_label, conc_col, rate_col, thermal=False):
    """Shared T1+T2 check for MWT and MWE. rate_col checks mwt_ad_ts/
    mwe_ad_ts's own RATE TS re-sync, not just the base-class one."""
    # ------------------------------------------------------------------ T1 --
    assert os.path.isfile(obs_fname), f"{pkg_label} obs file not found: {obs_fname}"
    tc = np.genfromtxt(obs_fname, names=True, delimiter=",")

    conc_obs = tc[conc_col]
    assert len(conc_obs) == 10, f"Expected 10 obs, got {len(conc_obs)}"
    # Steps 1-5: CONC_TEMP_EARLY (25.0); steps 6-10: CONC_TEMP_LATE (0.0)
    assert np.allclose(conc_obs[:5], CONC_TEMP_EARLY, rtol=1e-5), (
        f"{pkg_label} steps 1-5 conc/temp expected {CONC_TEMP_EARLY}, "
        f"got {conc_obs[:5]}"
    )
    assert np.allclose(conc_obs[5:], CONC_TEMP_LATE, atol=1e-10), (
        f"{pkg_label} steps 6-10 conc/temp expected {CONC_TEMP_LATE}, "
        f"got {conc_obs[5:]}"
    )

    # RATE obs is a mass/energy flow rate (Q_INJ * conc, scaled by RHOW*CPW
    # for the thermal/GWE case), not the conc/temp itself
    scale = RHOW * CPW if thermal else 1.0
    rate_early = Q_INJ * RATE_EARLY * scale
    rate_obs = tc[rate_col]
    assert np.allclose(rate_obs[:5], rate_early, rtol=1e-5), (
        f"{pkg_label} steps 1-5 rate expected {rate_early}, got {rate_obs[:5]}"
    )
    assert np.allclose(rate_obs[5:], 0.0, atol=1e-10 * max(scale, 1.0)), (
        f"{pkg_label} steps 6-10 rate expected 0.0, got {rate_obs[5:]}"
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
    """Shared T3 check for MWT and MWE."""
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


def check_output_mwt_ts(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_ts(
        os.path.join(ws, f"{gwtname}.mwt.bud"),
        os.path.join(ws, f"{gwtname}.mwt.obs.csv"),
        "MWT",
        "MWT1CONC",
        "MWT1RATE",
    )


def check_output_mwt_peraux(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_peraux(os.path.join(ws, f"{gwtname}.mwt.bud"), "MWT")


def check_output_mwe_ts(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_ts(
        os.path.join(ws, f"{gwename}.mwe.bud"),
        os.path.join(ws, f"{gwename}.mwe.obs.csv"),
        "MWE",
        "MWE1TEMP",
        "MWE1RATE",
        thermal=True,
    )


def check_output_mwe_peraux(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_peraux(os.path.join(ws, f"{gwename}.mwe.bud"), "MWE")


def _check_auxper(test, gwname, ext, budname):
    fname = test.workspace / f"{gwname}.{ext}.bud"
    bobj = flopy.utils.CellBudgetFile(fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records (5/period), got {len(aux_records)}"
    )
    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    period1 = aux1_vals[:5]
    period2 = aux1_vals[5:]
    assert np.allclose(period1, PERIOD_AUX1), (
        f"{budname} period 1 AUX1 expected {PERIOD_AUX1}, got {period1}"
    )
    assert np.allclose(period2, PERIOD_AUX1), (
        f"{budname} period 2 AUX1 expected {PERIOD_AUX1} (carried over from "
        f"period 1; period 2's PERIOD block doesn't mention AUXILIARY), "
        f"got {period2}"
    )


def check_output_mwt_auxper(test):
    _check_auxper(test, "gwt_" + test.name, "mwt", "MWT")


def check_output_mwe_auxper(test):
    _check_auxper(test, "gwe_" + test.name, "mwe", "MWE")


def check_output_mwt_auxovr(test):
    fname = test.workspace / f"gwt_{test.name}.mwt.bud"
    bobj = flopy.utils.CellBudgetFile(fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 15, (
        f"Expected 15 AUXILIARY budget records (5/period x 3 periods), "
        f"got {len(aux_records)}"
    )
    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    periods = [aux1_vals[0:5], aux1_vals[5:10], aux1_vals[10:15]]
    expected = [PERIOD_AUX1] * 3
    observed = [p[0] for p in periods]
    assert np.allclose(observed, expected), (
        f"AUX1 expected to keep the period-0 literal override {expected} "
        f"through periods 1/2 (never reissued), got {observed} -- a value "
        f"matching {PKGD_AUX_VALS} would mean the still-active PACKAGEDATA "
        f"AUX time series silently clobbered the override."
    )
    for p in periods:
        assert np.allclose(p, p[0]), "AUX1 changed within a single period"


def _check_conccont(test, gwname, ext, obscol, budname):
    fname = test.workspace / f"{gwname}.{ext}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")
    conc_obs = tc[obscol]
    assert len(conc_obs) == 15, (
        f"Expected 15 obs records (5/period, 3 periods), got {len(conc_obs)}"
    )
    period1, period2, period3 = conc_obs[:5], conc_obs[5:10], conc_obs[10:]
    assert np.allclose(period1, TS_VAL_P1), (
        f"{budname} period 1 CONCENTRATION/TEMPERATURE obs expected "
        f"{TS_VAL_P1}, got {period1}"
    )
    assert np.allclose(period2, TS_VAL_P2), (
        f"{budname} period 2 CONCENTRATION/TEMPERATURE obs expected "
        f"{TS_VAL_P2} (TS should still be tracked even though period 1's "
        f"PERIOD block reappears for RATE without repeating "
        f"CONCENTRATION/TEMPERATURE), got {period2}"
    )
    assert np.allclose(period3, TS_VAL_P3), (
        f"{budname} period 3 CONCENTRATION/TEMPERATURE obs expected "
        f"{TS_VAL_P3}, got {period3}"
    )


def check_output_mwt_conccont(test):
    _check_conccont(test, "gwt_" + test.name, "mwt", "MWT1CONC", "MWT")


def check_output_mwe_conccont(test):
    _check_conccont(test, "gwe_" + test.name, "mwe", "MWE1CONC", "MWE")


def _check_ratecont(test, gwname, ext, obscol, budname, thermal):
    fname = test.workspace / f"{gwname}.{ext}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")
    rate_obs = tc[obscol]
    assert len(rate_obs) == 15, (
        f"Expected 15 obs records (5/period, 3 periods), got {len(rate_obs)}"
    )
    period1, period2, period3 = rate_obs[:5], rate_obs[5:10], rate_obs[10:]
    scale = CPW * RHOW if thermal else 1.0
    exp1 = Q_INJ * TS_VAL_P1 * scale
    exp2 = Q_INJ * TS_VAL_P2 * scale
    exp3 = Q_INJ * TS_VAL_P3 * scale
    assert np.allclose(period1, exp1), (
        f"{budname} period 1 RATE obs expected {exp1}, got {period1}"
    )
    assert np.allclose(period2, exp2), (
        f"{budname} period 2 RATE obs expected {exp2} (TS should still be "
        f"tracked even though period 1's PERIOD block doesn't mention "
        f"RATE), got {period2}"
    )
    assert np.allclose(period3, exp3), (
        f"{budname} period 3 RATE obs expected {exp3}, got {period3}"
    )


def check_output_mwt_ratecont(test):
    _check_ratecont(test, "gwt_" + test.name, "mwt", "MWT1RATE", "MWT", thermal=False)


def check_output_mwe_ratecont(test):
    _check_ratecont(test, "gwe_" + test.name, "mwe", "MWE1RATE", "MWE", thermal=True)


def _check_continue(test, gwname, ext, budname):
    fname = test.workspace / f"{gwname}.{ext}.bud"
    bobj = flopy.utils.CellBudgetFile(fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 15, (
        f"Expected 15 AUXILIARY budget records (5/period, 3 periods), "
        f"got {len(aux_records)}"
    )
    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    period1, period2, period3 = aux1_vals[:5], aux1_vals[5:10], aux1_vals[10:]
    assert np.allclose(period1, TS_VAL_P1), (
        f"{budname} period 1 AUX1 expected {TS_VAL_P1}, got {period1}"
    )
    assert np.allclose(period2, TS_VAL_P2), (
        f"{budname} period 2 AUX1 expected {TS_VAL_P2} (TS should still be "
        f"tracked even though period 1's PERIOD block doesn't mention "
        f"AUXILIARY), got {period2}"
    )
    assert np.allclose(period3, TS_VAL_P3), (
        f"{budname} period 3 AUX1 expected {TS_VAL_P3}, got {period3}"
    )


def check_output_mwt_continue(test):
    _check_continue(test, "gwt_" + test.name, "mwt", "MWT")


def check_output_mwe_continue(test):
    _check_continue(test, "gwe_" + test.name, "mwe", "MWE")


_CHECK_FNS = {
    "mwt_01ts": check_output_mwt_ts,
    "mwt_01paux": check_output_mwt_peraux,
    "mwe_01ts": check_output_mwe_ts,
    "mwe_01paux": check_output_mwe_peraux,
    "mwtauxper": check_output_mwt_auxper,
    "mweauxper": check_output_mwe_auxper,
    "mwtauxovr": check_output_mwt_auxovr,
    "mwtconccont": check_output_mwt_conccont,
    "mweconccont": check_output_mwe_conccont,
    "mwtratecont": check_output_mwt_ratecont,
    "mweratecont": check_output_mwe_ratecont,
    "mwtcontinue": check_output_mwt_continue,
    "mwecontinue": check_output_mwe_continue,
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
