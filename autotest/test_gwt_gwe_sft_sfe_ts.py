"""
SFT/SFE time-series, PERIOD AUXILIARY, and cross-period persistence
integration tests.

Covers GWT-SFT and GWE-SFE in paired cases that share the same GWF+SFR model.

  PERIOD TS: CONCENTRATION/TEMPERATURE (CONSTANT status) and RAINFALL
             (sft_ad_ts/sfe_ad_ts's own package-specific re-sync) linked to
             a STEPWISE TS that steps at t=0.05. The obs file verifies that
             the per-timestep re-sync picks up the new value.

  PACKAGEDATA AUX TS: AUX1 carries TS name "aux1_ts"; the binary budget
             AUXILIARY term verifies per-timestep re-sync.

  PERIOD AUXILIARY literal: PERIOD AUXILIARY overrides PACKAGEDATA baseline
             (99.0 -> 55.0). Binary budget AUXILIARY must show 55.0.

  Cross-period persistence: a PERIOD INFLOW value linked to a TS keeps
             tracking that TS in a later period whose own PERIOD block
             reappears (for STATUS/CONCENTRATION) without repeating INFLOW.
             INFLOW is a package-specific (non-AUX) field, wired in
             gwt-sft.f90/gwe-sfe.f90 -- confirms the generalized
             apply_period_settings mechanism generalizes without
             per-package changes.

Cases:
  sft_01ts     — CONCENTRATION TS + PACKAGEDATA AUX TS for GWT-SFT
  sft_01paux   — PERIOD AUXILIARY literal for GWT-SFT
  sfe_01ts     — TEMPERATURE TS + PACKAGEDATA AUX TS for GWE-SFE
  sfe_01paux   — PERIOD AUXILIARY literal for GWE-SFE
  sftinflcont  — INFLOW TS continues across periods, GWT-SFT
  sfeinflcont  — INFLOW TS continues across periods, GWE-SFE
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "sft_01ts",
    "sft_01paux",
    "sfe_01ts",
    "sfe_01paux",
    "sftinflcont",
    "sfeinflcont",
]

TS_STEP = 0.05  # TS step-change time within the 0.1-day period
CONC_TEMP_EARLY = 25.0  # reach conc/temp (CONSTANT status) for steps 1-5
CONC_TEMP_LATE = 0.0  # reach conc/temp for steps 6-10
AUX1_EARLY = 99.0
AUX1_LATE = 0.0
PKGDATA_AUX1 = 99.0
PERIOD_AUX1 = 55.0
RAIN_EARLY = 25.0  # RAINFALL concentration/temperature, steps 1-5
RAIN_LATE = 0.0  # steps 6-10
Q_RAIN = 0.1  # GWF SFR reach 0 rainfall rate; RAINFALL obs = Q_RAIN * conc

# cross-period persistence: distinct value per period, 3 periods x 5 steps
TS_VAL_P1 = 10.0
TS_VAL_P2 = 20.0
TS_VAL_P3 = 30.0
Q_INFLOW = 1.0  # constant GWF SFR external inflow rate; obs = Q_INFLOW * conciflw

# GWE EST heat-capacity parameters
CPW = 4183.0
RHOW = 999.728
LHV = 2500.0
CPS = 800.0
RHOS = 2650.0
KTF = 1.0
RBTHCND = 0.1


def _grid_params():
    nlay, nrow, ncol = 1, 1, 7
    lx, lz = 7.0, 1.0
    delr = lx / ncol
    delc = 1.0
    delz = lz / nlay
    top = 0.0
    botm = [top - (k + 1) * delz for k in range(nlay)]
    idomain = np.full((nlay, nrow, ncol), 1)
    return nlay, nrow, ncol, delr, delc, top, botm, idomain


def _build_sfr_packagedata(ncol, delr, delc):
    rlen, rwid, rgrd, rtp, rbth, rhk, rman, ustrf, ndv = (
        delr,
        delc,
        1.0,
        0.0,
        0.1,
        0.0,
        1.0,
        1.0,
        0,
    )
    pak_data = []
    for irno in range(ncol):
        ncon = 2 if irno not in [0, ncol - 1] else 1
        pak_data.append(
            (
                irno,
                (0, 0, irno),
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
            )
        )
    con_data = []
    for irno in range(ncol):
        if irno == 0:
            con_data.append((irno, -(irno + 1)))
        elif irno == ncol - 1:
            con_data.append((irno, irno - 1))
        else:
            con_data.append((irno, irno - 1, -(irno + 1)))
    return pak_data, con_data


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
    """Build GWF+SFR model; aux_name is 'CONCENTRATION' or 'TEMPERATURE'."""
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
        stress_period_data=[[(0, 0, ncol - 1), 0.0, 0.0]],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=[[(0, 0, 0), 1.0, 0.0]],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.wel",
    )
    pak_data, con_data = _build_sfr_packagedata(ncol, delr, delc)
    flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        stage_filerecord=f"{gwfname}.sfr.stg",
        budget_filerecord=f"{gwfname}.sfr.bud",
        nreaches=ncol,
        packagedata=pak_data,
        connectiondata=con_data,
        perioddata=[(0, "INFLOW", 1.0), (0, "RAINFALL", Q_RAIN)],
        pname="SFR-1",
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
    """Build GWT model without SFT; returns (gwt, gwtname)."""
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
    flopy.mf6.ModflowGwtmst(gwt, porosity=1.0, filename=f"{gwtname}.sto")
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[("CHD-1", "AUX", "CONCENTRATION"), ("WEL-1", "AUX", "CONCENTRATION")],
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
    """Build GWE model without SFE; returns (gwe, gwename)."""
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
        porosity=1.0,
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
        sources=[("CHD-1", "AUX", "TEMPERATURE"), ("WEL-1", "AUX", "TEMPERATURE")],
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


def _sft_packagedata_ts(ncol):
    return [(irno, 0.0, "aux1_ts", 999.0, f"reach{irno + 1}") for irno in range(ncol)]


def _sft_packagedata_literal(ncol, aux1_val):
    return [(irno, 0.0, aux1_val, 999.0, f"reach{irno + 1}") for irno in range(ncol)]


def _sfe_packagedata_ts(ncol):
    return [
        (irno, 0.0, KTF, RBTHCND, "aux1_ts", 999.0, f"reach{irno + 1}")
        for irno in range(ncol)
    ]


def _sfe_packagedata_literal(ncol, aux1_val):
    return [
        (irno, 0.0, KTF, RBTHCND, aux1_val, 999.0, f"reach{irno + 1}")
        for irno in range(ncol)
    ]


def _build_sfr_packagedata_default(ncol, delr, delc):
    """SFR PACKAGEDATA with a trailing AUX column, for _build_gwf_default's
    auxiliary=[aux_name] SFR package."""
    rlen, rwid, rgrd, rtp, rbth, rhk, rman, ustrf, ndv = (
        delr,
        delc,
        1.0,
        0.0,
        0.1,
        0.0,
        1.0,
        1.0,
        0,
    )
    pak_data = []
    for irno in range(ncol):
        ncon = 2 if irno not in [0, ncol - 1] else 1
        pak_data.append(
            (
                irno,
                (0, 0, irno),
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
                0.0,
            )
        )
    con_data = []
    for irno in range(ncol):
        if irno == 0:
            con_data.append((irno, -(irno + 1)))
        elif irno == ncol - 1:
            con_data.append((irno, irno - 1))
        else:
            con_data.append((irno, irno - 1, -(irno + 1)))
    return pak_data, con_data


def _build_gwf_default(
    sim, name, aux_name, delr, delc, nlay, nrow, ncol, top, botm, idomain
):
    """Build GWF+SFR model with flopy's default IMS settings (no custom
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
        stress_period_data=[[(0, 0, ncol - 1), 0.0, 0.0]],
        pname="CHD-1",
        auxiliary=aux_name,
        filename=f"{gwfname}.chd",
    )
    pak_data, con_data = _build_sfr_packagedata_default(ncol, delr, delc)
    flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        budget_filerecord=f"{gwfname}.sfr.bud",
        nreaches=ncol,
        packagedata=pak_data,
        connectiondata=con_data,
        perioddata=[(0, "INFLOW", Q_INFLOW)],
        pname="SFR-1",
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
    """Build GWT model without SFT, flopy IMS defaults; returns (gwt, gwtname)."""
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
    flopy.mf6.ModflowGwtmst(gwt, porosity=1.0)
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
    """Build GWE model without SFE, flopy IMS defaults; returns (gwe, gwename)."""
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
        porosity=1.0,
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


def _inflow_perioddata(depvarkey):
    return {
        0: [
            (0, "STATUS", "ACTIVE"),
            (0, depvarkey, 0.0),
            (0, "INFLOW", "inflow_ts"),
        ],
        1: [
            (0, "STATUS", "ACTIVE"),
            (0, depvarkey, 0.0),
            # no INFLOW setting -- period 0's TS-linked value should
            # persist and keep tracking the TS
        ],
        2: [(0, "STATUS", "ACTIVE"), (0, depvarkey, 0.0)],
    }


def build_models_sftinflcont(name, ws):
    """Case sftinflcont: PERIOD INFLOW linked to a TS keeps tracking it in
    a later period whose own PERIOD block doesn't mention it, GWT-SFT."""
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
    sft_obs = {(gwtname + ".sft.obs.csv",): [("sft1-iflw", "ext-inflow", "reach1")]}
    sft = flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=[(irno, 0.0, 999.0, f"reach{irno + 1}") for irno in range(ncol)],
        reachperioddata=_inflow_perioddata("CONCENTRATION"),
        observations=sft_obs,
        pname="SFR-1",
        auxiliary=["aux1"],
    )
    sft.ts.initialize(
        filename="inflow.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="inflow_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_sfeinflcont(name, ws):
    """Case sfeinflcont: PERIOD INFLOW continues across periods, GWE-SFE."""
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
    sfe_obs = {(gwename + ".sfe.obs.csv",): [("sfe1-iflw", "ext-inflow", "reach1")]}
    sfe = flopy.mf6.ModflowGwesfe(
        gwe,
        boundnames=True,
        save_flows=True,
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=[
            (irno, 0.0, KTF, RBTHCND, 999.0, f"reach{irno + 1}") for irno in range(ncol)
        ],
        reachperioddata=_inflow_perioddata("TEMPERATURE"),
        observations=sfe_obs,
        flow_package_name="SFR-1",
        flow_package_auxiliary_name="TEMPERATURE",
        pname="SFR-1",
        auxiliary=["aux1"],
    )
    sfe.ts.initialize(
        filename="inflow.ts",
        timeseries=[
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ],
        time_series_namerecord="inflow_ts",
        interpolation_methodrecord="stepwise",
    )
    return sim, None


def build_models_sft_ts(name, ws):
    """Case sft_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWT-SFT."""
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

    sft_obs = {
        (gwtname + ".sft.obs.csv",): [
            ("sft1-conc", "CONCENTRATION", 1),
            ("sft1-rain", "RAINFALL", 1),
        ],
    }
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
        packagedata=_sft_packagedata_ts(ncol),
        reachperioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "CONCENTRATION", "conc_ts"),
            (0, "RAINFALL", "rain_ts"),
        ],
        observations=sft_obs,
        pname="SFR-1",
        auxiliary=["aux1", "aux2"],
    )

    sft.ts.initialize(
        filename=f"{gwtname}.sft.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, RAIN_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, RAIN_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, RAIN_LATE),
        ],
        time_series_namerecord=[("conc_ts", "aux1_ts", "rain_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_sft_peraux(name, ws):
    """Case sft_01paux: PERIOD AUXILIARY literal (T3) for GWT-SFT."""
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

    flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=_sft_packagedata_literal(ncol, PKGDATA_AUX1),
        reachperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "CONCENTRATION", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="SFR-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


def build_models_sfe_ts(name, ws):
    """Case sfe_01ts: PERIOD TS (T1) + PACKAGEDATA AUX TS (T2) for GWE-SFE."""
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

    sfe_obs = {
        (gwename + ".sfe.obs.csv",): [
            ("sfe1-temp", "TEMPERATURE", 1),
            ("sfe1-rain", "RAINFALL", 1),
        ],
    }
    sfe_obs["digits"] = 7
    sfe_obs["print_input"] = True
    sfe_obs["filename"] = gwename + ".sfe.obs"

    sfe = flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=_sfe_packagedata_ts(ncol),
        reachperioddata=[
            (0, "STATUS", "CONSTANT"),
            (0, "TEMPERATURE", "temp_ts"),
            (0, "RAINFALL", "rain_ts"),
        ],
        observations=sfe_obs,
        pname="SFR-1",
        auxiliary=["aux1", "aux2"],
    )

    sfe.ts.initialize(
        filename=f"{gwename}.sfe.ts",
        timeseries=[
            (0.0, CONC_TEMP_EARLY, AUX1_EARLY, RAIN_EARLY),
            (TS_STEP, CONC_TEMP_LATE, AUX1_LATE, RAIN_LATE),
            (0.1, CONC_TEMP_LATE, AUX1_LATE, RAIN_LATE),
        ],
        time_series_namerecord=[("temp_ts", "aux1_ts", "rain_ts")],
        interpolation_methodrecord=[("STEPWISE", "STEPWISE", "STEPWISE")],
    )

    return sim, None


def build_models_sfe_peraux(name, ws):
    """Case sfe_01paux: PERIOD AUXILIARY literal (T3) for GWE-SFE."""
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

    flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=_sfe_packagedata_literal(ncol, PKGDATA_AUX1),
        reachperioddata=[
            (0, "STATUS", "ACTIVE"),
            (0, "TEMPERATURE", 0.0),
            (0, "AUXILIARY", "aux1", PERIOD_AUX1),
        ],
        pname="SFR-1",
        auxiliary=["aux1", "aux2"],
    )

    return sim, None


_BUILD_FNS = {
    "sft_01ts": build_models_sft_ts,
    "sft_01paux": build_models_sft_peraux,
    "sfe_01ts": build_models_sfe_ts,
    "sfe_01paux": build_models_sfe_peraux,
    "sftinflcont": build_models_sftinflcont,
    "sfeinflcont": build_models_sfeinflcont,
}


def build_models(idx, test):
    name = cases[idx]
    return _BUILD_FNS[name](name, test.workspace)


def _check_ts(bud_fname, obs_fname, pkg_label, conc_col, rain_col, thermal=False):
    """Shared T1+T2 check for SFT and SFE. rain_col checks sft_ad_ts/
    sfe_ad_ts's own RAINFALL TS re-sync, not just the base-class one."""
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

    # RAINFALL obs is a mass (or, for GWE, energy) flow rate
    # (Q_RAIN * concentration/temperature), not the concentration itself;
    # GWE flow terms are additionally scaled by RHOW*CPW
    scale = RHOW * CPW if thermal else 1.0
    rain_obs = tc[rain_col]
    rain_early = Q_RAIN * RAIN_EARLY * scale
    assert np.allclose(rain_obs[:5], rain_early, rtol=1e-5), (
        f"{pkg_label} steps 1-5 rainfall expected {rain_early}, got {rain_obs[:5]}"
    )
    assert np.allclose(rain_obs[5:], 0.0, atol=1e-10 * max(scale, 1.0)), (
        f"{pkg_label} steps 6-10 rainfall expected 0.0, got {rain_obs[5:]}"
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
    """Shared T3 check for SFT and SFE."""
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


def check_output_sft_ts(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_ts(
        os.path.join(ws, f"{gwtname}.sft.bud"),
        os.path.join(ws, f"{gwtname}.sft.obs.csv"),
        "SFT",
        "SFT1CONC",
        "SFT1RAIN",
    )


def check_output_sft_peraux(test):
    ws, name = test.workspace, test.name
    gwtname = "gwt_" + name
    _check_peraux(os.path.join(ws, f"{gwtname}.sft.bud"), "SFT")


def check_output_sfe_ts(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_ts(
        os.path.join(ws, f"{gwename}.sfe.bud"),
        os.path.join(ws, f"{gwename}.sfe.obs.csv"),
        "SFE",
        "SFE1TEMP",
        "SFE1RAIN",
        thermal=True,
    )


def check_output_sfe_peraux(test):
    ws, name = test.workspace, test.name
    gwename = "gwe_" + name
    _check_peraux(os.path.join(ws, f"{gwename}.sfe.bud"), "SFE")


def _check_inflcont(test, gwname, ext, budname, thermal):
    fname = test.workspace / f"{gwname}.{ext}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")
    iflw_obs = tc["SFT1IFLW"] if "SFT1IFLW" in tc.dtype.names else tc["SFE1IFLW"]
    assert len(iflw_obs) == 15, (
        f"Expected 15 obs records (5/period, 3 periods), got {len(iflw_obs)}"
    )
    period1, period2, period3 = iflw_obs[:5], iflw_obs[5:10], iflw_obs[10:]
    scale = CPW * RHOW if thermal else 1.0
    exp1 = Q_INFLOW * TS_VAL_P1 * scale
    exp2 = Q_INFLOW * TS_VAL_P2 * scale
    exp3 = Q_INFLOW * TS_VAL_P3 * scale
    assert np.allclose(period1, exp1), (
        f"{budname} period 1 EXT-INFLOW obs expected {exp1}, got {period1}"
    )
    assert np.allclose(period2, exp2), (
        f"{budname} period 2 EXT-INFLOW obs expected {exp2} (TS should "
        f"still be tracked even though period 1's PERIOD block doesn't "
        f"mention INFLOW), got {period2}"
    )
    assert np.allclose(period3, exp3), (
        f"{budname} period 3 EXT-INFLOW obs expected {exp3}, got {period3}"
    )


def check_output_sftinflcont(test):
    _check_inflcont(test, "gwt_" + test.name, "sft", "SFT", thermal=False)


def check_output_sfeinflcont(test):
    _check_inflcont(test, "gwe_" + test.name, "sfe", "SFE", thermal=True)


_CHECK_FNS = {
    "sft_01ts": check_output_sft_ts,
    "sft_01paux": check_output_sft_peraux,
    "sfe_01ts": check_output_sfe_ts,
    "sfe_01paux": check_output_sfe_peraux,
    "sftinflcont": check_output_sftinflcont,
    "sfeinflcont": check_output_sfeinflcont,
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
