# Test the use of the sensible heat flux utility used in conjunction with the
# SFE advanced package.  This test is a single cell with a single reach.
# Channel flow characteristics are unrealistic: Manning's n is unrealistically
# low and slope is extremely high. These conditions result in an extremely high
# streamflow velocity that results in nearly all of the heat being added to the
# channel exiting at the outlet with very near negligle heat storage increases
# in the channel.  The result is a 1 deg C rise in temperature in the
# streamflow - an easy result to confirm in this test.

import os

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["sfe-shf"]

# Model units
length_units = "m"
time_units = "seconds"

# model domain and grid definition

nrow = 1
ncol = 1
nlay = 1
delr = 1.0
delc = 1.0
xmax = ncol * delr
ymax = nrow * delc

ibound = 1
top = 1.0
botm = 0.0
strthd = 0.0
chd_on = False

# model input parameters
k11 = 500.0
strt_gw_temp = 99.0

ss = 0.00001
sy = 0.20
hani = 1
laytyp = 1

# Package boundary conditions
sfr_evaprate = 0.0
rhk = 0.0
rwid = 1.0
strm_temp = 1.0
surf_Q_in = [
    [10.0],
]
# sensible heat flux parameter values
wspd = 20.0
tatm = 1189766.7  # unrealistically high to drive a 1 deg C rise in stream temperature


# Transport related parameters
porosity = sy  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rhow = 1000  # Density of water ($kg/m^3$)
rhos = 2650  # Density of the aquifer material ($kg/m^3$)
rhoa = 1.225  # Density of the atmosphere ($kg/m^3$)
Cpw = 4180  # Heat capacity of water ($J/kg/C$)
Cps = 880  # Heat capacity of the solids ($J/kg/C$)
Cpa = 717.0  # Heat capacity of the atmosphere ($J/kg/C$)
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)
c_d = 0.002  # Drag coefficient ($unitless$)
# Thermal conductivity of the streambed material ($W/m/C$)
K_therm_strmbed = 0.0
rbthcnd = 0.0001

# time params
steady = {0: True, 1: False}
transient = {0: False, 1: True}
nstp = [1]
tsmult = [1]
perlen = [1]

nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-3, 1e-4, 0.97

#
# MODFLOW 6 flopy GWF object
#


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    print(f"Building model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating time discretization
    tdis_rc = []
    for i in range(len(nstp)):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(
        sim, nper=len(nstp), perioddata=tdis_rc, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="newton",
    )

    # Instantiating solver
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="cooley",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])

    # Instantiate discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Instantiate node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
    )

    # Instantiate storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state=steady,
        transient=transient,
    )

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strthd)

    # Instantiate output control package
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate streamflow routing package
    # Determine the middle row and store in rMid (account for 0-base)
    rMid = 1
    # sfr data
    nreaches = ncol
    rlen = delr
    roughness = 1e-10
    rbth = 0.1
    strmbd_hk = rhk
    strm_up = 0.95
    strm_dn = 0.94
    # divide by 10 to further reduce slop
    slope = 0.04
    nconn = 0
    ustrf = 1.0
    ndv = 0
    strm_incision = 0.05

    packagedata = []
    for irch in range(nreaches):
        rp = [
            irch,
            (0, 0, 0),
            rlen,
            rwid,
            slope,
            top - strm_incision,
            rbth,
            strmbd_hk,
            roughness,
            nconn,
            ustrf,
            ndv,
        ]
        packagedata.append(rp)

    connectiondata = [0]

    sfr_perioddata = {}
    for t in np.arange(len(surf_Q_in[idx])):
        sfrbndx = []
        for i in np.arange(nreaches):
            if i == 0:
                sfrbndx.append([i, "INFLOW", surf_Q_in[idx][t]])
            # sfrbndx.append([i, "EVAPORATION", sfr_evaprate])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname}.sfr.obs.csv": [
            ("rch1_depth", "depth", 1),
            ("rch1_outf", "ext-outflow", 1),
            ("rch1_wetwidth", "wet-width", 1),
        ],
        "digits": 8,
        "print_input": True,
        "filename": gwfname + ".sfr.obs",
    }

    budpth = f"{gwfname}.sfr.cbc"
    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=1.0,
        budget_filerecord=budpth,
        mover=False,
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=sfr_perioddata,
        observations=sfr_obs,
        pname="SFR",
        filename=f"{gwfname}.sfr",
    )

    # --------------------------------------------------
    # Setup the GWE model for simulating heat transport
    # --------------------------------------------------
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename)

    # Instantiating solver for GWT
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

    # Instantiating DIS for GWE
    flopy.mf6.ModflowGwedis(
        gwe,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        pname="DIS",
        filename=f"{gwename}.dis",
    )

    # Instantiate Mobile Storage and Transfer package
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=porosity,
        heat_capacity_water=Cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=Cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # Instantiate Energy Transport Initial Conditions package
    flopy.mf6.ModflowGweic(gwe, strt=strt_gw_temp)

    # Instantiate Advection package
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")

    # Instantiate Dispersion package (also handles conduction)
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        ktw=0.5918,
        kts=0.2700,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    # [b/c at least one boundary back is active (SFR), ssm must be on]
    sourcerecarray = [[]]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{gwename}.ssm")

    # Instantiate Streamflow Energy Transport package
    sfepackagedata = []
    for irno in range(ncol):
        t = (irno, strm_temp, K_therm_strmbed, rbthcnd)
        sfepackagedata.append(t)

    sfeperioddata = []
    for irno in range(ncol):
        if irno == 0:
            sfeperioddata.append((irno, "INFLOW", strm_temp))

    # Instantiate SFE observation points
    sfe_obs = {
        f"{gwename}.sfe.obs.csv": [
            ("rch1_outftemp", "temperature", 1),
            ("rch1_outfener", "ext-outflow", 1),
        ],
        "digits": 8,
        "print_input": True,
        "filename": gwename + ".sfe.obs",
    }

    shf_filename = f"{gwename}.sfe.shf"
    sfe = flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=False,
        save_flows=True,
        print_input=False,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=sfepackagedata,
        reachperioddata=sfeperioddata,
        flow_package_name="SFR",
        observations=sfe_obs,
        pname="SFE",
        filename=f"{gwename}.sfe",
    )

    # Shf utility
    shf_spd = {}
    for kper in range(len(nstp)):
        spd = []
        for irno in range(ncol):
            spd.append([irno, "WSPD", wspd])
            spd.append([irno, "TATM", tatm])
        shf_spd[kper] = spd

    shf = flopy.mf6.ModflowUtlshf(
        sfe,
        print_input=True,
        density_air=rhoa,
        heat_capacity_air=Cpa,
        drag_coefficient=c_d,
        reachperioddata=shf_spd,
        filename=shf_filename,
    )

    # Instantiate Output Control package for transport
    flopy.mf6.ModflowGweoc(
        gwe,
        temperature_filerecord=f"{gwename}.ucn",
        saverecord=[("TEMPERATURE", "ALL")],
        temperatureprintrecord=[("COLUMNS", 3, "WIDTH", 20, "DIGITS", 8, "GENERAL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwename}.oc",
    )

    # Instantiate Gwf-Gwe Exchange package
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{gwename}.gwfgwe",
    )

    return sim, None


# sim, dum = build_models(0, r"c:\temp\_shf00")
# sim.write_simulation()


def check_output(idx, test):
    print("evaluating results...")
    msg0 = "Stream channel width less than 1.0, should be 1.0 m"

    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    # calc expected rise in temperature independent of mf6

    fpth = os.path.join(test.workspace, gwfname + ".sfr.obs.csv")
    assert os.path.isfile(fpth)
    df = pd.read_csv(fpth)
    calc_strm_wid = df.loc[0, "RCH1_WETWIDTH"].copy()
    # confirm stream width is 1.0 m
    assert np.isclose(calc_strm_wid, 1.0, atol=1e-9), msg0

    # confirm that the energy added to the stream results in a 1 deg C rise in temp
    # temperature gradient
    tgrad = tatm - strm_temp
    ener_per_sqm = c_d * rhoa * Cpa * wspd * tgrad
    ener_transfer = ener_per_sqm * (delr * calc_strm_wid)
    # calculate expected temperature rise based on energy transfer
    temp_rise = ener_transfer / (surf_Q_in[idx][0] * Cpw * rhow)

    fpth2 = os.path.join(test.workspace, gwename + ".sfe.obs.csv")
    assert os.path.isfile(fpth2)
    df2 = pd.read_csv(fpth2)


    # confirm 1 deg C rise in temp
    msg1 = (
        "The MF6 simulated rise in river temperature does not match \
        external calculations.  The calculated temperature rise \
        is: "
        + str(temp_rise)
    )

    assert np.isclose(df2.loc[0, "RCH1_OUTFTEMP"], strm_temp + temp_rise, atol=1e-6), (
        msg1
    )


# - No need to change any code below
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
