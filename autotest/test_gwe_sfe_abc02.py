# Test the use of the atmospheric boundary condition utility used in conjunction
# with the SFE advanced package.  This test uses four cells that each host a
# single reach. Channel flow characteristics are unrealistic: Manning's n is
# unrealistically low and slope is extremely high. These conditions result in
# an extremely high streamflow velocity that results in nearly all of the heat
# being added to, or subtracted from, the channel at the outlet with near-
# negligle heat storage increases (or decreases) in the channel.
#
# This test applies each of the atmospheric boundary conditions heat fluxes to
# only 1 of the reaches.  In other words, the shortwave radiation is applied to
# the first reach, longwave radiation to the second reach, sensible heat flux
# to the third, and latent heat flux to the fourth reach.  The idea is that the
# temperature change at the outlet of each reach, before it flows into the next
# downstream reach, is equal to externally calculated changes.

"""

     
               \
                \ Reach 1 (swr)
                 \
                  \  
   Reach 2 (lwr)   v   Reach 4 (precip)    Reach 5 (lhf)
            ------>+---------------------+--------------->
                   ^
                  /
                 /
                / Reach 3 (shf)
               /
"""

import math
import os

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["sfe-abc"]

DCTOK = 273.16

# Model units
length_units = "m"
time_units = "seconds"

# model domain and grid definition

nrow = 1
ncol = 1
nlay = 1
delr = 100.0
delc = 100.0
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
rlen = 1.0
# strm_temp = 11.0
strm_in_init = [11.0, 21.867108141, 25.0, 19.3310145906107, 19.27964737200225]
surf_Q_in = 10.0

# sensible and latent heat flux parameter values
wspd = [
    0.0,
    0.0,
    3919295.19947608,
    0.0,
    16118322.9664089,
]  # unrealistically high to drive observable temperature change
patm = [0.0, 0.0, 954.680843658077, 0.0, 954.680843658077]
tatm = [11.0, 40.9752, 11.0, 11.0, 1.0]  # used by lwr, shf, lhf
tatm = [t + DCTOK for t in tatm]
# shortwave radiation parameter values

# unrealistically high to drive a 0.1 deg C rise in stream temperature
solr = [
    4180000.0,
    0.0,
    0.0,
    0.0,
    0.0,
]
shd = [0.0, 0.0, 0.0, 0.0, 0.0]  # 100% shade "turns off" solar flux
swrefl = [0.0, 1.0, 0.0, 1.0, 0.0]
rh = [0.0, 20.0, 1.0, 0.0, 0.01]  # percent
lwrefl = 0.03  # Fogg et al 2023

# atmosphere composition adjustment factor (using dummy value to drive
# half a degree change)
atmc = [0.0, 9667.121567, 0.0, 0.0, 0.0]

# latent heat flux parameter values (these values are specified in the options
# block and are therefore constant across reaches
c_d = 0.0  # Drag coefficient ($unitless$)
wf_slope = 1.383e-08  # wind function slope ($1/mbar$)
wf_int = 3.445e-09  # wind function intercept ($m/s$)
emiss_water = 0.95  # Fogg et al 2023
emiss_riparian = 0.97  # Fogg et al 2023
stephan_boltzmann = 5.670374419e-08

# Transport related parameters
porosity = sy  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rhow = 1000.0  # Density of water ($kg/m^3$)
rhos = 2650.0  # Density of the aquifer material ($kg/m^3$)
rhoa = 1.225  # Density of the atmosphere ($kg/m^3$)
Cpw = 4180.0  # Heat capacity of water ($J/kg/C$)
Cps = 880.0  # Heat capacity of the solids ($J/kg/C$)
Cpa = 717.0  # Heat capacity of the atmosphere ($J/kg/C$)
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)

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
hclose, rclose, relax = 1e-10, 1e-10, 0.97

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

    # explicitly set connections
    conns = [(0, -3), (1, -3), (2, -3), (3, 0, 1, 2, -4), (4, 3)]
    nreaches = len(conns)

    packagedata = []
    for irch in range(nreaches):
        ncon = len(conns[irch]) - 1
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
            ncon,
            ustrf,
            ndv,
        ]
        packagedata.append(rp)

    sfr_perioddata = {}
    for t in np.arange(len(perlen)):
        sfrbndx = []
        for i in np.arange(nreaches):
            # only specify inflow for the first three reaches
            if i < 3:
                sfrbndx.append([i, "INFLOW", surf_Q_in])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname}.sfr.obs.csv": [
            ("rch1_depth", "depth", 1),
            ("rch1_outf", "ext-outflow", 1),
            ("rch1_wetwidth", "wet-width", 1),
            ("rch2_wetwidth", "wet-width", 2),
            ("rch3_wetwidth", "wet-width", 3),
            ("rch4_wetwidth", "wet-width", 4),
            ("rch5_wetwidth", "wet-width", 5),
            ("rch1_stg", "stage", 1),
            ("rch2_stg", "stage", 2),
            ("rch3_stg", "stage", 3),
            ("rch4_stg", "stage", 4),
            ("rch5_stg", "stage", 5),
        ],
        "digits": 15,
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
        connectiondata=conns,
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
        rcloserecord=f"{rclose} strict",
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
    for irno in range(len(conns)):
        t = (irno, strm_in_init[irno], K_therm_strmbed, rbthcnd)
        sfepackagedata.append(t)

    sfeperioddata = []
    sfeperioddata.append((0, "INFLOW", strm_in_init[0]))
    sfeperioddata.append((1, "INFLOW", strm_in_init[1]))
    sfeperioddata.append((2, "INFLOW", strm_in_init[2]))

    # Instantiate SFE observation points
    sfe_obs = {
        f"{gwename}.sfe.obs.csv": [
            ("rch5_evap", "surfevap", 5),
            ("rch1_outftemp", "temperature", 1),
            ("rch2_outftemp", "temperature", 2),
            ("rch3_outftemp", "temperature", 3),
            ("rch4_outftemp", "temperature", 4),
            ("rch5_outftemp", "temperature", 5),
            ("rch1_outfener", "ext-outflow", 1),
            ("rch2_outfener", "ext-outflow", 2),
            ("rch3_outfener", "ext-outflow", 3),
            ("rch4_outfener", "ext-outflow", 4),
            ("rch5_outfener", "ext-outflow", 5),
            ("rch1_swr", "swr", 1),
            ("rch2_swr", "swr", 2),
            ("rch3_swr", "swr", 3),
            ("rch4_swr", "swr", 4),
            ("rch5_swr", "swr", 5),
            ("rch1_lwr", "lwr", 1),
            ("rch2_lwr", "lwr", 2),
            ("rch3_lwr", "lwr", 3),
            ("rch4_lwr", "lwr", 4),
            ("rch5_lwr", "lwr", 5),
            ("rch1_lhf", "lhf", 1),
            ("rch2_lhf", "lhf", 2),
            ("rch3_lhf", "lhf", 3),
            ("rch4_lhf", "lhf", 4),
            ("rch5_lhf", "lhf", 5),
            ("rch1_shf", "shf", 1),
            ("rch2_shf", "shf", 2),
            ("rch3_shf", "shf", 3),
            ("rch4_shf", "shf", 4),
            ("rch5_shf", "shf", 5),
        ],
        "digits": 20,
        "print_input": True,
        "filename": gwename + ".sfe.obs",
    }

    abc_filename = f"{gwename}.sfe.abc"
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

    # Abc utility
    abc_spd = {}
    for kper in range(len(nstp)):
        spd = []
        for irno in range(len(conns)):
            spd.append([irno, "WSPD", wspd[irno]])
            spd.append([irno, "TATM", tatm[irno]])
            spd.append([irno, "SOLR", solr[irno]])
            spd.append([irno, "SHD", shd[irno]])
            spd.append([irno, "SWREFL", swrefl[irno]])
            spd.append([irno, "RH", rh[irno]])
            spd.append([irno, "ATMC", atmc[irno]])
            spd.append([irno, "PATM", patm[irno]])

        abc_spd[kper] = spd

    abc = flopy.mf6.ModflowUtlabc(
        sfe,
        print_input=True,
        density_air=rhoa,
        heat_capacity_air=Cpa,
        drag_coefficient=c_d,
        emissivity_water=emiss_water,
        emissivity_canopy=emiss_riparian,
        wind_func_slope=wf_slope,
        wind_func_int=wf_int,
        longwave_reflectance=lwrefl,
        reachperioddata=abc_spd,
        filename=abc_filename,
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


def calc_lwr_ener_transfer(ifno, updated_strm_temp):
    Ql_up = emiss_water * stephan_boltzmann * ((updated_strm_temp + DCTOK) ** 4)
    e_s = 6.1275 * math.exp(17.2693882 * ((tatm[ifno] - DCTOK) / (tatm[ifno] - 35.86)))
    e_a = (rh[ifno] / 100.0) * e_s
    emiss_air = (1.24 * (e_a / tatm[ifno]) ** (1.0 / 7.0)) * atmc[ifno]
    emiss_down = (1.0 - shd[ifno]) * emiss_air + shd[
        ifno
    ] * emiss_riparian  # calcs for emiss_riparian

    Ql_dn = emiss_down * stephan_boltzmann * (tatm[ifno] ** 4)
    # calc net longwave
    lwrflx = Ql_dn * (1.0 - lwrefl) - Ql_up

    return lwrflx


def calc_lhf_ener_transfer(ifno, updated_strm_temp, mf_strm_wid):
    L = 2499.64 - (2.51 * (updated_strm_temp - DCTOK))
    e_w = 6.1275 * math.exp(
        17.2693882 * ((updated_strm_temp - DCTOK) / (updated_strm_temp - 35.86))
    )
    e_s = 6.1275 * math.exp(17.2693882 * ((tatm[ifno] - DCTOK) / (tatm[ifno] - 35.86)))
    e_a = (rh[ifno] / 100) * e_s
    vap_press_deficit = e_w - e_a
    wind_function = wf_int + wf_slope * wspd[ifno]
    Ev = wind_function * vap_press_deficit
    lhf_ener_per_sqm = Ev * L * rhow

    ener_transfer = lhf_ener_per_sqm * rlen * mf_strm_wid

    return ener_transfer


def calc_shf_ener_transfer(ifno, lhflx, updated_strm_temp):
    e_w = 6.1275 * math.exp(
        17.2693882 * ((updated_strm_temp - DCTOK) / (updated_strm_temp - 35.86))
    )
    e_s = 6.1275 * math.exp(17.2693882 * ((tatm[ifno] - DCTOK) / (tatm[ifno] - 35.86)))
    e_a = (rh[ifno] / 100) * e_s
    br = 0.00061 * patm[ifno] * (updated_strm_temp - tatm[ifno]) / (e_w - e_a)

    shflx = br * lhflx

    return shflx


def check_output(idx, test):
    print("evaluating results...")
    msg0 = "Stream channel width less than 1.0, should be 1.0 m"
    msg1 = "SFE reach 1 should have roughly a 0.1 deg C rise from short wave"
    msg2 = "SFE reach 2 should have roughly a 0.1 deg C rise from long wave"
    msg3 = "SFE reach 3 not as expected"
    msg4 = "SFE reach 4 not reflective of simple mixing of reaches 1-3"
    msg5 = "SFE reach 5 evaporation rate from latent heat flux not 5 as expected"

    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    # calc expected rise in temperature independent of mf6

    fpth = os.path.join(test.workspace, gwfname + ".sfr.obs.csv")
    assert os.path.isfile(fpth)
    df = pd.read_csv(fpth)
    mf_strm_wid = df.loc[0, "RCH1_WETWIDTH"].copy()
    # confirm stream width is 1.0 m
    assert np.isclose(mf_strm_wid, 1.0, atol=1e-9), msg0

    # confirm that the stream temperature change fits with expectations
    # reach 1 - solar radiation
    fpth2 = os.path.join(test.workspace, gwename + ".sfe.obs.csv")
    assert os.path.isfile(fpth2)
    df2 = pd.read_csv(fpth2)
    # input such that a roughly 0.1 degC rise should result
    ifno = 0
    swr_ener_input = (
        (1.0 - shd[ifno]) * (1 - swrefl[ifno]) * solr[ifno] * mf_strm_wid * rlen
    )
    rch1_Tdelt = swr_ener_input / (Cpw * rhow * surf_Q_in)
    rch1Temp = df2.loc[0, "RCH1_OUTFTEMP"]
    sim_Tdelt = rch1Temp - strm_in_init[ifno]
    # atol needed to account for a small difference associated with a small
    # amount of energy storage change in the reach
    assert np.isclose(rch1_Tdelt, sim_Tdelt, atol=1e-4), msg1

    # reach 2 - longwave radiation
    chng = 1
    ifno = 1  # which is really 2 in 0-based!
    strt_strm_temp = strm_in_init[ifno]
    updated_strm_temp = strm_in_init[ifno]
    mf_strm_wid = df.loc[0, "RCH2_WETWIDTH"].copy()
    while chng > hclose:
        lwr_ener_transfer = calc_lwr_ener_transfer(ifno, updated_strm_temp)
        temp_change = lwr_ener_transfer / (surf_Q_in * Cpw * rhow)
        updated_temp = strm_in_init[ifno] + temp_change
        chng = abs(updated_strm_temp - updated_temp)
        updated_strm_temp = updated_temp

    rch2Temp = df2.loc[0, "RCH2_OUTFTEMP"].copy()
    sim_Tdelt = rch2Temp - strm_in_init[ifno]
    rch2_Tdelt = updated_strm_temp - strm_in_init[ifno]
    assert np.isclose(rch2_Tdelt, sim_Tdelt, atol=1e-4), msg2

    # reach 3 - a focus on sensible heat flux (though other fluxes present)
    chng = 1
    ifno = 2  # which is really 3 in 0-based!
    strt_strm_temp = strm_in_init[ifno]
    updated_strm_temp = strm_in_init[ifno]
    mf_strm_wid = df.loc[0, "RCH3_WETWIDTH"].copy()
    # latent calcs come first
    while chng > hclose:
        lwr_ener_transfer = calc_lwr_ener_transfer(ifno, updated_strm_temp)
        lhf_ener_transfer = calc_lhf_ener_transfer(
            ifno, updated_strm_temp + DCTOK, mf_strm_wid
        )
        shf_ener_transfer = calc_shf_ener_transfer(
            ifno, lhf_ener_transfer, updated_strm_temp + DCTOK
        )
        temp_change = (lwr_ener_transfer - lhf_ener_transfer + shf_ener_transfer) / (
            surf_Q_in * Cpw * rhow
        )
        updated_temp = strt_strm_temp + temp_change
        chng = abs(updated_strm_temp - updated_temp)
        updated_strm_temp = updated_temp

    rch3Temp = df2.loc[0, "RCH3_OUTFTEMP"].copy()
    assert np.isclose(updated_strm_temp, rch3Temp), msg3

    # reach 4 - simple mixing of reaches 1, 2, & 3, ABC boundary fluxes kept to a min
    rch4Temp = df2.loc[0, "RCH4_OUTFTEMP"].copy()
    rch4_expected_Temp = (
        rch1Temp * surf_Q_in + rch2Temp * surf_Q_in + rch3Temp * surf_Q_in
    ) / (3 * surf_Q_in)
    assert np.isclose(rch4Temp, rch4_expected_Temp, atol=1e-5), msg4

    # reach 5 - a focus on latent heat flux (though other fluxes present)
    # an expected evaporation rate of 5 mm/day
    rch5_simEvap = df2.loc[0, "RCH5_EVAP"].copy()
    assert np.isclose(rch5_simEvap, 5.0, atol=1e-8), msg5


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
