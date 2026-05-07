# Test the use of the atmospheric boundary condition utility used in conjunction
# with the SFE advanced package.  This test is focused on the use of time series
# in setting up ABC input.  There is no "check_output" routine, the test is 
# simply that the time series input is read in and the model runs to completion.


import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sfe-abc-ts"]

DCTOK = 273.16

# Model units
length_units = "m"
time_units = "seconds"

# time params
transient = {0: True, 1: True, 2: True, 3: True}
nstp = [1, 1, 1, 1]
tsmult = [1, 1, 1, 1]
perlen = [1, 1, 1, 1]

# model domain and grid definition
nrow = 1
ncol = 5
nlay = 1
delr = 10.0
delc = 10.0

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
strm_temp = 11.0
surf_Q_in = 10.0
# ABC parameter values
# wind speed
wspd = np.random.uniform(0.5, 5.0, size=[5, 5])
wspd_cols = ["wnd_rch1", "wnd_rch2", "wnd_rch3", "wnd_rch4", "wnd_rch5"]
wspd_data = []
for t in range(len(nstp) + 1):
    wspd_data.append((t, wspd[t, 0], wspd[t, 1], wspd[t, 2], wspd[t, 3], wspd[t, 4]))

# atmospheric temperature
tatm = np.random.uniform(10.0, 15.0, size=[5, 5])
tatmK = tatm + DCTOK
tatm_cols = ["tatm_rch1", "tatm_rch2", "tatm_rch3", "tatm_rch4", "tatm_rch5"]
tatm_data = []
for t in range(len(nstp) + 1):
    tatm_data.append(
        (t, tatmK[t, 0], tatmK[t, 1], tatmK[t, 2], tatmK[t, 3], tatmK[t, 4])
    )

# shortwave radiation parameter values
solr = np.random.uniform(800.0, 1200.0, size=[5, 5])
solr_cols = ["solr_rch1", "solr_rch2", "solr_rch3", "solr_rch4", "solr_rch5"]
solr_data = []
for t in range(len(nstp) + 1):
    solr_data.append((t, solr[t, 0], solr[t, 1], solr[t, 2], solr[t, 3], solr[t, 4]))

# shading from solar radiation
shd = np.random.uniform(0.3, 0.7, size=[5, 5])
shd_cols = ["shd_rch1", "shd_rch2", "shd_rch3", "shd_rch4", "shd_rch5"]
shd_data = []
for t in range(len(nstp) + 1):
    shd_data.append((t, shd[t, 0], shd[t, 1], shd[t, 2], shd[t, 3], shd[t, 4]))

# surface water reflectance
swrefl = np.random.uniform(0.03, 0.15, size=[5, 5])
swrefl_cols = ["swr_rch1", "swr_rch2", "swr_rch3", "swr_rch4", "swr_rch5"]
swrefl_data = []
for t in range(len(nstp) + 1):
    swrefl_data.append(
        (t, swrefl[t, 0], swrefl[t, 1], swrefl[t, 2], swrefl[t, 3], swrefl[t, 4])
    )

# relative humidity
rh = np.random.uniform(20.0, 60.0, size=[5, 5])
rh_cols = ["rh_rch1", "rh_rch2", "rh_rch3", "rh_rch4", "rh_rch5"]
rh_data = []
for t in range(len(nstp) + 1):
    rh_data.append((t, rh[t, 0], rh[t, 1], rh[t, 2], rh[t, 3], rh[t, 4]))

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
c_d = 0.0  # Drag coefficient ($unitless$) !!
wf_slope = 1.383e-08  # wind function slope ($1/mbar$)
wf_int = 3.445e-09  # wind function intercept ($m/s$)
# Thermal conductivity of the streambed material ($W/m/C$)
K_therm_strmbed = 0.0
rbthcnd = 0.0001


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
    roughness = 0.03
    rbth = 0.1
    strmbd_hk = rhk
    strm_up = [up for up in np.arange(0.95, 0.91, -0.01)]
    strm_dn = [dwn for dwn in np.arange(0.94, 0.90, -0.01)]
    # divide by 10 to further reduce slop
    slope = 0.001
    ustrf = 1.0
    ndv = 0
    strm_incision = 0.05

    packagedata = []
    for irch in range(nreaches):
        nconn = 1
        if 0 < irch < nreaches - 1:
            nconn += 1
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

    connectiondata = []
    for irch in range(nreaches):
        rc = [irch]
        if irch > 0:
            rc.append(irch - 1)
        if irch < nreaches - 1:
            rc.append(-(irch + 1))
        connectiondata.append(rc)

    sfrbndx = [0, "INFLOW", surf_Q_in]
    sfr_perioddata = {0: sfrbndx}

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

    # Instantiating solver for GWE
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

    # Instantiate Energy Storage and Transfer package
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

    # Instantiate Conduction package 
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
        "digits": 12,
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
    # use time series to specify boundary conditions
    abc_spd = []
    abc_spd.append([0, "WSPD", "wnd_rch1"])
    abc_spd.append([0, "TATM", "tatm_rch1"])
    abc_spd.append([0, "SOLR", "solr_rch1"])
    abc_spd.append([0, "SHD", "shd_rch1"])
    abc_spd.append([0, "SWREFL", "swr_rch1"])
    abc_spd.append([0, "RH", "rh_rch1"])
    abc_spd.append([2, "WSPD", "wnd_rch3"])
    abc_spd.append([2, "TATM", "tatm_rch3"])
    abc_spd.append([2, "SOLR", "solr_rch3"])
    abc_spd.append([2, "SHD", "shd_rch3"])
    abc_spd.append([2, "SWREFL", "swr_rch3"])
    abc_spd.append([2, "RH", "rh_rch3"])
    abc_spd.append([1, "WSPD", "wnd_rch2"])
    abc_spd.append([1, "TATM", "tatm_rch2"])
    abc_spd.append([1, "SOLR", "solr_rch2"])
    abc_spd.append([1, "SHD", "shd_rch2"])
    abc_spd.append([1, "SWREFL", "swr_rch2"])
    abc_spd.append([1, "RH", "rh_rch2"])
    abc_spd.append([3, "WSPD", "wnd_rch4"])
    abc_spd.append([3, "TATM", "tatm_rch4"])
    abc_spd.append([3, "SOLR", "solr_rch4"])
    abc_spd.append([3, "SHD", "shd_rch4"])
    abc_spd.append([3, "SWREFL", "swr_rch4"])
    abc_spd.append([3, "RH", "rh_rch4"])
    abc_spd.append([4, "WSPD", "wnd_rch5"])
    abc_spd.append([4, "TATM", "tatm_rch5"])
    abc_spd.append([4, "SOLR", "solr_rch5"])
    abc_spd.append([4, "SHD", "shd_rch5"])
    abc_spd.append([4, "SWREFL", "swr_rch5"])
    abc_spd.append([4, "RH", "rh_rch5"])

    abc = flopy.mf6.ModflowUtlabc(
        sfe,
        print_input=True,
        density_air=rhoa,
        heat_capacity_air=Cpa,
        drag_coefficient=c_d,
        wind_func_slope=wf_slope,
        wind_func_int=wf_int,
        reachperioddata=abc_spd,
        filename=abc_filename,
    )

    fname1 = f"{gwename}.sfe.abc.wspd.ts"
    abc.ts.initialize(
        filename=fname1,
        timeseries=wspd_data,  # wspd_df.reset_index().to_numpy(),
        time_series_namerecord=wspd_cols,  # wspd_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(wspd_cols),
    )

    fname2 = f"{gwename}.sfe.abc.tatm.ts"
    abc.ts.append_package(
        filename=fname2,
        timeseries=tatm_data,  # tatm_df.reset_index().to_numpy(),
        time_series_namerecord=tatm_cols,  # tatm_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(tatm_cols),
    )

    fname3 = f"{gwename}.sfe.abc.solr.ts"
    abc.ts.append_package(
        filename=fname3,
        timeseries=solr_data,  # solr_df.reset_index().to_numpy(),
        time_series_namerecord=solr_cols,  # solr_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(solr_cols),
    )

    fname4 = f"{gwename}.sfe.abc.shd.ts"
    abc.ts.append_package(
        filename=fname4,
        timeseries=shd_data,  # shd_df.reset_index().to_numpy(),
        time_series_namerecord=shd_cols,  # shd_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(shd_cols),
    )

    fname5 = f"{gwename}.sfe.abc.swrefl.ts"
    abc.ts.append_package(
        filename=fname5,
        timeseries=swrefl_data,  # swrefl_df.reset_index().to_numpy(),
        time_series_namerecord=swrefl_cols,  # swrefl_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(swrefl_cols),
    )

    fname6 = f"{gwename}.sfe.abc.rh.ts"
    abc.ts.append_package(
        filename=fname6,
        timeseries=rh_data,  # rh_df.reset_index().to_numpy(),
        time_series_namerecord=rh_cols,  # rh_df.columns.tolist(),
        interpolation_methodrecord=["linearend"] * len(rh_cols),
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
        xfail=[False],
    )
    test.run()
