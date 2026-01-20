"""
Test problem for GWE: Specifically, this tests conduction between an APT feature
hosted in a drycell with and without NEWTON activated.  The MODFLOW 6 simulation
includes 2 GWF and 2 GWE models.  Neither the flow nor the energy transport
models are connected.  However, there are two GWF-GWE exchanges.  Simulated cell
temperatures for the GWE models are expected to be the same regardless of
whether or not the Newton option is activate in the GWF name file OPTIONS block.

 Model configuration:

                             // <- SFR/SFE channel
                       +----//----+
                      /   ~//~   /|
 ~: Represents       /   ~//~   / |  <- Dry cell
    conduction      /   ~//~   /  |
                   +----//----+   +
                  /   ~//~   /|  /|
                 /   ~//~   / | / |
                /   ~//~   /  |/  |
               +----//----+   +   +
   Dry cell -> |   //     |  /|  /
               |          | / | /
               |          |/  |/
               +----------+   +
 water table   |          |  /
           |   |          | /
          -----+----------+---
               +----------+



"""

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework


# Monotonicity function
def isMonotonic(A):
    x, y = [], []
    x.extend(A)
    y.extend(A)
    x.sort()
    y.sort(reverse=True)
    if np.all(x == A) or np.all(y == A):
        return True
    return False


# Base simulation and model name and workspace

scheme = "UPSTREAM"
# scheme = "TVD"

cases = ["sfecnd"]


# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison

nrow = 2
ncol = 1
nlay = 2
top = 2
bot = np.array([[[1.0], [1.0]], [[0.0], [0.0]]], dtype=float)
strthd = 0.1  # Starting head ($m$)
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)
ss = 1e-6  # Specific storage
sy = 0.20  # Specific Yield
prsity = 0.20  # Porosity
nper = 4  # Number of periods
perlen = [1, 1000, 1, 1000]  # Simulation time ($days$)
nstp = [1, 10, 1, 10]  # 10 day transient time steps
steady = {0: False}
transient = {0: True}

# sfr data
nreaches = 2
rlen = 1.0
rwid = 0.1
roughness = 0.01
rbth = 0.1
rhk = 0.0
slope = 0.001
ustrf = 1.0
ndv = 0
nconn = 1
sfr_packagedata = []
for irno in range(nreaches):
    rp = [
        irno,
        (0, irno, 0),
        rlen,
        rwid,
        slope,
        top,
        rbth,
        rhk,
        roughness,
        nconn,
        ustrf,
        ndv,
    ]
    sfr_packagedata.append(rp)

# sfe data
strm_temp = 18.0  # ($C$)
K_therm_strmbed = 2.0  # ($W/m/C$)
rbthcnd = 0.0001  # ($m$)

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

icelltype = 1  # Cell conversion type (>1: unconfined)

# GWE related parameters
strt_temp = 4.0
dispersivity = 0.0  # dispersion (remember, 1D model)
ktw = 0.5918
kts = 0.2700
rhos = 1500.0
rhow = 1000.0
cps = 760.0
cpw = 4183.0
lhv = 2454.0

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-10, 1.0
ttsmult = 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen[i], nstp[i], ttsmult))

# ### Create MODFLOW 6 GWE MT3DMS Example 1 Boundary Conditions
#
# No GWF, only Heat conduction simulated


def add_gwf_model(sim, gwfname, newton=False):
    # Instantiating MODFLOW 6 groundwater flow model
    if newton:
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            newtonoptions="NEWTON",
            save_flows=True,
            model_nam_file=f"{gwfname}.nam",
        )
    else:
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            save_flows=True,
            model_nam_file=f"{gwfname}.nam",
        )

    # Instantiating MODFLOW 6 solver for flow model
    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwfname])

    # Instantiating MODFLOW 6 discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=bot,
        idomain=1,
        pname="DIS",
        filename=f"{gwfname}.dis",
    )

    # Instantiating MODFLOW 6 storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=ss,
        sy=sy,
        iconvert=iconvert,
        steady_state=steady,
        transient=transient,
        pname="STO",
        filename=f"{gwfname}.sto",
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        pname="NPF",
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=strthd,
        filename=f"{gwfname}.ic",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC",
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # sfr file
    connectiondata = [[0, -1], [1, 0]]

    sfr_perioddata = [[0, "inflow", 1.0]]

    pname = "SFR-" + gwfname[-1]
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        save_flows=True,
        nreaches=nreaches,
        packagedata=sfr_packagedata,
        connectiondata=connectiondata,
        perioddata=sfr_perioddata,
        pname=pname,
        filename=f"{gwfname}.sfr",
    )
    fname = f"{gwfname}.sfr.obs"
    sfr_obs = {
        f"{fname}.sfrobs": [
            ("inflow", "ext-inflow", 1),
            ("outflow", "ext-outflow", 2),
        ]
    }
    sfr.obs.initialize(filename=fname, print_input=True, continuous=sfr_obs)

    return sim


def add_gwe_model(sim, gwename):
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    gwe.name_file.save_flows = True

    imsgwe = flopy.mf6.ModflowIms(
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
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=bot,
        idomain=1,
        pname="DIS",
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, pname="IC", filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, pname="ADV", filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=ktw * 86400,
        kts=kts * 86400,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 source/sink mixing package for dealing with
    # auxiliary temperature specified in WEL boundary package.
    sourcerecarray = [[]]
    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sourcerecarray,
        pname="SSM",
        filename=f"{gwename}.ssm",
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        pname="OC",
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate Streamflow Energy Transport package
    sfe_packagedata = []
    for irno in range(nreaches):
        t = (irno, strm_temp, K_therm_strmbed, rbthcnd)
        sfe_packagedata.append(t)

    sfe_perioddata = []
    sfe_perioddata.append((0, "INFLOW", strm_temp))
    flwpckname = "SFR-" + gwename[-1]

    flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=sfe_packagedata,
        reachperioddata=sfe_perioddata,
        flow_package_name=flwpckname,
        pname="SFE",
        filename=f"{gwename}.sfe",
    )

    return sim


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname1 = "gwf-" + name + "nwt1"
    gwfname2 = "gwf-" + name + "non2"
    gwename1 = "gwe-" + name + "nwt1"
    gwename2 = "gwe-" + name + "non2"

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    # Build two flow models, one with NWT, one without
    sim = add_gwf_model(sim, gwfname1, newton=True)
    sim = add_gwf_model(sim, gwfname2, newton=False)

    # Add GWE models for each of the flow models above
    sim = add_gwe_model(sim, gwename1)
    sim = add_gwe_model(sim, gwename2)

    # Add the flow-transport exchanges
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname1,
        exgmnameb=gwename1,
        pname="GWFGWE1",
        filename=f"{gwename1}.gwfgwe1",
    )

    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname2,
        exgmnameb=gwename2,
        pname="GWFGWE2",
        filename=f"{gwename2}.gwfgwe2",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwename1 = "gwe-" + name + "nwt1"
    gwename2 = "gwe-" + name + "non2"

    fpth_nwt = os.path.join(test.workspace, f"{gwename1}.ucn")
    fpth_non = os.path.join(test.workspace, f"{gwename2}.ucn")

    try:
        # load temperatures
        tobj_nwt = flopy.utils.HeadFile(
            fpth_nwt, precision="double", text="TEMPERATURE"
        )
        temp_nwt = tobj_nwt.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth_nwt}"'

    try:
        # load temperatures
        tobj_non = flopy.utils.HeadFile(
            fpth_non, precision="double", text="TEMPERATURE"
        )
        temp_non = tobj_non.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth_non}"'

    # Ensure constant temperatures are initiated properly in the 1st and 3rd
    # stress periods, which are separated by period of "turning off" the
    # constant temperature boundary
    msg0 = (
        "Grid cell temperatures for the Newton and non-Newton base are "
        "different and should NOT be"
    )
    assert np.allclose(temp_nwt, temp_non), msg0

    msg1 = (
        "All layer 2 cells should be strictly colder than the cell above "
        "owing to a conductive flux of energy from the stream into the "
        "upper layer, but that is not the case"
    )
    assert np.all(np.diff(temp_nwt.squeeze(), axis=1) < 0), msg1

    msg2 = (
        "Each cell's temperatures should be monotonically increasing "
        "through time, but are not"
    )
    temp_diff = np.diff(temp_nwt.squeeze(), axis=0)

    # Because of the time discretization scheme, need to check each cell in two chunks
    # Cell ID: (0, 0, 0)
    assert isMonotonic(temp_diff[:, 0, 0][::-1][0:10]), msg2
    assert isMonotonic(temp_diff[:, 0, 0][::-1][10:-1]), msg2
    # Cell ID: (0, 1, 0)
    assert isMonotonic(temp_diff[:, 0, 1][::-1][0:10]), msg2
    assert isMonotonic(temp_diff[:, 0, 1][::-1][10:-1]), msg2
    # Cell ID: (1, 0, 0)
    assert isMonotonic(temp_diff[:, 1, 0][::-1][0:10]), msg2
    assert isMonotonic(temp_diff[:, 1, 0][::-1][10:-1]), msg2
    # Cell ID: (1, 1, 0)
    assert isMonotonic(temp_diff[:, 1, 1][::-1][0:10]), msg2
    assert isMonotonic(temp_diff[:, 1, 1][::-1][10:-1]), msg2


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
