"""
Test problem for GWE designed to fail when neighboring gwe models
have different values of DENSITY_WATER specified.  For this case, it's
not clear what value should be used in the interface model.  As such,
model should exit with error msg informing the user of the issue.

Problem set up:

           "Left" Model                              "Right" Model
           ------------                              -------------

     +-------+-------+-------+                 +-------+-------+-------+
     |       |       |       |                 |       |       |       |
     |       |       |       | <- connected -> |       |       |       |
     |       |       |       |                 |       |       |       |
     +-------+-------+-------+                 +-------+-------+-------+
     |       |       |       |                 |       |       |       |
     |       |       |       | <- connected -> |       |       |       |
     |       |       |       |                 |       |       |       |
     +-------+-------+-------+                 +-------+-------+-------+

                          --->  Direction of flow  --->
"""

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

scheme = "UPSTREAM"
# scheme = "TVD"

cases = ["good_spec", "bad_spec"]
rhow = [999.0, 1001.0]

# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWF & GWE params

nrow = 1
ncol = 3
nlay = 2
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
top = 20.0  # Top of model
k11 = 100.0  # Horizontal hydraulic conductivity ($m/d$)
ss = 1e-6  # Specific storage
sy = 0.20  # Specific Yield
prsity = 0.20  # Porosity
nper = 2  # Number of periods
perlen = [1, 100]  # Simulation time ($days$)
nstp = [1, 10]  # 10 day transient time steps
steady = {0: True, 1: False}
transient = {0: False, 1: True}

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
botm = []
botm.append(np.ones((nrow, ncol), dtype=float) * 10)
botm.append(np.zeros((nrow, ncol), dtype=float))
botm = np.array(botm)

# GWE related parameters
cpw = 4183.0
lhv = 2454.0
cps = 760.0
rhos = 1500.0

# Head input
left_hd = 15.0
right_hd = 2.0
strt_hd1 = np.ones((nrow, ncol), dtype=float) * 11.0
strt_hd2 = np.ones((nrow, ncol), dtype=float) * 11.0
strt_hd1[0] = strt_hd2[0] = left_hd
strt_hd1[-1] = strt_hd2[-1] = right_hd
strt_hd = np.array([strt_hd1, strt_hd2])
strt_temp = 10.0

chd_left = {}
chd_dat_left = [
    [(0, 0, 0), left_hd],
    [(1, 0, 0), left_hd],
]
chd_left.update({0: chd_dat_left})

chd_right = {}
(chd_right.update({0: [[(1, 0, ncol - 1), right_hd]]}),)

ctp_left = {}
ctp_dat_left = [
    [(0, 0, 0), strt_temp],
    [(1, 0, 0), strt_temp],
]
ctp_left.update({0: ctp_dat_left})


dispersivity = 0.0  # dispersion (remember, 1D model)

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-8, 1.0
ttsmult = 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen[i], nstp[i], ttsmult))

# ### Generate MODFLOW 6 Example test model
#


def build_gwf_model(sim, gwfname, side="left"):
    pckg_suffix = "-1"
    if side == "right":
        pckg_suffix = "-2"

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions="UNDER_RELAXATION",
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = True

    xorigin = 0.0
    if side == "right":
        xorigin = ncol * delr

    # Instantiating MODFLOW 6 discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        xorigin=xorigin,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol)),
        pname="DIS" + pckg_suffix,
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
        pname="STO" + pckg_suffix,
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
        pname="NPF" + pckg_suffix,
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 constant head package
    if side == "left":
        chd_mf6 = chd_left
    elif side == "right":
        chd_mf6 = chd_right

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_mf6,
        pname="CHD" + pckg_suffix,
        filename=f"{gwfname}.chd",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf, strt=strt_hd, pname="IC" + pckg_suffix, filename=f"{gwfname}.ic"
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC" + pckg_suffix,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwf.set_model_relative_path(gwfname)
    return gwf


# Build GWE model
def build_gwe_model(idx, sim, gwename, side="left"):
    # Set package suffix
    pckg_suffix = "-1"
    if side == "right":
        pckg_suffix = "-2"

    # Set origin information
    xorigin = 0.0
    if side == "right":
        xorigin = ncol * delr

    # Instantiate GWE model
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    gwe.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        xorigin=xorigin,
        length_units=length_units,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        pname="DIS" + pckg_suffix,
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe, strt=strt_temp, pname="IC" + pckg_suffix, filename=f"{gwename}.ic"
    )

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(
        gwe,
        scheme=scheme,
        pname="ADV" + pckg_suffix,
        filename=f"{gwename}.adv",
    )

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=False,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        kts=0.2700 * 86400,
        pname="CND" + pckg_suffix,
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    rhow_mf6 = rhow[0]
    if idx > 0 and side == "right":
        # Set a parameter value that should trip a failure
        rhow_mf6 = rhow[idx]

    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow_mf6,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST" + pckg_suffix,
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 heat transport source-sink mixing package
    flopy.mf6.ModflowGwessm(
        gwe, sources=[[]], pname="SSM" + pckg_suffix, filename=f"{gwename}.ssm"
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        pname="OC" + pckg_suffix,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate a constant temperature in 1 of the dry cells
    if side == "left":
        flopy.mf6.ModflowGwectp(
            gwe,
            stress_period_data=ctp_left,
            pname="CTP" + pckg_suffix,
            filename=f"{gwename}.ctp",
        )

    return gwe


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwename1 = "gwe-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    # left model
    gwf1 = build_gwf_model(sim, "gwfleft", side="left")

    # right model
    gwf2 = build_gwf_model(sim, "gwfright", side="right")

    # Add the exchange data
    exgdata = [
        ((0, 0, ncol - 1), (0, 0, 0), 1, delr / 2, delr / 2, delc, 0.0, delr),
        ((1, 0, ncol - 1), (1, 0, 0), 1, delr / 2, delr / 2, delc, 0.0, delr),
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata),
        exgmnamea=gwf1.name,
        exgmnameb=gwf2.name,
        exchangedata=exgdata,
        xt3d=True,
        print_flows=True,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="exchng.gwfgwf",
        dev_interfacemodel_on=True,
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
        filename=f"{gwf1.name}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf1.name, gwf2.name])

    # left gwe model
    gwe1 = build_gwe_model(idx, sim, "gweleft", side="left")

    # right gwe model
    gwe2 = build_gwe_model(idx, sim, "gweright", side="right")

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
        filename=f"{gwe1.name}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe1.name, gwe2.name])

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea="gwfleft",
        exgmnameb="gweleft",
        pname="GWFGWEl",
        filename=f"{gwe1.name}.gwfgwe",
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea="gwfright",
        exgmnameb="gweright",
        pname="GWFGWEr",
        filename=f"{gwe2.name}.gwfgwe",
    )

    # Create GWE GWE exchange
    flopy.mf6.ModflowGwegwe(
        sim,
        exgtype="GWE6-GWE6",
        gwfmodelname1=gwf1.name,
        gwfmodelname2=gwf2.name,
        adv_scheme=scheme,
        nexg=len(exgdata),
        exgmnamea=gwe1.name,
        exgmnameb=gwe2.name,
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="exchng.gwegwe",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # Look for successful termination in simulation lst file
    srch_str = "Normal termination of simulation"
    fpth = os.path.join(test.workspace, "mfsim.lst")
    proceed = False
    with open(fpth) as f:
        lines = f.readlines()
        if any(srch_str in line for line in lines):
            proceed = True

    if proceed:
        # read transport results from GWE model
        gwename1 = "gweleft" + ".ucn"
        gwename2 = "gweright" + ".ucn"

        fpth1 = os.path.join(test.workspace, gwename1)
        fpth2 = os.path.join(test.workspace, gwename2)

        # load temperatures
        tobj1 = flopy.utils.HeadFile(fpth1, precision="double", text="TEMPERATURE")
        tobj2 = flopy.utils.HeadFile(fpth2, precision="double", text="TEMPERATURE")
        temps1 = tobj1.get_alldata()
        temps2 = tobj2.get_alldata()
        temps_all = np.concatenate((temps1, temps2), axis=3)
        np.all(temps_all == strt_temp)
        msg0 = "Model finished, but output not uniformly 10.0"
        assert np.allclose(temps_all, 10.0), msg0
    else:
        msg1 = "Model expected to fail when idx > 0"
        print(msg1)


@pytest.mark.developmode
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets):
    xfail = ["bad" in cases[ct] for ct in np.arange(len(cases))]
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=xfail[idx],
    )
    test.run()
