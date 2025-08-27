"""
Test problem for GWE

Test the energy "flowing" through a dry cell. The test checks for
some of the flow-through energy being left behind and warming the
cell it passes through.  Based on the model appearing in the
MT3D-USGS documentation, pages 13-14.  Dry cell is in layer 1, row 1
column 4.

     +-------+-------+-------+-------+-------+-------+
  -> |   ->  |       |       | DRY   |       |       |
     |       |  ->   |       |   CELL|       |       |
     |       |       |  ->  -+-> |   |       |       |
     +-------+-------+-------+---+---+-------+-------+
     |       |       |       |   v   |       |       |
  -> |   ->  |  ->   |  ->   |  ->   |  ->   |  ->   | ->
     |       |       |       |       |       |       |
     +-------+-------+-------+-------+-------+-------+

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

cases = ["drycell0"]

# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

nrow = 1
ncol = 6
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
rhow = 1000.0
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

chd_data = {}
chd_data[0] = [
    [(0, 0, 0), left_hd],
    [(1, 0, 0), left_hd],
    [(1, 0, ncol - 1), right_hd],
]
chd_mf6 = chd_data

dispersivity = 0.0  # dispersion (remember, 1D model)

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-9, 1.0
ttsmult = 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen[i], nstp[i], ttsmult))

# ### Generate MODFLOW 6 Example test model
#


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwename1 = "gwe-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        newtonoptions="UNDER_RELAXATION",
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
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol)),
        pname="DIS-1",
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
        pname="NPF-1",
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_mf6,
        filename=f"{gwfname}.chd",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt_hd, filename=f"{gwfname}.ic")

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC-1",
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # ----------------------------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ----------------------------------------------------
    gwe1 = flopy.mf6.ModflowGwe(
        sim, modelname=gwename1, model_nam_file=f"{gwename1}.nam"
    )
    gwe1.name_file.save_flows = True
    imsgwe1 = flopy.mf6.ModflowIms(
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
        filename=f"{gwename1}.ims",
    )
    sim.register_ims_package(imsgwe1, [gwe1.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe1,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        pname="DIS-2",
        filename=f"{gwename1}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe1, strt=strt_temp, pname="IC-2", filename=f"{gwename1}.ic"
    )

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(
        gwe1, scheme=scheme, pname="ADV-2", filename=f"{gwename1}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe1,
        xt3d_off=False,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        # ktw=0.0,
        kts=0.2700 * 86400,
        # kts=0.0,
        pname="CND-2",
        filename=f"{gwename1}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe1,
        save_flows=True,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST-2",
        filename=f"{gwename1}.est",
    )

    # Instantiating MODFLOW 6 heat transport source-sink mixing package
    flopy.mf6.ModflowGwessm(
        gwe1, sources=[[]], pname="SSM-2", filename=f"{gwename1}.ssm"
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe1,
        pname="OC-2",
        budget_filerecord=f"{gwename1}.cbc",
        temperature_filerecord=f"{gwename1}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate a constant temperature in 1 of the dry cells
    ctmpspd = {
        0: [[(0, 0, 0), strt_temp], [(1, 0, 0), strt_temp]],
        1: [[(0, 0, 0), strt_temp + 10], [(1, 0, 0), strt_temp + 10]],
    }
    flopy.mf6.ModflowGwectp(
        gwe1,
        stress_period_data=ctmpspd,
        pname="CTP-2",
        filename=f"{gwename1}.ctmp",
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename1,
        pname="GWFGWE1",
        filename=f"{gwename1}.gwfgwe1",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwename = "gwe-" + name

    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    try:
        # load temperatures
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        conc1 = cobj.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth}"'

    # Check that the two perpetually dry cells:
    # 1) after warming begins (2nd stress period onward), should be greater
    #    than their initial condition
    msg0 = "There should be warming in a dry cell via conduction"
    assert np.all(conc1[1:, 0, 0, 4] > 10.0), msg0
    assert np.all(conc1[1:, 0, 0, 5] > 10.0), msg0

    msg1 = (
        "Cell at 1, 1, 5 should be warmer than the cell at 1, 1, 6 "
        "throughout the simulation by virtue of it being physically "
        "upstream"
    )
    assert np.all(conc1[1:, 0, 0, 4] > conc1[1:, 0, 0, 5]), msg1

    # 2) monotonically increase from being in contact with the warmer
    #    water passing by the cells below
    msg2 = (
        "Perpetually dry cell should be steadily warming as a result of it "
        "being in contact with warming cell below it."
    )
    assert np.all(np.diff(conc1[:, 0, 0, 5]) > 0), msg2
    assert np.all(np.diff(conc1[:, 0, 0, 4]) > 0), msg2

    # Because this is a steady flow problem, the largest incremental warming
    # increments happen in the first two stress periods.  After that, the amount
    # of warming from time step to time step decreases.
    msg3 = (
        "After the first two stress periods, the relative amount of "
        "warming in layer 1, row 1, and column 4 should slow, but isn't"
    )
    assert np.all(np.diff(np.diff(conc1[2:, 0, 0, 4])) < 0), msg3

    # The 'pass-through' cell (layer 1, row 1, column 4 - see diagram at top
    # of script) should be warming more than its two neighbors to the right.
    msg4 = "Pass through cell should be warming up at a higher rate than the dry cells."
    assert np.all(conc1[:, 0, 0, 3] > conc1[:, 0, 0, 4]), msg4

    # Pass through cell should not be as warm as the cell from which it
    # receives water, since that cell will have already robbed the water
    # passing through of some of its heat
    msg5 = "Pass through cell should not be as warm as its neighbor to the left"
    assert np.all(np.round(conc1[:, 0, 0, 3] - conc1[:, 0, 0, 2], 8) <= 0), msg5


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
