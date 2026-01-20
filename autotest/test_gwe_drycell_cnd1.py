"""
Test problem for GWE

Test the energy "flowing" between two dry cells via conduction
only using a temperature gradient

   ~: Represents conduction

 A) 1st model configuration

       +---------+---------+
       |         |~        |
       |         |~        |
       +---------+---------+

 B) 2nd model configuration

       +---------+
       |         |
       |  ~   ~  |
       +---------+
       |         |
       |         |
       +---------+

 C) 3rd model configuration

                 +---------+
       +---------+         |
       |         |~        |
       |         +---------+
       +---------+
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

cases = [
    "drycell2-a",  # 2-cell model, horizontally connected with tops and bots aligned
    "drycell2-b",  # 2-cell model, vertically connected
    "drycell2-c",  # 2-cell model, horizontally connected with staggered alignment
    #                              (reduced shared cell face area)
]

conn_types = (
    "horizontal",
    "vertical",
    "staggered",
)

dis_data = {
    conn_types[0]: {
        "nrow": 1,
        "ncol": 2,
        "nlay": 1,
        "top": np.ones(2, dtype=float),
        "bot": np.zeros(2, dtype=float),
        "strthd": np.zeros(2, dtype=float),
    },
    conn_types[1]: {
        "nrow": 1,
        "ncol": 1,
        "nlay": 2,
        "top": 2,
        "bot": np.array([[[1.0]], [[0.0]]], dtype=float),
        "strthd": np.zeros(2, dtype=float),
    },
    conn_types[2]: {
        "nrow": 2,
        "ncol": 1,
        "nlay": 1,
        "top": np.array([[[1.0], [1.5]]], dtype=float),
        "bot": np.array([[[0.0], [0.5]]], dtype=float),
        "strthd": -1 * np.ones(2, dtype=float),
    },
}

# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

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

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
strt_temp1 = 4.0
strt_temp2 = 34.0
dispersivity = 0.0  # dispersion (remember, 1D model)

# GWE related parameters
rhow = 1000.0
cpw = 4183.0
lhv = 2454.0
cps = 760.0
rhos = 1500.0

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


def build_models(idx, test):
    conn_type = conn_types[idx]

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
        linear_acceleration="CG",
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
        nlay=dis_data[conn_type]["nlay"],
        nrow=dis_data[conn_type]["nrow"],
        ncol=dis_data[conn_type]["ncol"],
        delr=delr,
        delc=delc,
        top=dis_data[conn_type]["top"],
        botm=dis_data[conn_type]["bot"],
        idomain=np.ones(
            (
                dis_data[conn_type]["nlay"],
                dis_data[conn_type]["nrow"],
                dis_data[conn_type]["ncol"],
            ),
            dtype=int,
        ),
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

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=dis_data[conn_type]["strthd"],
        filename=f"{gwfname}.ic",
    )

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
        filename=f"{gwename1}.ims",
    )
    sim.register_ims_package(imsgwe1, [gwe1.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe1,
        nogrb=True,
        nlay=dis_data[conn_type]["nlay"],
        nrow=dis_data[conn_type]["nrow"],
        ncol=dis_data[conn_type]["ncol"],
        delr=delr,
        delc=delc,
        top=dis_data[conn_type]["top"],
        botm=dis_data[conn_type]["bot"],
        idomain=1,
        pname="DIS-2",
        filename=f"{gwename1}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe1, strt=strt_temp1, pname="IC-2", filename=f"{gwename1}.ic"
    )

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(
        gwe1, scheme=scheme, pname="ADV-2", filename=f"{gwename1}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe1,
        xt3d_off=True,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        kts=0.2700 * 86400,
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
        0: [
            [(0, 0, 0), strt_temp1],
            [
                (
                    dis_data[conn_type]["nlay"] - 1,
                    dis_data[conn_type]["nrow"] - 1,
                    dis_data[conn_type]["ncol"] - 1,
                ),
                strt_temp1 + 10,
            ],
        ],
        1: [],
        2: [
            [(0, 0, 0), strt_temp2],
            [
                (
                    dis_data[conn_type]["nlay"] - 1,
                    dis_data[conn_type]["nrow"] - 1,
                    dis_data[conn_type]["ncol"] - 1,
                ),
                strt_temp2 - 10,
            ],
        ],
        3: [],
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

    # All indices are 0 based
    # initialize
    idxl = 0
    idxr = 0
    idxc = 0
    # override depending on scenario
    if idx == 0:
        idxc = 1

    if idx == 1:
        idxl = 1

    if idx == 2:
        idxr = 1

    fpth = os.path.join(test.workspace, f"{gwename}.ucn")

    try:
        # load temperatures
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        conc1 = cobj.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth}"'

    # Ensure constant temperatures are initiated properly in the 1st and 3rd
    # stress periods, which are separated by period of "turning off" the
    # constant temperature boundary
    msg0 = "Grid cell temperatures do not reflect user-specified difference"
    assert conc1[0, 0, 0, 0] + 10.0 == conc1[0, idxl, idxr, idxc], msg0
    assert conc1[11, 0, 0, 0] - 10.0 == conc1[11, idxl, idxr, idxc], msg0

    # After running transient stress period, temperatures in grid cells
    # should equilibrate through the process of conduction only (there
    # is no gwf flow)
    msg1 = "Grid cell temperatures should have equilabrated via conduction"
    assert np.isclose(conc1[10, 0, 0, 0], conc1[10, idxl, idxr, idxc]), msg1
    assert np.isclose(conc1[21, 0, 0, 0], conc1[21, idxl, idxr, idxc]), msg1

    # Ensure that as the cells equilibrate, they do so in a monotonic manner
    msg2 = "There should be a monotonic increase as the 2 cells equilibrate"
    msg3 = "There should be a monotonic decrease as the 2 cells equilibrate"
    assert isMonotonic(np.diff(conc1[1:11, 0, 0, 0])), msg2
    assert isMonotonic(np.diff(conc1[1:11, idxl, idxr, idxc])), msg3
    assert isMonotonic(np.diff(conc1[12:, 0, 0, 0])), msg3
    assert isMonotonic(np.diff(conc1[12:, idxl, idxr, idxc])), msg2

    # Ensure that the equilibrated temperature is half the starting difference
    # between the cells
    msg4 = (
        "The final equilibrated cell temperature does not split the "
        "difference of the starting temperature"
    )
    initTdiff1 = abs(conc1[0, 0, 0, 0] - conc1[0, idxl, idxr, idxc])
    initTdiff2 = abs(conc1[11, 0, 0, 0] - conc1[11, idxl, idxr, idxc])
    assert np.isclose(
        conc1[10, 0, 0, 0],
        conc1[0, 0, 0, 0] + initTdiff1 / 2,
    ), msg4
    assert np.isclose(
        conc1[21, 0, 0, 0],
        conc1[11, 0, 0, 0] - initTdiff2 / 2,
    ), msg4


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
