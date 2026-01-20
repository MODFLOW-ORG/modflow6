"""
Test problem for GWE

Test the obs utility by printing out the energy "flowing" between the
middle two cells.

Check energy flow including:
 1)  advection and conduction
 2)  conduction only

2 different model configurations:

Horizontal:

       +---------+---------+---------+---------+
       |         |         |         |         |
       |         |         |         |         |
       +---------+---------+---------+---------+

Vertical:

       +---------+
       |         |
       |         |
       +---------+
       |         |
       |         |
       +---------+
       |         |
       |         |
       +---------+
       |         |
       |         |
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
    "hrz-obs",  # 4-cell model, horizontally connected, check adv & cnd
    "hrz-badobs",  # Specify a 'concentration observation in a gwe model - should fail
    "vert-obs",  # 4-cell model, vertically connected, check adv & cnd
]

conn_types = (
    "horizontal",
    "horizontal",
    "vertical",
)

dis_data = {
    conn_types[0]: {
        "nrow": 1,
        "ncol": 4,
        "nlay": 1,
        "top": np.ones(4, dtype=float),
        "bot": np.zeros(4, dtype=float),
        "strthd": np.ones(4, dtype=float),
    },
    conn_types[1]: {
        "nrow": 1,
        "ncol": 4,
        "nlay": 1,
        "top": np.ones(4, dtype=float),
        "bot": np.zeros(4, dtype=float),
        "strthd": np.ones(4, dtype=float),
    },
    conn_types[2]: {
        "nrow": 1,
        "ncol": 1,
        "nlay": 4,
        "top": 4,
        "bot": np.array([[[3.0]], [[2.0]], [[1.0]], [[0.0]]], dtype=float),
        "strthd": np.ones(4, dtype=float) * 4.0,
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
perlen = [1, 1, 1, 1]  # Simulation time ($days$)
nstp = [1, 1, 1, 1]  # 10 day transient time steps
steady = {0: False}
transient = {0: True}

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
strt_temp = 4.0
dispersivity = 0.0  # dispersion (remember, 1D model)

# GWE related parameters
rhow = 1000.0
cpw = 4183.0
cps = 760.0
rhos = 1500.0
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


def build_gwf_model(sim, gwfname, idx, head1=2.0, head2=2.0):
    conn_type = conn_types[idx]

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
        pname="DIS" + gwfname[-2:],
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
        pname="NPF" + gwfname[-2:],
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=dis_data[conn_type]["strthd"],
        filename=f"{gwfname}.ic",
    )

    # Setup heads for advection case
    chd_spd = {}
    if dis_data[conn_type]["nlay"] == 1:
        chd1 = (0, 0, 0)
        chd2 = (0, 0, dis_data[conn_type]["ncol"] - 1)
        chd_spd.update({0: [[chd1, head1], [chd2, head2]]})
    elif dis_data[conn_type]["ncol"] == 1:
        chd1 = (0, 0, 0)
        chd2 = (dis_data[conn_type]["nlay"] - 1, 0, 0)
        chd_spd.update({0: [[chd1, head1], [chd2, head2]]})

    # Instantiate constant head
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chd_spd[0]),
        stress_period_data=chd_spd,
        save_flows=False,
        print_flows=True,
        pname="CHD" + gwfname[-2:],
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC" + gwfname[-2:],
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return gwf


def build_gwe_model(sim, gwename, idx):
    conn_type = conn_types[idx]

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
        nlay=dis_data[conn_type]["nlay"],
        nrow=dis_data[conn_type]["nrow"],
        ncol=dis_data[conn_type]["ncol"],
        delr=delr,
        delc=delc,
        top=dis_data[conn_type]["top"],
        botm=dis_data[conn_type]["bot"],
        idomain=1,
        pname="DIS",
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 energy transport initial temperature
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, pname="IC", filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, pname="ADV", filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 energy transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        kts=0.2700 * 86400,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport energy storage and transfer package
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

    # Instantiate source/sink mixing package
    flopy.mf6.ModflowGwtssm(
        gwe, sources=[[]], print_flows=True, filename=f"{gwename}.ssm"
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

    # Instantiate a constant temperature in the entry & exit cells
    if dis_data[conn_type]["nlay"] == 1:
        ctp1 = (0, 0, 0)
        ctp2 = (0, 0, dis_data[conn_type]["ncol"] - 1)
    elif dis_data[conn_type]["ncol"] == 1:
        ctp1 = (0, 0, 0)
        ctp2 = (dis_data[conn_type]["nlay"] - 1, 0, 0)

    ctp_spd = {0: [[ctp1, strt_temp + 2.0], [ctp2, strt_temp - 2.0]]}

    flopy.mf6.ModflowGwectp(
        gwe,
        stress_period_data=ctp_spd,
        pname="CTP",
        filename=f"{gwename}.ctmp",
    )

    # Test observations
    obs_txt = "temperature"
    if "bad" in gwename:
        obs_txt = "concentration"

    if dis_data[conn_type]["nlay"] == 1:
        obs_col1 = 1
        obs_col2 = 2
        obs_lay1 = 0
        obs_lay2 = 0
    elif dis_data[conn_type]["ncol"] == 1:
        obs_col1 = 0
        obs_col2 = 0
        obs_lay1 = 1
        obs_lay2 = 2

    obs_data0 = [
        ("temp1", obs_txt, (obs_lay2, 0, obs_col2), (obs_lay1, 0, obs_col1)),
        ("flow1", "flow-ja-face", (obs_lay2, 0, obs_col2), (obs_lay1, 0, obs_col1)),
    ]
    obs_recarray = {f"{gwename}.obs.csv": obs_data0}

    flopy.mf6.ModflowUtlobs(
        gwe,
        pname="gwe_obs",
        filename=f"{gwename}.obs",
        digits=15,
        print_input=True,
        continuous=obs_recarray,
    )

    return gwe


def add_flow_transport_exchange(sim, gwf, gwe):
    # Instantiate flow-transport exchange for each gwf-gwe model pair
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwf.name,
        exgmnameb=gwe.name,
        pname="GWFGWE",
        filename=f"{gwf.name}.gwfgwe",
    )


def build_models(idx, test):
    conn_type = conn_types[idx]

    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    gwf1 = build_gwf_model(sim, gwfname + "-1", idx, 10.0, 7.0)
    gwf2 = build_gwf_model(sim, gwfname + "-2", idx, 4.0, 4.0)

    # Instantiating MODFLOW 6 GWE model
    gwe1 = build_gwe_model(sim, gwename + "-1", idx)
    gwe2 = build_gwe_model(sim, gwename + "-2", idx)

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    add_flow_transport_exchange(sim, gwf1, gwe1)
    add_flow_transport_exchange(sim, gwf2, gwe2)

    return sim, None


def check_output(idx, test):
    print("evaluating results...")
    name = cases[idx]

    msg0 = "Simulation with advection should be greater than just conduction"
    msg1 = "Results should be monotonic, but are not"

    with open(test.workspace / "mfsim.lst", "r") as f:
        lines = f.readlines()
        error_count = 0
        for line in lines:
            if "error report" in line.lower():
                error_count += 1

    if "bad" in name:
        assert error_count > 0, "Model should've exited with an error msg"
        # xfail will catch the expected error for case 2, so
        # return out of function
        return None
    else:
        assert error_count == 0

    name = cases[idx]
    gwename1 = "gwe-" + name + "-1"
    gwename2 = "gwe-" + name + "-2"

    # load the gwt observation file
    fname1 = gwename1 + ".obs.csv"
    fname2 = gwename2 + ".obs.csv"
    fname1 = os.path.join(test.workspace, fname1)
    fname2 = os.path.join(test.workspace, fname2)
    gweobs1 = np.genfromtxt(fname1, names=True, delimiter=",", deletechars="")
    gweobs2 = np.genfromtxt(fname2, names=True, delimiter=",", deletechars="")

    # Compare results
    advAndDsp_temp = []
    advAndDsp_enerFlux = []
    justDsp_temp = []
    justDsp_enerFlux = []
    for i in np.arange(len(gweobs1)):
        advAndDsp_temp.append(gweobs1[i][1])
        advAndDsp_enerFlux.append(gweobs1[i][2])
        justDsp_temp.append(gweobs2[i][1])
        justDsp_enerFlux.append(gweobs2[i][1])

    # Ensure that the case with advection is strictly greater than
    # with conduction only
    assert np.all(np.greater(advAndDsp_temp, justDsp_temp)), msg0
    assert np.all(np.greater(advAndDsp_enerFlux, justDsp_enerFlux)), msg0

    # Ensure both the adv/cnd and conduction-only results are Monotonic
    assert isMonotonic(advAndDsp_temp), msg1
    assert isMonotonic(justDsp_temp), msg1
    assert isMonotonic(advAndDsp_enerFlux), msg1
    assert isMonotonic(justDsp_enerFlux), msg1


# - No need to change any code below
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
