"""
Test problem for GWE

Test conduction from a partially saturated group of cells into their dry
neighbors.  Referring to this test as a flowing through problem.


  Profile view of columns w/ approximate water table profile shown
     +-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
     |-------|_______|_______|       |       |       |       |       |       |       |       |       |
     |       |       |       |-------|_______|       |       |       |       |       |       |       |
     |       |       |       |       |       |-------|_______|       |       |       |       |       |
     |       |       |       |       |       |       |       |-------|       |       |       |       |
     |       |       |       |       |       |       |       |       |-------|       |       |       |
     |       |       |       |       |       |       |       |       |       |-------|       |       |
     |       |       |       |       |       |       |       |       |       |       |-------|       |
     |       |       |       |       |       |       |       |       |       |       |       |-------|
     |       |       |       |       |       |       |       |       |       |       |       |       |
     |       |       |       |       |       |       |       |       |       |       |       |       |
     |       |       |       |       |       |       |       |       |       |       |       |       |
     |       |       |       |       |       |       |       |       |       |       |       |       |
     +-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+

         ^
         |
         v

  Profile view of rows
     +-------+-------+
     |       |       |
     |       |       | <- dry cell that's checked for warming owing to increasing
     |  v    +-------+                          temperature of neighboring cell
     |-------| <- water table
     |       |
     |       |
     |       |
     |       |
     |       |
     |       |
     |       |
     |       |
     +-------+
"""  # noqa

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace


def my_ceil(a, precision=0):
    return np.round(a + 0.5 * 10 ** (-precision), precision)


def my_floor(a, precision=0):
    return np.round(a - 0.5 * 10 ** (-precision), precision)


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


scheme = "UPSTREAM"
# scheme = "TVD"

cases = ["drycl-cnduct"]

# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

nlay = 1  # Number of layers
ncol = 12  # Number of columns
nrow = 2  # Number of rows
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
top = 10.0  # Top of the model ($m$)
k11 = 4.61426  # Horizontal hydraulic conductivity ($m/d$)
ss = 1e-6  # Specific storage
sy = 0.20  # Specific Yield
nper = 11  # Number of periods
perlen = 1  # Simulation time ($days$)
nstp = 1  # One day time steps
steady = {0: True, 1: False}
transient = {0: False, 1: True}


# Set some static model parameter values
k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

bot_r1 = np.zeros(ncol).tolist()
# This is the head solution for row 1, so round up to
# ensure neighboring cells remain dry
r2 = [
    6.99850,
    6.78111,
    6.55650,
    6.32391,
    6.08242,
    5.83092,
    5.56805,
    5.29211,
    5.00092,
    4.69161,
    4.36032,
    4.00150,
]
bot_r2 = [my_ceil(val, precision=1) for val in r2]
botm = np.array([bot_r1, bot_r2])

strtdry = [my_floor(val, precision=2) for val in r2]
strt = np.array([strtdry, strtdry])  # Starting head ($m$)
icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
prsity = 0.20  # Porosity
strt_conc = np.ones((nlay, nrow, ncol), dtype=float) * 4.0
strt_temp = np.ones((nlay, nrow, ncol), dtype=float) * 4.0
dispersivity = 0.1  # dispersion (remember, 1D model)
cpw = 4183.0
cps = 760.0
rhow = 1000.0
rhos = 1500.0
lhv = 2454.0
ktw = 0.5918 * 86400
kts = 0.2700 * 86400

# Parameter equivalents:
# ----------------------
# "Distribution Coefficient"
Kd = cps / (cpw * rhow)
# "Bulk Density"
rhob = (1 - prsity) * rhos
# "Molecular Diffusion"
Kt_bulk = prsity * ktw + (1 - prsity) * kts
Dm = Kt_bulk / (prsity * rhow * cpw)

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-10, 1.0
ttsmult = 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen, nstp, ttsmult))

# ### Create MODFLOW 6 GWE MT3DMS Example 1 Boundary Conditions
#
# Constant head cells are specified on both ends of the model

# Scenario 1a (GHB to GHB)
ghbcond = 1000.0
ghb_conc = 4.0
ghb_temp = 4.0
ghb_conc_warmup = 30.0
ghb_temp_warmup = 30.0
# left boundary:  cellid, elv, cond, conc, temp
# right boundary: cellid, elv, cond, conc, temp
ghbspd = {
    # Steady state stress period
    0: [
        [(0, 0, 0), 7.0, ghbcond, ghb_conc, ghb_temp],
        [(0, 0, ncol - 1), 4.0, ghbcond, ghb_conc, ghb_temp],
    ],
    # First transient stress period
    1: [
        [(0, 0, 0), 7.0, ghbcond, ghb_conc_warmup, ghb_temp_warmup],
        [(0, 0, ncol - 1), 4.0, ghbcond, ghb_conc, ghb_temp],
    ],
}


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name
    gwename = "gwe-" + name

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
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # Instantiate left and right general-head boundaries for driving flow
    flopy.mf6.ModflowGwfghb(
        gwf,
        maxbound=len(ghbspd[0]),
        stress_period_data=ghbspd,
        auxiliary=["CONCENTRATION", "TEMPERATURE"],
        save_flows=True,
        pname="GHB-1",
        filename=f"{gwfname}.ghb",
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
    # Instantiating MODFLOW 6 GWT model
    # ----------------------------------------------------
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True
    imsgwt = flopy.mf6.ModflowIms(
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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
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
        filename=f"{gwtname}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGwtic(gwt, strt=strt_conc, pname="IC-2", filename=f"{gwtname}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGwtadv(
        gwt, scheme=scheme, pname="ADV-2", filename=f"{gwtname}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    if dispersivity != 0:
        flopy.mf6.ModflowGwtdsp(
            gwt,
            xt3d_off=True,
            alh=dispersivity,
            ath1=dispersivity,
            atv=dispersivity,
            diffc=Dm,
            pname="DSP-2",
            filename=f"{gwtname}.dsp",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=prsity,
        first_order_decay=False,
        decay=None,
        decay_sorbed=None,
        sorption="linear",
        bulk_density=rhob,
        distcoef=Kd,
        save_flows=True,
        pname="MST-2",
        filename=f"{gwtname}.mst",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    global pname
    srctype = "AUX"
    auxname = "CONCENTRATION"
    pname = ["GHB-1"]
    # Inpput to SSM is: <pname> <srctype> <auxname>
    sources = [[itm, srctype, auxname] for itm in pname]

    flopy.mf6.ModflowGwtssm(
        gwt, sources=sources, pname="SSM-2", filename=f"{gwtname}.ssm"
    )

    # Instantiate MODFLOW 6 solute transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        pname="OC-2",
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        pname="GWFGWT",
        filename=f"{gwtname}.gwfgwt",
    )

    # ---------------------------------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ---------------------------------------------------------
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
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
        botm=botm,
        idomain=1,
        pname="DIS-3",
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, pname="IC-3", filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(
        gwe, scheme=scheme, pname="ADV-3", filename=f"{gwename}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    if dispersivity != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            xt3d_off=False,
            alh=dispersivity,
            ath1=dispersivity,
            ktw=ktw,
            kts=kts,
            pname="CND-3",
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
        pname="EST-3",
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    srctype = "AUX"
    auxname = "TEMPERATURE"
    pname = ["GHB-1"]
    # Inpput to SSM is: <pname> <srctype> <auxname>
    sources = [[itm, srctype, auxname] for itm in pname]

    flopy.mf6.ModflowGwessm(
        gwe, sources=sources, pname="SSM-3", filename=f"{gwename}.ssm"
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        pname="OC-3",
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        pname="GWFGWE",
        filename=f"{gwename}.gwfgwe",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name
    gwename = "gwe-" + name

    # Pull gw heads from GWF
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        # load temperatures
        hobj = flopy.utils.HeadFile(fpth, precision="double", text="HEAD")
        hds = hobj.get_alldata()
    except:
        assert False, f'could not load head data from "{fpth}"'

    # Pull parameter-equivalent calculated temperatures from GWT
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        # load temperatures
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc1 = cobj.get_alldata()
    except:
        assert False, f'could not load concentration data from "{fpth}"'

    # Pull natively-calculated temperature output from GWE
    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    try:
        # load temperatures
        tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        temp1 = tobj.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth}"'

    # Check heads satisfy problem set up (i.e., all of row 2 is dry)
    hdsr1 = hds[0, 0, 0, :]
    assert np.all(hdsr1 < bot_r2), (
        "heads in row 1 should be below bottom elevation of row 2"
    )
    assert np.all(hds[0, 0, 1, :] < 0), "row 2 is not dry"

    # Starting temperatures after steady flow period should be 4.0 degrees
    (
        np.all(np.isclose(temp1[0], strt_temp, atol=1e-10)),
        "Steady state temperatures not as expected",
    )
    # Same with concentrations in non-dry cells
    assert np.all(np.isclose(conc1[0, :, 0], strt_conc[:, 0])), (
        "Steady state concentrations not as expected"
    )
    # Unlike GWE, GWT will not keep dry cells active and output should reflect this
    assert np.all(conc1[:, :, 1, :] < 0), (
        "Concentrations should be set to 'inactive' (-1.e+30) in row 2"
    )

    # Starting in the transient stress period, the water entering in the
    # 'through' row is warmed to 30.0 C.  First check that this is the case
    assert np.all(temp1[-1] > temp1[0]), (
        "Transient period temperature increase does not appear to have kicked in"
    )
    # Dry cell temperatures are only warmed through conduction with neighboring
    # 'through' cells and shouldn't be as warm as wet cells at the end of the
    # warming period
    assert np.all(temp1[-1, :, 0] > temp1[-1, 0, 1]), (
        "Cells with water should be warmer than dry cells"
    )

    # None of the cells should reach the temperature of the incoming water
    # during the simulation period owing to the thermal bleeding that occurs
    # to the adjacent dry cells.
    assert np.all(temp1[-1] < ghb_temp_warmup), (
        "Cells should not reach the temperature of the incoming water"
    )

    # An increasing amount of thermal bleeding should occur moving in the
    # direction of flow since the thickness of the dry cells increases in the
    # downstream direction
    assert np.all(np.diff(temp1[-1, :, 0, :-1]) < 0), (
        "Temperature change in the downstream direction should be negative"
    )
    assert isMonotonic(np.diff(temp1[-1, :, 0, :-1])), (
        "A monotonic increase in the amount of heat lost to neighboring dry "
        "cells is expected"
    )

    # Check temporal changes in temperature in the most upstream and downstream
    # dry cells. Cell bottoms in row 2 were calculated using a rounding function
    # to the nearest tenth place and as a result a monotonic tapering-off of
    # the temperature increase in the in the interior cells cannot be
    # guaranteed because of conduction with not only the wet cell in row 1, but
    # also conduction among its two dry neighbors (also with variable thickness
    colid = 0
    m_arr = np.diff(temp1[1:, 0, 1, colid])
    assert isMonotonic(m_arr), (
        "Temperatures should be monotonically tapering-off in their "
        "relative temperature increase with time in the upstream-most"
        "dry cell"
    )

    colid = 11
    m_arr = np.diff(temp1[1:, 0, 1, colid])
    assert isMonotonic(m_arr), (
        "Temperatures should be monotonically tapering-off in their "
        "relative temperature increase with time in the downstream-most"
        "dry cell"
    )

    # Run a few checks between GWE and its GWT counterpart
    # In GWT there is no solute interaction with a dry cell, so concentrations
    # in the dry cell should remain inactive (i.e., no "molecular diffusion")
    # and greater than their GWE counterpart temperatures since there is
    # no "retardation" of concentration (temperature) owing to conduction
    assert np.all(conc1[:, :, 1, :] < 0), (
        "The dry cells should never have a non-inactive concentration value"
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
