"""
Two-Dimensional Heat Transport (GWE) in a Radial Flow Field with Comparison
to a MODFLOW 6 GWT model

The purpose of this script is to test the new heat transport model developed
for MODFLOW 6.  To that end, this problem uses the setup of the fifth MT3DMS
test problem but adapts it for heat.  Moreover, this example also attempts to
adapt the fifth MT3DMS problem to a GWE-GWE example by splitting the model
domain diagonally in half.  As such, MODFLOW 6 is setup using the new GWE
model with input parameters entered in their native units.  The equivalent
values are calculated for "tricking" MT3DMS into heat transport.

Cheap depiction of model in plan view follows (not to scale):

       +---------------------------------------------------------------------------+
       |                                                                           |
       +---+                                                                       |
           |                                                                       |
           +---+                                                                   |
               |                                                                   |
   +---+       +---+                                                               |
   |   |           |                                                               |
   |   +---+       +---+                                                           |
   |       |           |                                                           |
   |       +---+       +---+              "GWx-ur"                                 |
   |           |           |                                                       |
   |           +---+       +---+                                                   |
   |               |           |                                                   |
   |               +---+       +---+                                               |
   |                   |           |                                               |
   |                   +---+       +---+                                           |
   |                       |           |                                           |
   |                       +---+       +---+                                       |
   |                           |           | * <-- Inj. well with 100 deg C        |
   |                           +---+       +---+      water entering simulation>   |
   |                               |           |      Inj. well creates a radial   |
   |                               +---+       +---+  (outward) flow field.        |
   |                                   |           |                               |
   |                                   +---+       +---+                           |
   |                                       |           |                           |
   |                                       +---+       +---+                       |
   |                                           |           |                       |
   |                                           +---+       +---+                   |
   |                                               |           |                   |
   |                                               +---+       +---+               |
   |                                                   |           |               |
   |                                                   +---+       +---+           |
   |                                  "GWx-ll"             |           |           |
   |                                                       +---+       +---+       |
   |                                                           |           |       |
   |                                                           +---+       +---+   |
   |                                                               |           |   |
   |                                                               +---+       +---+
   |                                                                   |
   |                                                                   +---+
   |                                                                       |
   |                                                                       +---+
   |                                                                           |
   +---------------------------------------------------------------------------+

"""

import os
import sys

sys.path.append(os.path.join("..", "common"))

# Imports

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pytest
from framework import TestFramework

# from figspecs import USGSFigure


# Set figure properties specific to this problem
plotModel = plotSave = False
figure_size = (6, 4.5)

cases = ["gwegwe-gwtgwt"]

# Model units

length_units = "meters"
time_units = "days"

# Table

nlay = 1  # Number of layers
nrow = 31  # Number of rows
ncol = 31  # Number of columns
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
delz = 1.0  # Layer thickness ($m$)
top = 0.0  # Top of the model ($m$)
prsity = 0.3  # Porosity
perlen = 27  # Simulation time ($days$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)
qwell = 100.0  # Volumetric injection rate ($m^3/d$)
twell = 100.0  # Concentration of injected water ($mg/L$)
al = 10.0  # Longitudinal dispersivity ($m$)
trpt = 1.0  # Ratio of transverse to longitudinal dispersitivity

dmcoef = 3.2519e-7
rhob = 1110.0
sp2 = 0.0  # read, but not used in this problem
kd = 1.8168e-4

cpw = 4183.0
cps = 760.0
rhow = 1000.0
rhos = 1500.0
lhv = 2454.0

# Additional model input

perlen = [27]
nper = len(perlen)
nstp = [27]
tsmult = [1.0]
strt_temp = 0.0
dt0 = 0.3
ath1 = al * trpt
botm = [top - delz]  # Model geometry
k33 = k11  # Vertical hydraulic conductivity ($m/d$)
icelltype = 0
mixelm = 0
strt = np.zeros((nlay, nrow, ncol), dtype=float)

# Active model domain
ibound_mf2k5 = np.ones((nlay, nrow, ncol), dtype=int) * -1
ibound_mf2k5[:, 1 : nrow - 1, 1 : ncol - 1] = 1
icbund = 1
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain_ll = idomain.copy()
for rw in np.arange(idomain_ll.shape[1]):  # For each row
    idomain_ll[0, rw, rw:] = 0

idomain_ur = 1 - idomain_ll

# Work up the exchanges
# Exchange data for GWF-GWF and GWE-GWE
exgdata = []

# Stair-step down the shared interface
for i in np.arange(nrow):
    for j in np.arange(ncol):
        # Checking between rows (i.e., the vertical (y-axis) direction)
        if i < (nrow - 1):
            # Check to see if two touching cells in adjacent models are both active
            # Check
            if idomain_ur[0, i, j] > 0 and idomain_ll[0, i + 1, j] > 0:
                exgdata.append(((0, i, j), (0, i + 1, j), 1, 5, 5, 10, 270.0, 10.0))
        if j < (ncol - 1):
            if idomain_ur[0, i, j + 1] > 0 and idomain_ll[0, i, j] > 0:
                exgdata.append(((0, i, j + 1), (0, i, j), 1, 5, 5, 10, 180.0, 10.0))


# Boundary conditions
# MF2K5 pumping info:

welspd = {0: [[0, 15, 15, qwell]]}  # Well pumping info for MF2K5
spd = {0: [0, 15, 15, twell, -1]}  # Well pupming info for MT3DMS

# MF6 pumping information

#              (k,  i,  j),  flow,  temperature
spd_mf6 = {0: [[(0, 15, 15), qwell, twell]]}

# MF6 constant head boundaries:

chdspd_ur = []
chdspd_ll = []
# Loop through the left & right sides.
for i in np.arange(nrow):
    if i == 0:
        # A special case, this one cell on the left boundary happens to fall in the
        # in the "ur" model
        chdspd_ur.append([(0, i, 0), strt[0, i, 0]])
    else:
        chdspd_ll.append([(0, i, 0), strt[0, i, 0], 0.0])

    chdspd_ur.append([(0, i, ncol - 1), strt[0, i, ncol - 1]])

# Loop through the top & bottom while omitting the corner cells
for j in np.arange(1, ncol - 1):
    chdspd_ur.append([(0, 0, j), strt[0, 0, j]])
    chdspd_ll.append([(0, nrow - 1, j), strt[0, nrow - 1, j], 0.0])

chdspd_ur = {0: chdspd_ur}
chdspd_ll = {0: chdspd_ll}

# Solver settings

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-8, 1e-8, 1.0
percel = 1.0  # HMOC parameters
itrack = 3
wd = 0.5
dceps = 1.0e-5
nplane = 1
npl = 0
nph = 16
npmin = 2
npmax = 32
dchmoc = 1.0e-3
nlsink = nplane
npsink = nph

# Static temporal data used by TDIS file

tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))

# Model names
gwfname_up = "gwf-ur"
gwfname_lo = "gwf-ll"
gwename_up = "gwe-ur"
gwename_lo = "gwe-ll"
gwtname_up = "gwt-ur"
gwtname_lo = "gwt-ll"

# Functions to build, write, and run models
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_models(idx, test):
    # -----------
    #  MODFLOW 6
    # -----------

    ws = test.workspace
    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    # add both solutions to the simulation
    add_flow(sim)
    add_energy(sim)
    add_transport(sim)

    # Instantiating MODFLOW 6 flow-energy transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname_up,
        exgmnameb=gwename_up,
        filename="upper.gwfgwe",
    )
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname_lo,
        exgmnameb=gwename_lo,
        filename="lower.gwfgwe",
    )

    # Next, instantiate MODFLOW 6 flow-solute transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname_up,
        exgmnameb=gwtname_up,
        filename="upper.gwfgwt",
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname_lo,
        exgmnameb=gwtname_lo,
        filename="lower.gwfgwt",
    )

    return sim, None


# Instantiate the upper and lower flow models
def add_flow(sim):
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
        filename="gwfsolver.ims",
    )

    gwf_upper = add_upper_gwfmodel(sim)
    gwf_lower = add_lower_gwfmodel(sim)

    sim.register_ims_package(imsgwf, [gwf_upper.name, gwf_lower.name])

    # Add the exchange data
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata),
        exgmnamea=gwf_upper.name,
        exgmnameb=gwf_lower.name,
        exchangedata=exgdata,
        xt3d=False,
        print_flows=True,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="exchng.gwfgwf",
    )


# Create the upper GWF model
def add_upper_gwfmodel(sim):
    mname = gwfname_up

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=mname,
        save_flows=True,
        model_nam_file=f"{mname}.nam",
    )

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
        idomain=idomain_ur,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        filename=f"{mname}.npf",
    )

    # Instantiating MODFLOW 6 storage package
    # (steady flow conditions, so no actual storage, using to print values in .lst file)
    flopy.mf6.ModflowGwfsto(gwf, ss=0, sy=0, filename=f"{mname}.sto")

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{mname}.ic")

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd_ur),
        stress_period_data=chdspd_ur,
        save_flows=False,
        pname="CHD-1",
        filename=f"{mname}.chd",
    )

    # Instantiate the wel package
    flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=spd_mf6,
        save_flows=False,
        auxiliary="TEMPERATURE",
        pname="WEL-1",
        filename=f"{mname}.wel",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname}.hds",
        budget_filerecord=f"{mname}.bud",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


# Create the lower GWF model
def add_lower_gwfmodel(sim):
    mname = gwfname_lo

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=mname,
        save_flows=True,
        model_nam_file=f"{mname}.nam",
    )

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
        idomain=idomain_ll,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{mname}.ic")

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        filename=f"{mname}.npf",
    )

    # Instantiating MODFLOW 6 storage package
    # (steady flow conditions, so no actual storage, using to print values in .lst file)
    flopy.mf6.ModflowGwfsto(gwf, ss=0, sy=0, filename=f"{mname}.sto")

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd_ll),
        stress_period_data=chdspd_ll,
        auxiliary="TEMPERATURE",
        save_flows=False,
        pname="CHD-1",
        filename=f"{mname}.chd",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname}.hds",
        budget_filerecord=f"{mname}.bud",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


# Create the upper and lower GWE models
def add_energy(sim):
    # create iterative model solution and register the gwt model with it
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
        filename="gwesolver.ims",
    )

    # Set the advection scheme, it is needed by both gwe model instantiation
    # and gwegwe exchange
    if mixelm >= 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()

    # Add transport models
    gwe_upper = add_upper_gwemodel(sim, scheme)
    gwe_lower = add_lower_gwemodel(sim, scheme)

    sim.register_ims_package(imsgwe, [gwe_upper.name, gwe_lower.name])

    # Create energy transport to energy transport coupling
    assert exgdata is not None
    flopy.mf6.ModflowGwegwe(
        sim,
        exgtype="GWE6-GWE6",
        gwfmodelname1=gwfname_up,
        gwfmodelname2=gwfname_lo,
        adv_scheme=scheme,
        nexg=len(exgdata),
        exgmnamea=gwe_upper.name,
        exgmnameb=gwe_lower.name,
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="exchng.gwegwe",
    )

    return sim


# Create the outer GWT model
def add_upper_gwemodel(sim, scheme):
    # Instantiating MODFLOW 6 groundwater heat transport package
    mname = gwename_up

    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=mname,
        model_nam_file=f"{mname}.nam",
    )
    gwe.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain_ur,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 heat transport initial temperature
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, filename=f"{mname}.ic")

    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, filename=f"{mname}.adv")

    # Instantiating MODFLOW 6 heat transport dispersion package
    if al != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            alh=al,
            ath1=ath1,
            ktw=0.5918,
            kts=0.2700,
            filename=f"{mname}.cnd",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST-UP",
        filename=f"{mname}.est",
    )

    # Instantiating MODFLOW 6 heat transport source-sink mixing package
    sourcerecarray = [("WEL-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{mname}.ssm")

    # Instantiating MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{mname}.cbc",
        temperature_filerecord=f"{mname}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    return gwe


# Create the lower GWE model
def add_lower_gwemodel(sim, scheme):
    # Instantiating MODFLOW 6 groundwater heat transport package
    mname = gwename_lo

    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=mname,
        model_nam_file=f"{mname}.nam",
    )
    gwe.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain_ll,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 heat transport initial temperature
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, filename=f"{mname}.ic")

    # Instantiating MODFLOW 6 heat transport advection package (lower model)
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, filename=f"{mname}.adv")

    # Instantiating MODFLOW 6 heat transport dispersion package
    if al != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            alh=al,
            ath1=ath1,
            ktw=0.5918,
            kts=0.2700,
            filename=f"{mname}.cnd",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST-LO",
        filename=f"{mname}.est",
    )

    # Instantiating MODFLOW 6 heat transport source-sink mixing package
    sourcerecarray = [("CHD-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{mname}.ssm")

    # Instantiating MODFLOW 6 heat transport output control package
    # flopy.mf6.ModflowGweoc(
    #    gwe,
    #    budget_filerecord="{}.cbc".format(mname),
    #    temperature_filerecord="{}.ucn".format(mname),
    #    temperatureprintrecord=[
    #        ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
    #    ],
    #    saverecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    #    printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    # )
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{mname}.cbc",
        temperature_filerecord=f"{mname}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    return gwe


# Add an equivalent solute transport model for comparing to
def add_transport(sim):
    # create iterative model solution and register the gwt model with it
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
        filename="gwtsolver.ims",
    )

    # Set the advection scheme, it is needed by both gwe model instantiation
    # and gwegwe exchange
    if mixelm >= 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()

    # Add transport models
    gwt_upper = add_upper_gwtmodel(sim, scheme)
    gwt_lower = add_lower_gwtmodel(sim, scheme)

    sim.register_ims_package(imsgwt, [gwt_upper.name, gwt_lower.name])

    # Create transport-transport coupling
    assert exgdata is not None
    flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1=gwfname_up,
        gwfmodelname2=gwfname_lo,
        adv_scheme=scheme,
        nexg=len(exgdata),
        exgmnamea=gwt_upper.name,
        exgmnameb=gwt_lower.name,
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="exchng.gwtgwt",
    )

    return sim


# Create the upper GWT model
def add_upper_gwtmodel(sim, scheme):
    # Instantiating MODFLOW 6 groundwater solute transport
    # model for comparison with GWE
    mname = gwtname_up

    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=mname,
        model_nam_file=f"{mname}.nam",
    )
    gwt.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain_ur,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 heat transport initial temperature
    flopy.mf6.ModflowGwtic(gwt, strt=strt_temp, filename=f"{mname}.ic")

    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{mname}.adv")

    # Instantiating MODFLOW 6 heat transport dispersion package
    if al != 0:
        flopy.mf6.ModflowGwtdsp(
            gwt,
            alh=al,
            ath1=ath1,
            diffc=dmcoef,
            filename=f"{mname}.dsp",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGwtmst(
        gwt,
        sorption="LINEAR",
        porosity=prsity,
        bulk_density=1110.0,
        distcoef=kd,
        filename=f"{mname}.mst",
    )

    # Instantiating MODFLOW 6 source-sink mixing package transport
    sourcerecarray = [("WEL-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{mname}.ssm")

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{mname}.cbc",
        concentration_filerecord=f"{mname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return gwt


# Create the lower GWT model
def add_lower_gwtmodel(sim, scheme):
    # Instantiating MODFLOW 6 groundwater solute transport package
    mname = gwtname_lo

    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=mname,
        model_nam_file=f"{mname}.nam",
    )
    gwt.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain_ll,
        filename=f"{mname}.dis",
    )

    # Instantiating MODFLOW 6 solute transport initial concentration (temperature)
    flopy.mf6.ModflowGwtic(gwt, strt=strt_temp, filename=f"{mname}.ic")

    # Instantiating MODFLOW 6 solute transport advection package (lower model)
    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{mname}.adv")

    # Instantiating MODFLOW 6 heat transport dispersion package
    if al != 0:
        flopy.mf6.ModflowGwtdsp(
            gwt,
            alh=al,
            ath1=ath1,
            diffc=dmcoef,
            filename=f"{mname}.dsp",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGwtmst(
        gwt,
        sorption="LINEAR",
        porosity=prsity,
        bulk_density=rhob,
        distcoef=kd,
        filename=f"{mname}.mst",
    )

    # Instantiating MODFLOW 6 solute transport source-sink mixing package
    sourcerecarray = [("CHD-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{mname}.ssm")

    # Instantiating MODFLOW 6 solute transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{mname}.cbc",
        concentration_filerecord=f"{mname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return gwt


def fict_model(out_pth):
    # Instantiate the MODFLOW model
    mf = flopy.modflow.Modflow(
        modelname="dummy_mod", model_ws=out_pth, exe_name="mfnwt"
    )

    # Instantiate discretization package
    # units: itmuni=4 (days), lenuni=2 (m)
    flopy.modflow.ModflowDis(
        mf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        nper=nper,
        nstp=nstp,
        perlen=perlen,
        itmuni=4,
        lenuni=2,
    )

    return mf


# Function to plot the model results
def check_output(idx, test):
    print("evaluating results...")

    out_pth = test.workspace
    # read transport results from GWE and GWT models
    fpth = os.path.join(out_pth, gwename_up + ".ucn")
    turobj = flopy.utils.HeadFile(fpth, text="temperature")
    temp_ur = turobj.get_alldata()

    fpth = os.path.join(out_pth, gwename_lo + ".ucn")
    tllobj = flopy.utils.HeadFile(fpth, text="temperature")
    temp_ll = tllobj.get_alldata()

    fpth = os.path.join(out_pth, gwtname_up + ".ucn")
    curobj = flopy.utils.HeadFile(fpth, text="concentration")
    conc_ur = curobj.get_alldata()

    fpth = os.path.join(out_pth, gwtname_lo + ".ucn")
    cllobj = flopy.utils.HeadFile(fpth, text="concentration")
    conc_ll = cllobj.get_alldata()

    # Stitch together the temperatures from the upper-right and lower-left models
    tm = k = 0
    stitched_temps = np.zeros((nrow, ncol))
    for i in np.arange(nrow):
        for j in np.arange(ncol):
            if idomain_ur[k, i, j] > 0:
                stitched_temps[i, j] = temp_ur[tm, k, i, j]
            elif idomain_ll[k, i, j] > 0:
                stitched_temps[i, j] = temp_ll[tm, k, i, j]

    # Stitch together the "concentrations" (which represent temperatures)
    # from the upper-right and lower-left models
    tm = k = 0
    stitched_conc = np.zeros((nrow, ncol))
    for i in np.arange(nrow):
        for j in np.arange(ncol):
            if idomain_ur[k, i, j] > 0:
                stitched_conc[i, j] = conc_ur[tm, k, i, j]
            elif idomain_ll[k, i, j] > 0:
                stitched_conc[i, j] = conc_ll[tm, k, i, j]

    diff_ans = np.array(
        [
            1.72851147e-06,
            1.02569324e-05,
            6.03854917e-05,
            3.24474242e-04,
            1.55996452e-03,
            6.58075914e-03,
            2.38852818e-02,
            7.31492711e-02,
            1.85545249e-01,
            3.83629779e-01,
            6.39985972e-01,
            8.62111344e-01,
            9.55210100e-01,
            9.08020311e-01,
            7.95058119e-01,
            7.17466620e-01,
            7.95058119e-01,
            9.08020311e-01,
            9.55210110e-01,
            8.62111344e-01,
            6.39985972e-01,
            3.83629779e-01,
            1.85545249e-01,
            7.31492710e-02,
            2.38852818e-02,
            6.58075913e-03,
            1.55996452e-03,
            3.24474242e-04,
            6.03854915e-05,
            1.02569324e-05,
            1.72851146e-06,
        ]
    )

    diff = stitched_temps - stitched_conc
    assert diff.max() <= diff_ans.max(), "Max difference should be <= 0.95"
    assert np.allclose(diff[18, :], diff_ans), (
        "The difference between a GWE and equivalent GWT solution has "
        "changed. Some investigation required."
    )

    # Ensure that all 6 models are present within the simulation
    mf6_sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    models = []
    for mname in mf6_sim.model_names:
        models.append(mf6_sim.get_model(mname))
    assert len(models) == 6, (
        "Unexpected number of models encountered while loading simulation"
    )

    if plotModel:
        # Create figure for scenario
        sim_name = cases[idx]
        mod = mf6_sim.get_model(mname)

        plt.rcParams["lines.dashed_pattern"] = [5.0, 5.0]

        fig = plt.figure(figsize=figure_size, dpi=300, tight_layout=True)
        ax2 = fig.add_subplot(1, 1, 1, aspect="equal")

        # temp1 = stitched_conc[:, :]
        # temp2 = stitched_temps[:, :]

        levels = [0.2, 5, 20, 50, 90]

        mod4plot = fict_model(out_pth)
        mm = flopy.plot.PlotMapView(model=mod4plot)
        mm.plot_grid(color=".5", alpha=0.2)
        mm.plot_ibound(ibound=idomain_ur, color_noflow="grey", alpha=0.1)
        cs1 = mm.contour_array(stitched_conc, levels=levels, colors="r")
        plt.clabel(cs1, inline=1, fontsize=10)
        cs2 = mm.contour_array(
            stitched_temps, levels=levels, colors="k", linestyles=":"
        )
        plt.clabel(cs2, inline=1, fontsize=10)
        labels = ["GWT", "GWE"]
        lines = [cs1.collections[0], cs2.collections[0]]

        plt.xlabel("Distance Along X-Axis, in meters")
        plt.ylabel("Distance Along Y-Axis, in meters")
        title = "Comparison of GWT and GWE isoconcentration lines"
        ax2.legend(lines, labels, loc="upper right")

        # Add labels to the plot
        # style = dict(size=10, color='black')
        # ax2.text(235, 140, "Location of x-section \nshown in Fig. x", **style)

        # save figure
        if plotSave:
            fpth = os.path.join(out_pth, f"{sim_name + '-planView.png'}")
            fig.savefig(fpth)


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
