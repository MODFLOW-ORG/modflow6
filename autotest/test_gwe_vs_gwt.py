"""
Two-Dimensional Transport in a Radial Flow Field Comparison of
MODFLOW 6 GWT with GWE

The purpose of this script is to test the new heat transport model developed
for MODFLOW 6.  To that end, this problem uses the setup of MT3DMS example problem
#5 but adapts it for heat. MODFLOW 6 is setup using GWT and the new GWE
model with input parameters entered in their native units.  The equivalent
values are calculated for "tricking" GWT into heat transport for comparison
between the two.
"""

# Imports
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace
cases = ["t_vs_c"]

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
c0 = 100.0
dt0 = 0.3
ath1 = al * trpt
botm = [top - delz]  # Model geometry
k33 = k11  # Vertical hydraulic conductivity ($m/d$)
icelltype = 0
mixelm = -1
strt = np.zeros((nlay, nrow, ncol), dtype=float)

# Active model domain
idomain = np.ones((nlay, nrow, ncol), dtype=int)
icbund = 1

# Boundary conditions
# MF2K5 pumping info:

welspd = {0: [[0, 15, 15, qwell]]}  # Well pumping info for MF2K5
spd = {0: [0, 15, 15, twell, -1]}  # Well pumping info for MT3DMS

# MF6 pumping information

#              (k,  i,  j),  flow,   conc
spd_mf6 = {0: [[(0, 15, 15), qwell, c0]]}

# MF6 constant head boundaries:

chdspd = []
# Loop through the left & right sides.
for i in np.arange(nrow):
    chdspd.append([(0, i, 0), strt[0, i, 0]])
    chdspd.append([(0, i, ncol - 1), strt[0, i, ncol - 1]])
# Loop through the top & bottom while omitting the corner cells
for j in np.arange(1, ncol - 1):
    chdspd.append([(0, 0, j), strt[0, 0, j]])
    chdspd.append([(0, nrow - 1, j), strt[0, nrow - 1, j]])

chdspd = {0: chdspd}

# Solver settings

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0
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

## Functions to build, write, and run models and plot MT3DMS Example 10 Problem results
#
# MODFLOW 6 flopy simulation object (sim) is returned if building the model


def build_models(idx, test):
    # MODFLOW 6
    ws = test.workspace
    name = cases[idx]
    gwfname = "gwf_" + name
    gwename = "gwe_" + name
    gwtname = "gwt_" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=ws,
        exe_name="mf6",
    )

    # Instantiating MODFLOW 6 time discretization
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

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
    sim.register_ims_package(imsgwf, [gwf.name])

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
        idomain=idomain,
        pname="DIS-1",
        filename=f"{gwfname}.dis",
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        pname="NPF-1",
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 storage package (steady flow conditions,
    # so no actual storage, using to print values in .lst file)
    flopy.mf6.ModflowGwfsto(gwf, ss=0, sy=0, pname="STO-1", filename=f"{gwfname}.sto")

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt, pname="IC-1", filename=f"{gwfname}.ic")

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
        filename=f"{gwfname}.chd",
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
        filename=f"{gwfname}.wel",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.bud",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # ----------------------------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ----------------------------------------------------

    # Instantiating MODFLOW 6 groundwater heat transport package
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    gwe.name_file.save_flows = True

    # create iterative model solution and register the gwe model with it
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

    # Instantiating MODFLOW 6 heat transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 heat transport initial temperatures
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, pname="IC-1", filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 heat transport advection package
    if mixelm >= 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 heat transport conduction package
    if al != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            alh=al,
            ath1=ath1,
            ktw=0.5918,
            kts=0.2700,
            filename=f"{gwename}.cnd",
        )

    # Instantiating MODFLOW 6 heat transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 heat transport source-sink mixing package
    sourcerecarray = [("WEL-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sourcerecarray,
        pname="SSM-1",
        filename=f"{gwename}.ssm",
    )

    # Instantiating MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )

    # ----------------------------------------------------
    # Instantiating MODFLOW 6 GWT model
    # ----------------------------------------------------

    # Instantiating MODFLOW 6 groundwater transport package to simulate heat
    # transport the now "old-fashion way"
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

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
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGwtic(gwt, strt=strt_temp, pname="IC-1", filename=f"{gwtname}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGwtadv(
        gwt, scheme=scheme, pname="ADV-1", filename=f"{gwtname}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    if al != 0:
        flopy.mf6.ModflowGwtdsp(
            gwt,
            alh=al,
            ath1=ath1,
            diffc=dmcoef,
            pname="DSP-1",
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
        distcoef=kd,
        filename=f"{gwtname}.mst",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    sourcerecarray = [("WEL-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, pname="SSM-1", filename=f"{gwtname}.ssm"
    )

    # Instantiating MODFLOW 6 transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE and GWT models
    name = cases[idx]
    gwename = "gwe_" + name
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    try:
        # load temperatures
        tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        temps = tobj.get_alldata()
    except:
        assert False, f'could not load temperature data from "{fpth}"'

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        # load temperatures (though stored as "concentrations")
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_alldata()
    except:
        assert False, f'could not load concentration data from "{fpth}"'

    # There are some differences between the GWE and GWT solutions.
    # For the time being, the solutions are close in that they have
    # the same shape; however, as the energy spreads outward in a
    # radial flow field simulated by a rectalinear grid, there appears
    # to be a halo that forms when differencing the GWE and GWT solutions.
    # The maximum of that difference occurs along row 18 (0-based).  The
    # entire row containing the max difference is stored here for comparison

    diff_ans = np.array(
        [
            1.23070342e-08,
            1.45018687e-07,
            1.58995667e-06,
            1.57489411e-05,
            1.37313132e-04,
            1.02408164e-03,
            6.32680308e-03,
            3.12276340e-02,
            1.18381212e-01,
            3.32011273e-01,
            6.72832216e-01,
            9.91435902e-01,
            1.11198602e00,
            1.02310179e00,
            8.56220661e-01,
            7.61168469e-01,
            8.56220661e-01,
            1.02310179e00,
            1.11198602e00,
            9.91435902e-01,
            6.72832216e-01,
            3.32011273e-01,
            1.18381212e-01,
            3.12276340e-02,
            6.32680308e-03,
            1.02408164e-03,
            1.37313132e-04,
            1.57489411e-05,
            1.58995667e-06,
            1.45018687e-07,
            1.23070339e-08,
        ]
    )

    diff = temps[0, 0] - conc[0, 0]
    assert diff.max() <= diff_ans.max(), "Max difference should be <= 1.112"
    assert np.allclose(diff[18, :], diff_ans), (
        "The difference between a GWE and equivalent GWT solution has "
        "changed. Some investigation required."
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
