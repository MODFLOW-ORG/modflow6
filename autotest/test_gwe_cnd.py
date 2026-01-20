"""
Test problem for GWE

One-Dimensional Transport in a Uniform Flow Field.
The purpose of this script is to test the new heat transport model developed
for MODFLOW 6.  To that end, this problem uses the setup of the first MT3DMS
test problem but adapts it for heat. MODFLOW 6 is setup using the new GWE
model with input parameters entered in their native units.

It may be possible to find a 1D heat transport analytical solution in the
future.
"""

# Imports

import os

import numpy as np
import pytest

try:
    pass
except:
    msg = "Error. Pymake package is not available.\n"
    msg += "Try installing using the following command:\n"
    msg += " pip install https://github.com/modflowpy/pymake/zipball/master"
    raise Exception(msg)

try:
    import flopy
except:
    msg = "Error. FloPy package is not available.\n"
    msg += "Try installing using the following command:\n"
    msg += " pip install flopy"
    raise Exception(msg)


from framework import TestFramework

# Base simulation and model name and workspace

viscosity_on = [False]
cases = ["cnd01"]

# Model units

length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

nper = 1  # Number of periods
nlay = 1  # Number of layers
ncol = 101  # Number of columns
nrow = 1  # Number of rows
delr = 10.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
top = 0.0  # Top of the model ($m$)
botm = -1.0  # Layer bottom elevations ($m$)
prsity = 0.25  # Porosity
perlen = 2000  # Simulation time ($days$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
laytyp = 1
nstp = 100.0
dt0 = perlen / nstp
Lx = (ncol - 1) * delr
v = 0.24
q = v * prsity
h1 = q * Lx
strt = np.zeros((nlay, nrow, ncol), dtype=float)
strt[0, 0, 0] = h1  # Starting head ($m$)
l = 1000.0  # Needed for plots
icelltype = 1  # Cell conversion type
ibound = np.ones((nlay, nrow, ncol), dtype=int)
ibound[0, 0, 0] = -1
ibound[0, 0, -1] = -1

# Set some static transport related model parameter values

mixelm = 0  # FD
rhob = 1110.0
sp2 = 0.0  # read, but not used in this problem
kd = 1.8168e-4
strt_temp = np.zeros((nlay, nrow, ncol), dtype=float)
dispersivity = 1.0
dmcoef = 3.2519e-7  # Molecular diffusion coefficient

# Set some static heat transport related model parameter values
cpw = 4183.0
rhow = 1000.0
lhv = 2454.0
cps = 760.0
rhos = 1500.0

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0
ttsmult = 1.0
dceps = 1.0e-5  # HMOC parameters in case they are invoked
nplane = 1  # HMOC
npl = 0  # HMOC
nph = 4  # HMOC
npmin = 0  # HMOC
npmax = 8  # HMOC
nlsink = nplane  # HMOC
npsink = nph  # HMOC

# Static temporal data used by TDIS file

tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))

# ### Create MODFLOW 6 GWE MT3DMS Example 1 Boundary Conditions
#
# Constant head cells are specified on both ends of the model

chdspd = [[(0, 0, 0), h1], [(0, 0, ncol - 1), 0.0]]
c0 = 40.0
ctpspd = [[(0, 0, 0), c0]]


def build_models(idx, test):
    # Base MF6 GWE model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    sim_ws = os.path.join(ws, name)
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
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # Instantiating VSC
    if viscosity_on[idx]:
        # Instantiate viscosity (VSC) package
        vsc_filerecord = f"{gwfname}.vsc.bin"
        vsc_pd = [(0, 0.0, 20.0, gwename, "temperature")]
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=8.904e-4,
            viscosity_filerecord=vsc_filerecord,
            thermal_formulation="nonlinear",
            thermal_a2=10.0,
            thermal_a3=248.37,
            thermal_a4=133.16,
            nviscspecies=len(vsc_pd),
            packagedata=vsc_pd,
            pname="vsc",
            filename=f"{gwfname}.vsc",
        )

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
        filename=f"{gwfname}.chd",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiating MODFLOW 6 groundwater transport package
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
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 transport advection package
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 transport dispersion package
    if dispersivity != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            xt3d_off=True,
            alh=dispersivity,
            ath1=dispersivity,
            ktw=0.5918,
            kts=0.2700,
            filename=f"{gwename}.cnd",
        )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        heat_capacity_water=cpw,
        density_water=rhow,
        porosity=prsity,
        heat_capacity_solid=cps,
        density_solid=rhos,
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 transport constant concentration package
    flopy.mf6.ModflowGwectp(
        gwe,
        maxbound=len(ctpspd),
        stress_period_data=ctpspd,
        save_flows=False,
        pname="CTP-1",
        filename=f"{gwename}.ctp",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    flopy.mf6.ModflowGwessm(gwe, sources=[[]], filename=f"{gwename}.ssm")

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
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
        filename=f"{name}.gwfgwe",
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
        assert False, f'could not load concentration data from "{fpth}"'

    # This is the answer
    c_ans = [
        4.00000000e01,
        3.99999983e01,
        3.99999898e01,
        3.99999566e01,
        3.99998462e01,
        3.99995197e01,
        3.99986427e01,
        3.99964775e01,
        3.99915230e01,
        3.99809477e01,
        3.99597839e01,
        3.99198995e01,
        3.98488519e01,
        3.97288247e01,
        3.95359427e01,
        3.92403042e01,
        3.88070317e01,
        3.81985089e01,
        3.73777505e01,
        3.63125911e01,
        3.49801399e01,
        3.33708033e01,
        3.14911723e01,
        2.93652158e01,
        2.70334931e01,
        2.45504338e01,
        2.19800532e01,
        1.93907148e01,
        1.68496655e01,
        1.44180473e01,
        1.21469471e01,
        1.00748333e01,
        8.22648357e00,
        6.61329449e00,
        5.23470060e00,
        4.08034410e00,
        3.13261741e00,
        2.36924164e00,
        1.76562010e00,
        1.29679741e00,
        9.38944408e-01,
        6.70362685e-01,
        4.72056032e-01,
        3.27947150e-01,
        2.24829892e-01,
        1.52144844e-01,
        1.01654320e-01,
        6.70766201e-02,
        4.37223104e-02,
        2.81598160e-02,
        1.79249349e-02,
        1.12795213e-02,
        7.01828727e-03,
        4.31895689e-03,
        2.62924728e-03,
        1.58374083e-03,
        9.44125798e-04,
        5.57133590e-04,
        3.25507431e-04,
        1.88330495e-04,
        1.07925092e-04,
        6.12700035e-05,
        3.44648666e-05,
        1.92125906e-05,
        1.06157638e-05,
        5.81494908e-06,
        3.15821246e-06,
        1.70101068e-06,
        9.08679391e-07,
        4.81524218e-07,
        2.53159103e-07,
        1.32068539e-07,
        6.83748562e-08,
        3.51353218e-08,
        1.79225415e-08,
        9.07652498e-09,
        4.56413759e-09,
        2.27913640e-09,
        1.13033292e-09,
        5.56823550e-10,
        2.72491770e-10,
        1.32483548e-10,
        6.40015158e-11,
        3.07244529e-11,
        1.46584603e-11,
        6.95098705e-12,
        3.27643160e-12,
        1.53530190e-12,
        7.15261898e-13,
        3.31325318e-13,
        1.52616350e-13,
        6.99104644e-14,
        3.18504005e-14,
        1.44329547e-14,
        6.50576657e-15,
        2.91728603e-15,
        1.30145909e-15,
        5.77678170e-16,
        2.55141072e-16,
        1.12178999e-16,
        5.01900830e-17,
    ]

    msg = "gwe temperatures do not match stored concentrations"
    assert np.allclose(conc1[-1, 0, 0, :], c_ans, atol=1e-5), msg


# - No need to change any code below
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
