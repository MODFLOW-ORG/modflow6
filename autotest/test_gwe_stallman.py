# ## Test problem for GWE
#
# One-Dimensional Stallman Problem
# Compare MF6-GWE simulation results with analytical solution


# Imports
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

viscosity_on = [False]
cases = ["stallman"]

# Model units

length_units = "meters"
time_units = "seconds"

# Table MODFLOW 6 GWE parameters

nper = 600  # Number of periods
nstp = 6  # Number of time steps
perlen = 525600  # Simulation time length ($s$)
nlay = 120  # Number of layers
nrow = 1  # Number of rows
ncol = 1  # Number of columns
system_length = 60.0  # Length of system ($m$)
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv_str = "ranges from 0.1 to 1"  # Layer thickness
top = 60.0  # Top of the model ($m$)
hydraulic_conductivity = 1.0e-4  # Hydraulic conductivity ($m s^{-1}$)
porosity = 0.35  # Porosity (unitless)
alphal = 0.0  # Longitudinal dispersivity ($m$)
alphat = 0.0  # Transverse dispersivity ($m$)
diffc = 1.02882e-06  # Diffusion coefficient ($m s^{-1}$)
T_az = 10  # Ambient temperature ($^o C$)
dT = 5  # Temperature variation ($^o C$)
bulk_dens = 2630  # Bulk density ($kg/m^3$)
kd = 0.000191663  # Distribution coefficient (unitless)
ktw = 0.58
kts = 2
cpw = 4174.0
cps = 800.0
rhow = 1000.0
rhos = bulk_dens
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)

# Stress period input
per_data = []
for k in range(nper):
    per_data.append((perlen, nstp, 1.0))
per_mf6 = per_data

# Geometry input
tp = top
botm = []
for i in range(nlay):
    if i == 0:
        botm.append(59.9)
    elif i == 119:
        botm.append(0.0)
    else:
        botm.append(60 - i * 0.5)

# Head input
chd_data = {}
# for k in range(nper):
chd_data[0] = [[(0, 0, 0), 60.000000], [(119, 0, 0), 59.701801]]
chd_mf6 = chd_data

# Initial temperature input
strt_temp = T_az * np.ones((nlay, 1, 1), dtype=np.float32)

# Boundary temperature input
cnc_data = {}
for k in range(nper):
    cnc_temp = T_az + dT * np.sin(2 * np.pi * k * perlen / 365 / 86400)
    cnc_data[k] = [[(0, 0, 0), cnc_temp]]
cnc_mf6 = cnc_data

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-8, 1e-8, 0.97


# Analytical solution for Stallman analysis (Stallman 1965, JGR)
def Stallman(T_az, dT, tau, t, c_rho, darcy_flux, ko, c_w, rho_w, zbotm, nlay):
    zstallman = np.zeros((nlay, 2))
    K = np.pi * c_rho / ko / tau
    V = darcy_flux * c_w * rho_w / 2 / ko
    a = ((K**2 + V**4 / 4) ** 0.5 + V**2 / 2) ** 0.5 - V
    b = ((K**2 + V**4 / 4) ** 0.5 - V**2 / 2) ** 0.5
    for i in range(len(zstallman)):
        zstallman[i, 0] = zbotm[i]
        zstallman[i, 1] = (
            dT
            * np.exp(-a * (-zstallman[i, 0]))
            * np.sin(2 * np.pi * t / tau - b * (-zstallman[i, 0]))
            + T_az
        )
    return zstallman


#
# MODFLOW 6 (sim) flopy objects returned if building the model
#


def build_models(idx, test):
    # Base MF6 GWE model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    # sim_ws = os.path.join(ws, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=per_mf6, time_units=time_units)

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
        nogrb=True,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{gwfname}.dis",
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=0,
        k=hydraulic_conductivity,
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=top, filename=f"{gwfname}.ic")

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
        stress_period_data=chd_mf6,
        filename=f"{gwfname}.chd",
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # ----------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ----------------------------------

    # Instantiating MODFLOW 6 groundwater heat transport package
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
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(gwe, strt=strt_temp, filename=f"{gwename}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(gwe, scheme="TVD", filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        alh=alphal,
        ath1=alphat,
        ktw=ktw,
        kts=kts,
        pname="CND-1",
        filename=f"{gwename}.dsp",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=porosity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        filename=f"{gwename}.mst",
    )

    # Instantiating MODFLOW 6 transport constant concentration package
    flopy.mf6.ModflowGwectp(
        gwe,
        stress_period_data=cnc_mf6,
        pname="CTP-1",
        filename=f"{gwename}.tmp",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    flopy.mf6.ModflowGwessm(gwe, sources=[[]], filename=f"{gwename}.ssm")

    # Instantiate MODFLOW 6 heat transport output control package
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
    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwename = "gwe-" + name

    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    try:
        # load temperatures
        tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        times = tobj.get_times()
        sim_temps = tobj.get_data(totim=times[540])
    except:
        assert False, f'could not load concentration data from "{fpth}"'

    # Prepare to compare the results of MF6-GWE with analytical solution
    zconc = np.zeros((nlay, 2))
    for i in range(nlay):
        if i != (nlay - 1):
            zconc[i + 1, 0] = -(60 - botm[i])
        zconc[i, 1] = sim_temps[i][0][0]

    # Analytical solution - Stallman analysis
    tau = 365 * 86400
    t = 283824000.0
    c_w = 4174
    rho_w = 1000
    c_r = 800
    rho_r = 2630
    c_rho = c_r * rho_r * (1 - porosity) + c_w * rho_w * porosity
    darcy_flux = 5.00e-07
    ko = 1.503
    analytical_temps = Stallman(
        T_az, dT, tau, t, c_rho, darcy_flux, ko, c_w, rho_w, zconc[:, 0], nlay
    )

    # plt.plot(zconc[:, 1], zconc[:, 0], "k--", linewidth=0.5, label="MF6-GWE")
    # plt.plot(zanal[:, 1], zanal[:, 0], "bo", mfc="none", label="Analytical")
    # plt.xlim(T_az - dT, T_az + dT)
    # plt.ylim(-top, 0)
    # plt.ylabel("Depth (m)")
    # plt.xlabel("Temperature (deg C)")
    # plt.legend()
    # plt.savefig("stallman.png")

    msg = "gwe temperatures do not match analytical temperatures"
    assert np.allclose(zconc[:, 1], analytical_temps[:, 1], atol=1e-1), msg

    return


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
