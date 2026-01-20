"""Uses general-head and drain boundaries on the left and right
sides of the model domain, respectively, to drive flow from left to
right.  Tests that head-dependent boundary conditions are properly
accounting for viscosity when VSC is active.  Similar to gwf-vsc01-bnd
but employs head-dependent boundary on the left and right side of the
model
"""

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Setup scenario input
cases = ["no-vsc02-bnd", "vsc02-bnd", "no-vsc02-k"]
hyd_cond = [1205.49396942506, 864.0]  # Hydraulic conductivity (m/d)
viscosity_on = [False, True, False]
hydraulic_conductivity = [hyd_cond[0], hyd_cond[1], hyd_cond[1]]

# Model units

length_units = "cm"
time_units = "seconds"

# Table of model parameters

nper = 1  # Number of periods
nstp = 500  # Number of time steps
perlen = 0.5  # Simulation time length ($d$)
nlay = 1  # Number of layers
nrow = 10  # Number of rows
ncol = 80  # Number of columns
system_length = 2.0  # Length of system ($m$)
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delv = 1.0  # Layer thickness
top = 1.0  # Top of the model ($m$)
initial_temperature = 35.0  # Initial temperature (unitless)
porosity = 0.26  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rho_water = 1000  # Density of water ($kg/m^3$)
rho_solids = 2650  # Density of the aquifer material ($kg/m^3$)
C_p_w = 4180  # Heat Capacity of water ($J/kg/C$)
C_s = 880  # Heat capacity of the solids ($J/kg/C$)
D_m = K_therm / (porosity * rho_water * C_p_w)
rhob = (1 - porosity) * rho_solids  # Bulk density ($kg/m^3$)
K_d = C_s / (rho_water * C_p_w)  # Partitioning coefficient ($m^3/kg$)
inflow = 5.7024  # ($m^3/d$)

botm = [top - k * delv for k in range(1, nlay + 1)]

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-6, 0.97


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    print(f"Building model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating time discretization
    tdis_ds = ((perlen, nstp, 1.0),)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_ds, time_units=time_units)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    # Instantiating solver
    ims = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])

    # Instantiating DIS
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
    )

    # Instantiating NPF
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=0,
        k=hydraulic_conductivity[idx],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    # Instantiating VSC
    if viscosity_on[idx]:
        # Instantiate viscosity (VSC) package
        vsc_filerecord = f"{gwfname}.vsc.bin"
        vsc_pd = [(0, 0.0, 20.0, gwtname, "temperature")]
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

    # Instantiating GHB
    ghbcond = hydraulic_conductivity[idx] * delv * delc / (0.5 * delr)
    ghbspd = [[(0, i, 0), top + 3, ghbcond, initial_temperature] for i in range(nrow)]
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        pname="GHB-1",
        auxiliary="temperature",
    )

    # Instantiating DRN
    drnspd = [
        [(0, i, ncol - 1), top, 1.2 * ghbcond, initial_temperature] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfdrn(
        gwf,
        stress_period_data=drnspd,
        pname="DRN-1",
        auxiliary="temperature",
    )

    # Instatiatingi OC
    head_filerecord = f"{gwfname}.hds"
    budget_filerecord = f"{gwfname}.bud"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Setup the GWT model for simulating heat transport
    # -------------------------------------------------
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    # Instantiating solver for GWT
    imsgwt = flopy.mf6.ModflowIms(
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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwtname])

    # Instantiating DIS for GWT
    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Instantiating MST for GWT
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        sorption="linear",
        bulk_density=rhob,
        distcoef=K_d,
        pname="MST-1",
        filename=f"{gwtname}.mst",
    )

    # Instantiating IC for GWT
    flopy.mf6.ModflowGwtic(gwt, strt=initial_temperature)

    # Instantiating ADV for GWT
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # Instantiating DSP for GWT
    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=D_m)

    # Instantiating SSM for GWT
    sourcerecarray = [
        ("GHB-1", "AUX", "TEMPERATURE"),
        ("DRN-1", "AUX", "TEMPERATURE"),
    ]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # Instantiating OC for GWT
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # Instantiating GWF/GWT Exchange
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return sim, None


def check_output(idx, test):
    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name

    fname = gwfname + ".bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    budobj = flopy.utils.CellBudgetFile(fname, precision="double")
    outbud = budobj.get_data(text="             GHB")

    # Establish known answer:
    stored_ans = 452.5316256451224

    if idx == 0:
        no_vsc_bud_last = np.array(outbud[-1].tolist())
        sim_val_1 = no_vsc_bud_last[:, 2].sum()

        # Ensure latest simulated value hasn't changed from stored answer
        assert np.allclose(sim_val_1, stored_ans, atol=1e-4), (
            "Flow in the " + cases[0] + " test problem (doesn't simulate "
            "viscosity) has changed,\n should be "
            + str(stored_ans)
            + " but instead is "
            + str(sim_val_1)
        )

    elif idx == 1:
        with_vsc_bud_last = np.array(outbud[-1].tolist())
        sim_val_2 = with_vsc_bud_last[:, 2].sum()

        # Ensure latest simulated value hasn't changed from stored answer
        assert np.allclose(sim_val_2, stored_ans, atol=1e-4), (
            "Flow in the "
            + cases[1]
            + " test problem (simulates viscosity) has changed,\n should be "
            + str(stored_ans)
            + " but instead is "
            + str(sim_val_2)
        )

    elif idx == 2:
        no_vsc_low_k_bud_last = np.array(outbud[-1].tolist())
        sim_val_3 = no_vsc_low_k_bud_last[:, 2].sum()

        # Ensure the flow leaving model 3 is less than what leaves model 2
        assert abs(stored_ans) > abs(sim_val_3), (
            "Exit flow from model "
            + cases[1]
            + " should be greater than flow existing "
            + cases[2]
            + ", but it is not."
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
