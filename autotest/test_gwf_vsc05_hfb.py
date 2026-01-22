"""
Uses constant head and general-head boundaries on the left and right
sides of a 10 row by 10 column by 1 layer model to drive flow from left to
right.  Tests that a hydraulic flow barrier accounts for changes in
viscosity when temperature is simulated. Barrier is between middle two
columns, but only cuts across the bottom 5 rows.
Model 1: VSC inactive, uses a higher speified K that matches what the VSC
         package will come up with
Model 2: VSC active, uses a lower K so that when VSC is applied, resulting
         K's match model 1 and should result in the same flows across the
         model domain
Model 3: VSC inactive, uses the lower K of model 2 and checks that flows
         in model 3 are indeed lower than in model 2 when turning VSC off.
         Model simulates hot groundwater with lower viscosity resulting in
         more gw flow through the model domain.Flows that are checked are
         the row-wise flows between columns 5 and 6 (e.g., cell 5 to 6, 15
         to 16, etc.)
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["no-vsc05-hfb", "vsc05-hfb", "no-vsc05-k", "hfb-idomain"]
hyd_cond = [1205.49396942506, 864.0]  # Hydraulic conductivity (m/d)
viscosity_on = [False, True, False, False]
hydraulic_conductivity = [hyd_cond[0], hyd_cond[1], hyd_cond[1], hyd_cond[0]]
icell_inactive = 3

# Model units

length_units = "centimeters"
time_units = "seconds"

# Table of model parameters

nper = 1  # Number of periods
nstp = 10  # Number of time steps
perlen = 10  # Simulation time length ($d$)
nlay = 1  # Number of layers
nrow = 10  # Number of rows
ncol = 10  # Number of columns
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
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
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

    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    if idx == icell_inactive:
        sim.simulation_data.verify_data = False
        idomain[0, 7, 4] = 0

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
        idomain=idomain,
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

    # Instantiating CHD (leftside, "inflow" boundary)
    chdspd = [[(0, i, 0), 2.0, initial_temperature] for i in range(nrow)]
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        pname="CHD-1",
        auxiliary="temperature",
    )

    # Instantiating GHB (rightside, "outflow" boundary)
    ghbcond = hydraulic_conductivity[idx] * delv * delc / (0.5 * delr)
    ghbspd = [
        [(0, i, ncol - 1), top, ghbcond, initial_temperature] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        pname="GHB-1",
        auxiliary="temperature",
    )

    # Instantiate Hydraulic Flow-Barrier (HFB) package
    # Barrier present between middle two columns of the model domain, but only
    # in rows 6-10.  Remember that the hydraulic characteristic is the barrier
    # hydraulic conductivity divided by the width of the horizontal-flow
    # barrier.  Assuming a barrier width of 10 cm (0.1 m) and desire to have
    # the barrier's K be 1/10th of the aquifer hydraulic conductivity.
    hfbspd = []
    K = 0.1 * hydraulic_conductivity[idx]
    for i in np.arange(5, 10, 1):
        hfbspd.append(((0, i, 4), (0, i, 5), K))
    flopy.mf6.ModflowGwfhfb(
        gwf,
        print_input=True,
        maxhfb=len(hfbspd),
        stress_period_data=hfbspd,
        pname="HFB-1",
    )

    # Instantiating OC
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
        idomain=idomain,
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
        ("CHD-1", "AUX", "TEMPERATURE"),
        ("GHB-1", "AUX", "TEMPERATURE"),
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
    sim1 = flopy.mf6.MFSimulation.load(sim_ws=test.workspace, load_only=["dis"])
    gwf = sim1.get_model(gwfname)

    # Get grid data
    grdname = gwfname + ".dis.grb"
    bgf = flopy.mf6.utils.MfGrdFile(os.path.join(test.workspace, grdname))
    ia, ja = bgf.ia, bgf.ja

    fname = gwfname + ".bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    budobj = flopy.utils.CellBudgetFile(fname, precision="double")
    outbud = budobj.get_data(text="    FLOW-JA-FACE")[-1].squeeze()

    # Establish known answer for the "with viscosity" variant:
    stored_ans = np.array(
        [
            [4, 5, 131.03196884892344],
            [14, 15, 133.1834658429856],
            [24, 25, 139.31716925610493],
            [34, 35, 156.14497435040056],
            [44, 45, 209.1055337693415],
            [54, 55, 36.91267872240113],
            [64, 65, 46.16474722642168],
            [74, 75, 51.2708505192076],
            [84, 85, 54.04369740428511],
            [94, 95, 55.27469944201896],
        ]
    )

    # Look at flow entering the left face for the cells in the 6th (1-based) column
    cells = [gwf.modelgrid.get_node([(0, i, 5)])[0] for i in np.arange(nrow)]

    vals_to_store = []  # Will always have 10 vals, 1 per row

    # Note that the layer, row, column indices will be zero-based
    for celln in cells:
        for ipos in range(ia[celln] + 1, ia[celln + 1]):
            cellm = ja[ipos]
            if cellm == celln - 1:
                vals_to_store.append([cellm, celln, outbud[ipos]])

    if idx == 0:
        no_vsc_bud_last = np.array(vals_to_store)

        # Ensure with and without VSC simulations give nearly identical flow results
        # for each cell-to-cell exchange between columns 5 and 6
        assert np.allclose(no_vsc_bud_last[:, 2], stored_ans[:, 2], atol=1e-3), (
            "Flow in models "
            + cases[0]
            + " and the established answer should be approximately "
            "equal, but are not."
        )

    elif idx == 1:
        with_vsc_bud_last = np.array(vals_to_store)

        assert np.allclose(with_vsc_bud_last[:, 2], stored_ans[:, 2], atol=1e-3), (
            "Flow in models "
            + cases[1]
            + " and the established answer should be approximately "
            "equal, but are not."
        )

    elif idx == 2:
        no_vsc_low_k_bud_last = np.array(vals_to_store)

        # Ensure the cell-to-cell flow between columns 5 and 6 in model
        # 3 is less than what's in the "with viscosity" model
        assert np.less(no_vsc_low_k_bud_last[:, 2], stored_ans[:, 2]).all(), (
            "Exit flow from model the established answer "
            "should be greater than flow existing " + cases[2] + ", but it is not."
        )
    elif idx == 3:
        cell_inactive_vals = np.array(vals_to_store)

        # Ensure HFB inactive cell (1,8,5) cell-to-cell flows
        # have been excluded
        assert np.allclose(
            cell_inactive_vals[:, 0], [4, 14, 24, 34, 44, 54, 64, 84, 94]
        )
        assert np.allclose(
            cell_inactive_vals[:, 1], [5, 15, 25, 35, 45, 55, 65, 85, 95]
        )


# - No need to change any code below
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
