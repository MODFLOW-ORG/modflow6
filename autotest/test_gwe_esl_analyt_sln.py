r"""
An analytical solution provided by Carslaw & Jaeger (1947) and discussed in
accompanying Techniques & Methods report.

Energy is added to the right hand side boundary using the energy source loading
(ESL) package.  Basic model set up is below, with a slab of unit thickness
(1.0 m) that is 100 m "deep" with energy being loaded on right side.
Temperature will begin to rise on the right and propagate to the left. There are
no sinks in this first example.  The titles that follow, for example
"Section 43, case x" refer to specific analytical solutions found in Carslaw &
Jaeger (1947)

Section 43, case i:
-------------------

       | <--------------------------   10 m   --------------------------> |

       +------------------------------------------------------------------+
       |                      Initial temperature = T_0                   | <-- *ESL
       +------------------------------------------------------------------+
       ^                                      * ESL: Energy Source Loading Boundary
       |
       No heat-flow boundary


Section 43, case ii:
--------------------

       +------------------------------------------------------------------+
       |                     Initial temperature = 0.0                    | <-- *ESL
       +------------------------------------------------------------------+
       ^
       |
       Specified temperature boundary, T_0


Section 43, case iii:
---------------------

        +------------------------------------------------------------------+
CTP ->  |                                                                  | <- CTP = T_0
  = T_0 +------------------------------------------------------------------+
         \-------------------------------------------------------------/
                                   |
              Uniform, constant heat production throughout the slab


  Specified temperature boundary, T_0

"""  # noqa

import math
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Parameters that vary by scenario
cases = ["esl_casei", "esl_caseii", "esl_caseiii"]
perlen = {0: [100, 900], 1: [100, 9900], 2: [100, 900]}
nstp = [100, 900]
tsmult = [
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
]

T_0 = 0.0  # Initial temperature in all scenarios.
# Additionally serves as the CTP bnd temperature in scenarios ii and iii

xt3d = [True]

# Parameters for tdis package

# Model parameters

# Slab thickness
el = 10.0  # meters

# Cell dimensions
nlay, nrow, ncol = 1, 1, 1000
delc = delz = 1.0
delr = el / ncol

top = 1.0
laytyp = 1
strt = 0.0
ss = 0.0
sy = 0.1
botm = [0.0]
strt = 0.0
hnoflo = 1e30
hdry = -1e30
hk = 1.0
alh = 0.0
alv = 0.0
ath1 = 0.0
atv = 0.0
lhv = 2454.0

# Solver parameters
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0

# Boundary condition (same for all scenarios
chd_perdat = {0: [[(0, 0, 0), 0.0000000], [(0, 0, ncol - 1), 0.0000000]]}

# The following lists the parameter used in generating an analytical solution
# for comparing the MF6 solution to.

# Density of the solids
rhos = 2700.0  # kg/m^3

# Density of water
rhow = 1000.0  # kg/m^3

# Heat capacity of the solids
Cps = 703.7  # J / (kg * C)

# Heat capacity of water
Cpw = 4183.0  # J / (kg * C)

# "Bulk" thermal conductivity
# For this problem, K_t_bulk represents the thermal conductivity of the
# solid material only since the problem represents a dry slab
Kts = 0.2700 * 86400  # * 1e8  # J / (day * m * C)
Ktw = 0.5918 * 86400  # J / (day * m * C)

# Amount of saturation in a cell (cells are dry in this example)
Sw = 0.0  # dimensionless

# Define porosity
theta = 0.1  # dimensionless

# Equation 4-6 (for now anyway)
K_t_bulk = Sw * theta * Ktw + (1 - theta) * Kts

# Eqn 7-4: Bulk specific heat on a per volume basis
rho_C_bulk = Sw * theta * rhow * Cpw + (1 - theta) * rhos * Cps

# Eqn 7-3: Bulk thermal diffusivity
D = K_t_bulk / rho_C_bulk


# Energy input to boundary (q_x term in the documentation)
def calc_ener_input(primer_val):
    ener_add_rate = delr * delc * delz * rho_C_bulk * primer_val
    return ener_add_rate


# Define function to solve analytical solution
# Function names derive from equation numbers in Techniques and Methods


def build_models(idx, test, ener_input):
    name = cases[idx]

    # Build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # Create tdis package
    tdis_rc = []
    for i in range(len(perlen[idx])):
        tdis_rc.append((perlen[idx][i], nstp[i], tsmult[idx][i]))

    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=len(perlen[idx]), perioddata=tdis_rc
    )

    # Create GWF model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = True

    # Create iterative model solution and register the gwf model with it
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

    # Discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        pname="DIS-GWF",
        filename=f"{gwfname}.dis",
    )

    # Initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=strt, pname="IC-HD", filename=f"{gwfname}.ic")

    # Node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
        filename=f"{gwfname}.npf",
    )

    # Instantiating MODFLOW 6 storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=ss,
        sy=sy,
        iconvert=1,
        steady_state=False,
        transient=True,
        pname="STO",
        filename=f"{gwfname}.sto",
    )

    # Constant head files
    chd = None
    #    chd = flopy.mf6.ModflowGwfchd(
    #        gwf,
    #        maxbound=len(chd_perdat[0]),
    #        stress_period_data=chd_perdat,
    #        save_flows=False,
    #        pname="CHD-1",
    #        filename=f"{gwfname}.chd"
    #    )

    # Output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # -----------------
    # Create GWE model
    # -----------------

    gwename = "gwe_" + name
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    gwe.name_file.save_flows = True

    # Create iterative model solution and register the gwt model with it
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

    # Initial conditions
    flopy.mf6.ModflowGweic(gwe, strt=T_0, pname="IC-1", filename=f"{gwename}.ic")

    # Advection
    flopy.mf6.ModflowGweadv(
        gwe, scheme="UPSTREAM", pname="ADV-E", filename=f"{gwename}.adv"
    )

    # Heat conduction
    flopy.mf6.ModflowGwecnd(
        gwe,
        alh=alh,
        alv=alv,
        ath1=ath1,
        atv=atv,
        ktw=Ktw,
        kts=Kts,
        pname="CND-1",
        filename=f"{gwename}.cnd",
    )

    flopy.mf6.ModflowGweest(
        gwe,
        porosity=theta,
        heat_capacity_water=Cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=Cps,
        density_solid=rhos,
        pname="EST-1",
        filename=f"{gwename}.est",
    )

    # Constant temperature
    # Note: Implementation of the CTP boundary depends on which analytical sln
    #       is in view
    #       See notes at top of script regarding scenarios
    if idx > 0:
        if idx == 1:
            ctp = {0: [[(0, 0, 0), T_0]]}
        elif idx == 2:
            ctp = {0: [[(0, 0, 0), T_0], [(0, 0, ncol - 1), T_0]]}
        flopy.mf6.ModflowGwectp(
            gwe,
            maxbound=len(ctp),
            stress_period_data=ctp,
            save_flows=True,
            pname="CTP-1",
            filename=f"{gwename}.ctp",
        )

    # Instantiate energy source loading (ESL) package
    if idx < 2:
        esrc = {0: [[(0, 0, ncol - 1), ener_input]]}

    elif idx == 2:
        esrcs = []
        for j in np.arange(ncol):
            esrcs.append([(0, 0, j), ener_input])
        esrc = {0: esrcs}

    flopy.mf6.ModflowGweesl(
        gwe,
        maxbound=len(esrc[0]),
        stress_period_data=esrc,
        save_flows=False,
        pname="ESL-1",
        filename=f"{gwename}.esl",
    )

    # Sources
    if chd is not None:
        flopy.mf6.ModflowGwessm(
            gwe, sources=[[]], pname="SSM-E", filename=f"{gwename}.ssm"
        )

    # Output control
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    # GWF GWE exchange
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )

    return sim, None


def eq7_24(x, t, l, D, T_0, ener_add_rate):
    # Compute corresponding x_hat term
    x_hat = x / l  # Dimensionless distance

    # Compute corresponding t_hat term
    t_hat = D * t / l**2  # Dimensionless time

    # Solve equation 7-24
    term1 = (1 / 2) * (x_hat**2 - 1 / 3)
    summation_terms = [
        ((-1) ** n)
        / n**2
        * math.exp(-1 * n**2 * math.pi**2 * t_hat)
        * math.cos(n * math.pi * x_hat)
        for n in np.arange(1, 1000)
    ]
    term2 = 2 / (math.pi**2) * np.sum(summation_terms)
    T = T_0 + ener_add_rate * l / K_t_bulk * (t_hat + term1 - term2)

    return T


def eq7_25(x, t, l, D, T_0, ener_add_rate):
    # Compute corresponding x_hat term
    x_hat = x / l  # Dimensionless distance

    # Compute corresponding t_hat term
    t_hat = D * t / l**2  # Dimensionless time

    # Solve equation 7-25
    summation_terms = [
        ((-1) ** n)
        / (2 * n + 1) ** 2
        * math.exp(-1 * (2 * n + 1) ** 2 * math.pi**2 * t_hat / 4)
        * math.sin((2 * n + 1) * math.pi * x_hat / 2)
        for n in np.arange(0, 1000)
    ]
    term1 = (8 / math.pi**2) * np.sum(summation_terms)

    T = T_0 + ener_add_rate * l / K_t_bulk * (x_hat - term1)

    return T


def eq7_26(x, t, el, D, T_0, ener_add_rate):
    # Compute corresponding x_hat term
    x_hat = x / el  # Dimensionless distance

    # Compute corresponding t_hat term
    t_hat = D * t / el**2  # Dimensionless time

    # Solve equation 7-26
    term1 = x_hat * (1 - x_hat)
    summation_terms = [
        1
        / (2 * n + 1) ** 3
        * math.exp(-1 * (2 * n + 1) ** 2 * math.pi**2 * t_hat)
        * math.sin((2 * n + 1) * math.pi * x_hat)
        for n in np.arange(0, 1000)
    ]
    term2 = (8 / math.pi**3) * np.sum(summation_terms)
    T = T_0 + 0.5 * ener_add_rate * el**2 / K_t_bulk * (term1 - term2)

    return T


def check_output(idx, test, ener_input):
    name = test.name
    gwename = "gwe_" + name

    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
    sim_temps = tobj.get_alldata()

    # Compare simulated output to analytical solutions (scenario dependent)
    if idx < 2:
        t_accumulate = 0.0
        area_input = delc * delz
        ener_flux = ener_input / area_input

        for sp, t in enumerate(perlen[idx]):
            # Time to solution
            t_accumulate += t  # days

            # Iterate over x which represents the cell centroid locations for which
            # analytical solution is desired (cell centroid locations)
            cell_centroids = []
            analytical_temps = []
            for x in np.arange(delr / 2, el + (delr / 2), delr):
                cell_centroids.append(x)
                if idx == 0:
                    T = eq7_24(x, t_accumulate, el, D, T_0, ener_flux)
                elif idx == 1:
                    T = eq7_25(x, t_accumulate, el, D, T_0, ener_flux)
                analytical_temps.append(T)

            analytical_temps = np.array(analytical_temps)
            assert np.allclose(analytical_temps, sim_temps[sp, 0, 0, :], atol=0.005), (
                "simulated solution is whacked"
            )
            # plt.plot(cell_centroids, analytical_temps, "r-", label="Analytical Solution")  # noqa
            # plt.plot(cell_centroids, sim_temps[sp, 0, 0, :], "b--", label="GWE")
            # plt.axhline(0.0, color='black')
            # plt.legend()
            # plt.show()

    elif idx == 2:
        t_accumulate = 0.0
        ener_src = ener_input / (delr * delc * delz)

        for sp, t in enumerate(perlen[idx]):
            # Time to solution
            t_accumulate += t  # days

            analytical_temps = []
            cell_centroids = []
            for x in np.arange(delr / 2, el + (delr / 2), delr):
                cell_centroids.append(x)
                T = eq7_26(x, t_accumulate, el - (delr * 2), D, T_0, ener_src)
                analytical_temps.append(T)

            analytical_temps = np.array(analytical_temps)
            if sp == 0:
                atol = 0.16
            else:
                atol = 0.47

            assert np.allclose(analytical_temps, sim_temps[sp, 0, 0, :], atol=atol), (
                "simulated solution is whacked"
            )

            # plt.plot(cell_centroids, analytical_temps, "r-", label="Analytical Solution")  # noqa
            # plt.plot(cell_centroids, sim_temps[sp, 0, 0, :], "b--", label="GWE")
            # plt.axhline(0.0, color='black')
            # plt.legend()
            # plt.show()


# - No need to change any code below
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    if idx < 2:
        ener_input = calc_ener_input(1.0)
    elif idx == 2:
        ener_input = calc_ener_input(0.1)
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t, ener_input),
        check=lambda t: check_output(idx, t, ener_input),
    )
    test.run()
