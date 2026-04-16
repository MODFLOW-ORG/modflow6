import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from scipy import optimize

cases = ["mst08a", "mst08b", "mst08c"]
sorption_type = ["linear", "freundlich", "langmuir"]

# Test description
#
# This test builds and runs a single-cell with Sorption and a Source.
# The MST package is tested with several equilibrium isotherm.
#
# We inject a known total mass over a finite time and validate the
# end-of-period concentration against the **analytical mass-balance**
# solution: $$M_{total} = \theta V C - f_m \rho_b V \overline{C}$$.
#
# This solution doesn't have an exact closed form for all isotherms,
# so we solve it algebraically using a root-finding method.

# Parameters
## Time parameters
nper = 1
perlen = [10.0]
nstp = [10]
tsmult = [1.0]

## Solver parameters
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0

## Model domain and grid definition
nlay, nrow, ncol = 1, 1, 1
delr, delc = 1.0, 1.0
top = 1.0
botm = 0.0

## Mass source parameters
Msrc = 0.05  # mass source term (M/T)

## MST parameters
porosity = 0.25  # porosity (dimensionless)
fm = 1  # fraction of mass in the solid phase (dimensionless)
rho_bulk = 1700.0  # bulk density of the solid phase (M/L^3)

## Sorption parameters
### Linear sorption parameters
Kd = 5.0e-4  # linear sorption coefficient (L^3/M)
# Freundlich sorption parameters
Kf = 1.0e-4  # Freundlich sorption coefficient (L^3M^(-1) )**n
n_exp = 0.8  # Freundlich exponent (dimensionless)
### Langmuir sorption parameters
Kl = 5.0  # Langmuir sorption coefficient (L^3/M)
Sbar = 1e-4  # Langmuir maximum solid phase concentration


# Algebraic solution for equilibrium concentration
# solved using a root-finding method
def algebraic_solution(sorption_type):
    Mtot = Msrc * perlen[0]  # total mass added over the period (M)
    V = delr * delc * (top - botm)  # volume of the cell (L^3)

    C0 = Mtot / (porosity * V)

    # Find exact concentration at equilibrium for the listed isotherm
    def C_equilibrium_linear():
        def f(C):
            return porosity * V * C + fm * rho_bulk * Kd * C * V - Mtot

        def fprime(C):
            return porosity * V + fm * rho_bulk * Kd * V

        c_exact = optimize.newton(f, C0, fprime=fprime)
        return c_exact

    # Find exact concentration at equilibrium for the Freundlich  isotherm
    def C_equilibrium_freundlich():
        def f(C):
            return porosity * V * C + fm * rho_bulk * Kf * C**n_exp * V - Mtot

        def fprime(C):
            return porosity * V + fm * rho_bulk * Kf * n_exp * C ** (n_exp - 1) * V

        c_exact = optimize.newton(f, C0, fprime=fprime)
        return c_exact

    # Find exact concentration at equilibrium for the Langmuir isotherm
    def C_equilibrium_langmuir():
        def f(C):
            return (
                porosity * V * C
                + fm * rho_bulk * (Kl * Sbar * C) / (1 + Kl * C) * V
                - Mtot
            )

        def fprime(C):
            return porosity * V + fm * rho_bulk * (Kl * Sbar) / ((1 + Kl * C) ** 2) * V

        c_exact = optimize.newton(f, C0, fprime=fprime)
        return c_exact

    if sorption_type == "linear":
        c_exact = C_equilibrium_linear()
    elif sorption_type == "freundlich":
        c_exact = C_equilibrium_freundlich()
    elif sorption_type == "langmuir":
        c_exact = C_equilibrium_langmuir()
    else:
        raise ValueError(f"Unknown sorption type: {sorption_type}")

    return c_exact


def build_models(idx, test):
    name = cases[idx]
    sorption = sorption_type[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")

    # create tdis package
    tdis_rc = [perlen[0], nstp[0], tsmult[0]]
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwt model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name)

    # create iterative model solution
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
    )

    # create discretization package
    dis = flopy.mf6.ModflowGwtdis(
        gwt, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # source term
    spd = {0: [((0, 0, 0), Msrc)]}
    src = flopy.mf6.ModflowGwtsrc(gwt, stress_period_data=spd)

    # create mst package
    if sorption == "linear":
        mst = flopy.mf6.ModflowGwtmst(
            gwt,
            sorption="linear",
            porosity=porosity,
            bulk_density=rho_bulk,
            distcoef=Kd,
        )
    elif sorption == "freundlich":
        mst = flopy.mf6.ModflowGwtmst(
            gwt,
            sorption="freundlich",
            porosity=porosity,
            bulk_density=rho_bulk,
            distcoef=Kf,
            sp2=n_exp,
        )
    elif sorption == "langmuir":
        mst = flopy.mf6.ModflowGwtmst(
            gwt,
            sorption="langmuir",
            porosity=porosity,
            bulk_density=rho_bulk,
            distcoef=Kl,
            sp2=Sbar,
        )
    else:
        raise ValueError(f"Unknown sorption type: {sorption}")

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


# Algebraic solution for equilibrium concentration


def check_output(idx, test):
    sorption = sorption_type[idx]
    c_exact = algebraic_solution(sorption)

    name = cases[idx]
    fpth = os.path.join(test.workspace, f"{name}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data(idx=-1)

    assert np.isclose(conc[0, 0, 0], c_exact, rtol=1e-6), "Concentrations do not match!"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
