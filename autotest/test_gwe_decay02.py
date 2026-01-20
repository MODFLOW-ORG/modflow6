"""
Test problem for decay of energy in EST package of GWE.  Compares energy loss
to that of an ESL boundary
"""

# Imports

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

cases = ["decay"]

# Model units
length_units = "meters"
time_units = "seconds"

rho_w = 1000  # kg/m^3
rho_s = 2500  # kg/m^3
n = 0.2  # -
gamma_w = -1000  # J/s/m^3, arbitrary value for zero-order aqueous heat production
gamma_s = -0.1  # J/s/kg, arbitrary value for zero-order solid heat production
c_w = 4000  # J/kg/degC
c_s = 1000  # J/kg/decC
T0 = 0  # degC

nrow = 1
ncol = 1
nlay = 1
delr = 1  # m
delc = 1  # m
top = 1  # m
botm = 0  # m

perlen = 86400  # s
nstp = 20


def add_gwe(sim, gwename, add_esl=False):
    gwe = flopy.mf6.ModflowGwe(
        sim,
        modelname=gwename,
        save_flows=True,
        model_nam_file=f"{gwename}.nam",
    )
    dis = flopy.mf6.ModflowGwedis(
        gwe, nrow=nrow, ncol=ncol, nlay=nlay, delr=delr, delc=delc, top=top, botm=botm
    )
    ic = flopy.mf6.ModflowGweic(gwe, strt=T0)
    if not add_esl:
        zero_order_decay_water = True
        zero_order_decay_solid = True
    else:
        zero_order_decay_water = False
        zero_order_decay_solid = False

        esl_amt = n * gamma_w + (1 - n) * gamma_s * rho_s
        esl_spd = {
            0: [[(0, 0, 0), -esl_amt]],
        }
        esl = flopy.mf6.ModflowGweesl(
            gwe,
            stress_period_data=esl_spd,
            pname="ESL",
            filename=f"{gwename}.esl",
        )

    est = flopy.mf6.ModflowGweest(
        gwe,
        zero_order_decay_water=zero_order_decay_water,
        zero_order_decay_solid=zero_order_decay_solid,
        density_water=rho_w,
        density_solid=rho_s,
        heat_capacity_water=c_w,
        heat_capacity_solid=c_s,
        porosity=n,
        decay_water=gamma_w,
        decay_solid=gamma_s,
    )

    oc = flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwe.name}.bud",
        temperature_filerecord=f"{gwe.name}.ucn",
        printrecord=[("BUDGET", "ALL"), ("TEMPERATURE", "ALL")],
        saverecord=[("BUDGET", "ALL"), ("TEMPERATURE", "ALL")],
    )

    return sim


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]
    gwename = "gwe-" + name

    print(f"Building MF6 model...{name}")

    sim = flopy.mf6.MFSimulation(
        sim_name="heat",
        sim_ws=ws,
        exe_name="mf6",
        version="mf6",
    )
    tdis = flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(perlen, nstp, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim, complexity="SIMPLE", inner_dvclose=0.001
    )  # T can not become negative in this model

    # add first GWE model
    sim = add_gwe(sim, gwename + "-1", add_esl=False)
    # add second GWE model
    sim = add_gwe(sim, gwename + "-2", add_esl=True)

    return sim, None


def check_output(idx, test):
    print("evaluating results...")
    ws = test.workspace

    msg = (
        "Differences detected between the simulated results for zeroth-order "
        "energy decay and the ESL Package.  The respective approaches are "
        "expected to give the same answer."
    )

    # read transport results from GWE model
    name = cases[idx]
    gwename1 = "gwe-" + name + "-1"
    gwename2 = "gwe-" + name + "-2"

    # Get the MF6 temperature output
    sim = test.sims[0]
    gwe1 = sim.get_model(gwename1)
    temp_ts1 = gwe1.output.temperature().get_ts((0, 0, 0))
    t1 = temp_ts1[:, 0]

    gwe2 = sim.get_model(gwename2)
    temp_ts2 = gwe2.output.temperature().get_ts((0, 0, 0))
    t2 = temp_ts2[:, 0]

    assert np.isclose(temp_ts1[-1, 1], temp_ts2[-1, 1], atol=1e-10), msg


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
