"""
Test problem for decay of energy in EST package of GWE.  Uses a single-cell
model.  Test contributed by Cas Neyens.
"""

# Imports

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

cases = ["decay-aqe", "decay-sld", "decay-both"]

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

parameters = {
    # aqueous
    "decay-aqe": {
        "zero_order_decay_water": True,
        "zero_order_decay_solid": False,
        "decay_water": gamma_w,
        "decay_solid": gamma_s,
    },
    # solid
    "decay-sld": {
        "zero_order_decay_water": False,
        "zero_order_decay_solid": True,
        "decay_water": gamma_w,
        "decay_solid": gamma_s,
    },
    # combined
    "decay-both": {
        "zero_order_decay_water": True,
        "zero_order_decay_solid": True,
        "decay_water": gamma_w,
        "decay_solid": gamma_s,
    },
}


def temp_z_decay(t, rho_w, rho_s, c_w, c_s, gamma_w, gamma_s, n, T0):
    t = np.atleast_1d(t)
    coeff = (-gamma_w * n - gamma_s * (1 - n) * rho_s) / (
        rho_w * c_w * n + rho_s * c_s * (1 - n)
    )
    return coeff * t + T0


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]
    gwename = "gwe-" + name

    decay_kwargs = parameters[name]

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
    est = flopy.mf6.ModflowGweest(
        gwe,
        zero_order_decay_water=decay_kwargs["zero_order_decay_water"],
        zero_order_decay_solid=decay_kwargs["zero_order_decay_solid"],
        density_water=rho_w,
        density_solid=rho_s,
        heat_capacity_water=c_w,
        heat_capacity_solid=c_s,
        porosity=n,
        decay_water=decay_kwargs["decay_water"],
        decay_solid=decay_kwargs["decay_solid"],
    )

    oc = flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwe.name}.bud",
        temperature_filerecord=f"{gwe.name}.ucn",
        printrecord=[("BUDGET", "ALL"), ("TEMPERATURE", "ALL")],
        saverecord=[("BUDGET", "ALL"), ("TEMPERATURE", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")
    msg = (
        "Differences detected between the simulated results for zeroth-order "
        "energy decay and the expected solution for decay specified in "
    )
    msg0 = msg + "the aqueous phase."
    msg1 = msg + "the solid phase."
    msg2 = msg + "both the aqueous and solid phases."

    # read transport results from GWE model
    name = cases[idx]
    gwename = "gwe-" + name

    # Get the MF6 temperature output
    sim = test.sims[0]
    gwe = sim.get_model(gwename)
    temp_ts = gwe.output.temperature().get_ts((0, 0, 0))
    t = temp_ts[:, 0]

    temp_analy_w = temp_z_decay(
        t, rho_w, rho_s, c_w, c_s, gamma_w, 0, n, T0
    )  # aqueous decay only
    temp_analy_s = temp_z_decay(
        t, rho_w, rho_s, c_w, c_s, 0, gamma_s, n, T0
    )  # aqueous decay only
    temp_analy_ws = temp_z_decay(
        t, rho_w, rho_s, c_w, c_s, gamma_w, gamma_s, n, T0
    )  # aqueous + solid decay

    print("temperature evaluation: " + str(temp_analy_w))

    if "aqe" in name:
        assert np.isclose(temp_ts[-1, 1], temp_analy_w[-1], atol=1e-10), msg0

    if "sld" in name:
        assert np.isclose(temp_ts[-1, 1], temp_analy_s[-1], atol=1e-10), msg1

    if "both" in name:
        assert np.isclose(temp_ts[-1, 1], temp_analy_ws[-1], atol=1e-10), msg2


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
