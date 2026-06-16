"""
Test hydraulic conductivity updates through the API, both via set_value
and via pointer mutation followed by an explicit on_value_changed signal.
"""

import os
from ctypes import c_char_p

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["api_k_setvalue", "api_k_ptr"]

# temporal discretization
nper = 1
tdis_rc = [(1.0, 1, 1.0)]

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 11

# cell spacing
delr = 50.0
delc = 1.0

# top of the aquifer
top = 10.0

# bottom of the aquifer
botm = 0.0

# hydraulic conductivity
initial_hk = 0.5
updated_hk = 25.0

# boundary heads
h1 = 10.0
h2 = 5.0

# recharge
recharge = 0.02

# build chd stress period data
chd_spd = {0: [[(0, 0, 0), h1], [(0, 0, ncol - 1), h2]]}

# solver data
nouter, ninner = 100, 300
hclose, rclose = 1e-9, 1e-3


def get_model(ws, name, hk):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="all",
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    flopy.mf6.ModflowGwfic(gwf, strt=np.linspace(h1, h2, num=ncol))
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=hk,
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name, hk=updated_hk)

    ws = os.path.join(test.workspace, "libmf6")
    sim_compare = get_model(ws, name, hk=initial_hk)

    return sim, sim_compare


def on_value_changed(mf6, var_address):
    mf6._execute_function(mf6.lib.on_value_changed, c_char_p(var_address.encode()))


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx].upper()
    if model_ws is None:
        model_ws = "."

    output_file_path = os.path.join(model_ws, "mfsim.stdout")

    try:
        mf6 = ModflowApi(exe, working_directory=model_ws)
    except Exception as e:
        print("Failed to load " + str(exe))
        print("with message: " + str(e))
        return False, open(output_file_path).readlines()

    try:
        mf6.initialize()
    except:
        return False, open(output_file_path).readlines()

    k11_tag = mf6.get_var_address("K11", name, "NPF")
    if cases[idx].endswith("setvalue"):
        new_k11 = mf6.get_value(k11_tag)
        new_k11.fill(updated_hk)
        mf6.set_value(k11_tag, new_k11)
    else:
        k11 = mf6.get_value_ptr(k11_tag)
        k11.fill(updated_hk)
        on_value_changed(mf6, k11_tag)

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()

    while current_time < end_time:
        try:
            mf6.update()
        except:
            return False, open(output_file_path).readlines()
        current_time = mf6.get_current_time()

    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()