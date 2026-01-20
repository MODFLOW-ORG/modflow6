"""
Test to make sure that recharge is passed to the highest active layer and
verify that recharge is in the highest active layer by looking at the
individual budget terms.  For this test, there are two layers and five
columns.  The top layer is dry except for the middle cell.  Recharge is
applied to the top layer.  In the test a, IRCH is not specified.  In test b
IRCH is specified as 1, and in test c IRCH is specified as [2, 2, 1, 2, 2]
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_rch01"]

# recharge package name
rch_pname = "RCH-1"

# average recharge rate
avg_rch = 0.001

# calculate recharge rates
dx = 1 / 20
rad = np.arange(0, 1 + dx, dx) * 2.0 * np.pi
f = np.sin(rad)
rch_rates = avg_rch + f * avg_rch

# temporal discretization
nper = rch_rates.shape[0]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((1.0, 1, 1))

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 100

# cell spacing
delr = 50.0
delc = 1.0
area = delr * delc

# top of the aquifer
top = 25.0

# bottom of the aquifer
botm = 0.0

# hydraulic conductivity
hk = 50.0

# boundary heads
h1 = 20.0
h2 = 11.0

# build chd stress period data
chd_spd = {0: [[(0, 0, 0), h1], [(0, 0, ncol - 1), h2]]}

# build recharge spd
rch_spd = {}
for n in range(nper):
    rch_spd[n] = rch_rates[n]

# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def get_model(ws, name, rech):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="all",
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=top)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk)

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # recharge file
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rech, pname=rch_pname)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name, rech=rch_spd)

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    mc = get_model(ws, name, rech=0.0)

    return sim, mc


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

    # initialize the model
    try:
        mf6.initialize()
    except:
        return False, open(output_file_path).readlines()

    # time loop
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()

    # maximum outer iterations
    mxit_tag = mf6.get_var_address("MXITER", "SLN_1")
    max_iter = mf6.get_value(mxit_tag)

    # get copy of recharge array
    rch_tag = mf6.get_var_address("RECHARGE", name, rch_pname)
    new_recharge = mf6.get_value(rch_tag)

    # model time loop
    idx = 0
    while current_time < end_time:
        # get dt and prepare for non-linear iterations
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)

        # convergence loop
        kiter = 0
        mf6.prepare_solve()

        # update recharge
        new_recharge[:] = rch_spd[idx]
        mf6.set_value(rch_tag, new_recharge)

        while kiter < max_iter:
            has_converged = mf6.solve()
            kiter += 1

            if has_converged:
                msg = f"Component {1}" + f" converged in {kiter}" + " outer iterations"
                print(msg)
                break

        if not has_converged:
            return False, open(output_file_path).readlines()

        # finalize time step
        mf6.finalize_solve()

        # finalize time step and update time
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

        # increment counter
        idx += 1

    # cleanup
    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    # cleanup and return
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
