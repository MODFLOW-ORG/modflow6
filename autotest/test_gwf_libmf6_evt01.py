"""
Test bmi with a head-based pumping rate equivalent to
the evapotranspiration package in a non-bmi simulation.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_evt01"]

# et variables
et_max = 0.1
et_depth = 5.0

# temporal discretization
nper = 100
tdis_rc = []
for i in range(nper):
    tdis_rc.append((1.0, 1, 1))

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 1

# cell spacing
delr = 10.0
delc = 10.0
area = delr * delc

# top of the aquifer
top = 10.0

# bottom of the aquifer
botm = 0.0

# hydraulic conductivity
hk = 50.0

# starting head
strt = 15.0

# solver data
nouter, ninner = 100, 100
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def get_model(ws, name, bmi=False):
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
        under_relaxation="SIMPLE",
        under_relaxation_gamma=0.98,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )

    # create gwf model
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        newtonoptions=newtonoptions,
        modelname=name,
        print_input=True,
        save_flows=True,
    )

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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf, save_flows=True, iconvert=1, ss=1e-5, sy=0.2, transient={0: True}
    )

    # evapotranspiration
    if not bmi:
        evt = flopy.mf6.ModflowGwfevta(gwf, surface=top, rate=et_max, depth=et_depth)
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(0, 0, 0), 0.0]])

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
    sim = get_model(ws, name)

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    mc = get_model(ws, name, bmi=True)

    return sim, mc


def head2et_wellrate(h):
    if h > top:
        q = -et_max
    elif h < top - et_depth:
        q = 0.0
    else:
        f = (h - (top - et_depth)) / et_depth
        q = -f * et_max
    return q * area


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

    # get pointer to simulated heads
    head_tag = mf6.get_var_address("X", "LIBGWF_EVT01")
    head = mf6.get_value_ptr(head_tag)

    # maximum outer iterations
    mxit_tag = mf6.get_var_address("MXITER", "SLN_1")
    max_iter = mf6.get_value(mxit_tag)

    # get copy of well data
    well_tag = mf6.get_var_address("Q", name, "WEL_0")
    well = mf6.get_value(well_tag)

    # check NPF type
    package_type_tag = mf6.get_var_address("PACKAGE_TYPE", name, "NPF")
    package_type = mf6.get_value(package_type_tag)[0]
    assert package_type == "NPF", f"{package_type} /= NPF"

    # check wel type
    package_type_tag = mf6.get_var_address("PACKAGE_TYPE", name, "WEL_0")
    package_type = mf6.get_value(package_type_tag)[0]
    assert package_type == "WEL", f"{package_type} /= WEL"

    twell = np.zeros(ncol, dtype=np.float64)

    # model time loop
    idx = 0
    while current_time < end_time:
        # get dt and prepare for non-linear iterations
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)

        # convergence loop
        kiter = 0
        mf6.prepare_solve()

        while kiter < max_iter:
            # update well rate
            twell[:] = head2et_wellrate(head[0])
            well[:] = twell[:]
            mf6.set_value(well_tag, well)

            # solve with updated well rate
            has_converged = mf6.solve()
            kiter += 1

            if has_converged:
                msg = f"Component {1}" + f" converged in {kiter}" + " outer iterations"
                print(msg)
                break

        if not has_converged:
            print("model did not converge")
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
