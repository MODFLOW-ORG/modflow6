"""
Test the bmi which is used update to set the river stages to
the same values as they are in the non-bmi simulation.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_riv01"]

# temporal discretization
nper = 10
tdis_rc = []
for i in range(nper):
    tdis_rc.append((1.0, 1, 1))

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 10

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
h1 = 11.0
h2 = 11.0

# build chd stress period data
chd_spd = {0: [[(0, 0, 0), h1], [(0, 0, ncol - 1), h2]]}

strt = np.linspace(h1, h2, num=ncol)


# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97

# uniform river stage
riv_stage = 15.0
riv_stage2 = 20.0
riv_bot = 12.0
riv_cond = 35.0
riv_packname = "MYRIV"


def get_model(ws, name, riv_spd):
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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf, save_flows=True, iconvert=1, ss=0.0, sy=0.2, transient={0: True}
    )

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # riv package
    riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd, pname=riv_packname)

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

    # create river data
    rd = [[(0, 0, icol), riv_stage, riv_cond, riv_bot] for icol in range(1, ncol - 1)]
    rd2 = [[(0, 0, icol), riv_stage2, riv_cond, riv_bot] for icol in range(1, ncol - 1)]
    sim = get_model(ws, name, riv_spd={0: rd, 5: rd2})

    # build comparison model with zeroed values
    ws = os.path.join(test.workspace, "libmf6")
    rd_bmi = [[(0, 0, icol), 999.0, 999.0, 0.0] for icol in range(1, ncol - 1)]
    mc = get_model(ws, name, riv_spd={0: rd_bmi})

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

    # get copy of (multi-dim) array with river parameters
    stage_tag = mf6.get_var_address("STAGE", name, riv_packname)
    cond_tag = mf6.get_var_address("COND", name, riv_packname)
    rbot_tag = mf6.get_var_address("RBOT", name, riv_packname)
    new_stage = mf6.get_value(stage_tag)
    new_cond = mf6.get_value(cond_tag)
    new_rbot = mf6.get_value(rbot_tag)

    # model time loop
    idx = 0
    while current_time < end_time:
        # get dt
        dt = mf6.get_time_step()

        # prepare... and reads the RIV data from file!
        mf6.prepare_time_step(dt)

        # set the RIV data through the BMI
        # change cond and rbot data
        new_cond[:] = [riv_cond]
        new_rbot[:] = [riv_bot]
        mf6.set_value(cond_tag, new_cond)
        mf6.set_value(rbot_tag, new_rbot)
        # change stage data
        if current_time < 5:
            new_stage[:] = [riv_stage]
            mf6.set_value(stage_tag, new_stage)
        else:
            new_stage[:] = [riv_stage2]
            mf6.set_value(stage_tag, new_stage)

        kiter = 0
        mf6.prepare_solve()

        while kiter < nouter:
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
