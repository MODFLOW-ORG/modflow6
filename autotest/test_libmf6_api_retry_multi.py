"""
Like test_libmf6_api_retry.py, but with two models, each with its own solution,
both in the same solution group. Only the first is deliberately failed. Checks
that a failed solve in the first solution causes both solutions to be retried,
even if the second succeeds.
"""

import os
from ctypes import byref, c_bool

import flopy
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libmf6_api_retry_multi"]

model_names = ["model1", "model2"]

# ATS settings: dt0, dtmin, dtmax, dtadj, dtfailadj
dt0 = 1.0
dtmin = 1.0e-5
dtmax = 1.0
dtadj = 2.0
dtfailadj = 2.0

# solver data
nouter, ninner = 100, 300
hclose, rclose = 10e-9, 1e-3


def get_model(dir):
    # tdis
    nper = 1
    perlen = 1.0
    tdis_rc = [(perlen, 1, 1.0)]

    # model spatial discretization
    nlay, nrow, ncol = 1, 1, 5

    # cell spacing
    delr = 1.0
    delc = 1.0

    # top/bot of the aquifer
    tops = [0.0, -1.0]

    # hydraulic conductivity
    k11 = 1.0

    # boundary heads
    h_left = 0.0
    h_right = 5.0

    # initial head
    h_start = 0.0

    sim = flopy.mf6.MFSimulation(
        sim_name="sim", version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # enable adaptive time stepping so failed steps are retried
    ats_filerecord = "sim.ats"
    atsperiod = [(0, dt0, dtmin, dtmax, dtadj, dtfailadj)]
    tdis.ats.initialize(
        maxats=len(atsperiod),
        perioddata=atsperiod,
        filename=ats_filerecord,
    )

    chd_data = [[(0, 0, 0), h_left], [(0, 0, ncol - 1), h_right]]
    chd_spd = {0: chd_data}

    # two models/solutions in the same solution group
    for name in model_names:
        gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
        dis = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=tops[0],
            botm=tops[1:],
        )
        ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            save_specific_discharge=True,
            save_flows=True,
            icelltype=0,
            k=k11,
        )
        chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=f"{name}.hds",
            budget_filerecord=f"{name}.cbc",
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )

        ims = flopy.mf6.ModflowIms(
            sim,
            filename=f"{name}.ims",
            print_option="SUMMARY",
            outer_dvclose=hclose,
            outer_maximum=nouter,
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=rclose,
            linear_acceleration="CG",
        )
        sim.register_ims_package(ims, [name])

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    sim = get_model(ws)

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    sim_compare = get_model(ws)

    return sim, sim_compare


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

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

    # testing API
    comp_str = mf6.get_component_name()
    version_str = mf6.get_version()
    print(f"Loaded {comp_str} with version {version_str}")

    n_solutions = mf6.get_subcomponent_count()
    if n_solutions != len(model_names):
        print(f"Expected {len(model_names)} solutions, got {n_solutions}")
        return False, open(output_file_path).readlines()

    # time loop
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()

    # count how many time-step retries the ATS retry loop performed
    n_retries = 0
    # force only the first subcomponent to fail on the first attempt
    force_failure = True

    while current_time < end_time:
        dt = mf6.get_time_step()
        try:
            mf6.prepare_time_step(dt)
        except:
            return False, open(output_file_path).readlines()

        # prepare the ATS retry loop (resets the retry counter)
        mf6._execute_function(mf6.lib.prepare_retryloop)

        attempt = 0
        while True:
            # start of a (re)try, advances IDM when retrying
            mf6._execute_function(mf6.lib.start_retry)

            # solve every subcomponent (solution) in turn, as a modflowapi
            # user driving a multi-model (e.g. loosely-coupled) simulation
            # would; only subcomponent 1 is deliberately failed, and only on
            # the very first attempt, so ATS should still trigger a retry
            for isub in range(1, n_solutions + 1):
                mf6.prepare_solve(isub)

                kiter_max = 1 if (isub == 1 and force_failure) else nouter
                kiter = 0
                has_converged = False
                while kiter < kiter_max:
                    has_converged = mf6.solve(isub)
                    kiter += 1
                    if has_converged:
                        break

                mf6.finalize_solve(isub)

            force_failure = False

            finish_retry = c_bool(False)
            mf6._execute_function(mf6.lib.finish_retry, byref(finish_retry))

            attempt += 1
            if finish_retry.value:
                break

        # everything beyond the first attempt is an ATS retry
        n_retries += attempt - 1

        try:
            mf6.finalize_time_step()
        except:
            return False, open(output_file_path).readlines()
        current_time = mf6.get_current_time()

    # the deliberately failed first subcomponent must have triggered a retry,
    # even though the second subcomponent converged right away every time
    if n_retries < 1:
        print(f"Expected at least one ATS retry, got {n_retries}")
        return False, open(output_file_path).readlines()

    # finish
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
        build=lambda t: build_models(idx, t),
        targets=targets,
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()
