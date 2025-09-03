"""
This test runs the simulations in test_prt_budget.py, but uses the api
to run the PRT test.
"""

from pathlib import Path

import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from test_prt_budget import build_gwf_sim, build_gwt_sim, build_prt_sim, check_output

simname = "prt_libmf6"
cases = [simname]


def build_models(idx, test):
    gwf_sim = build_gwf_sim(test.workspace / "gwf", test.targets["libmf6"])
    gwt_sim = build_gwt_sim(
        test.workspace / "gwf", test.workspace / "gwt", test.targets["libmf6"]
    )
    prt_sim = build_prt_sim(
        test.workspace / "gwf",
        test.workspace / "prt",
        test.targets["libmf6"],
    )

    return gwf_sim, gwt_sim, prt_sim


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx].upper()
    if model_ws is None:
        model_ws = Path(".")

    output_file_path = model_ws / "mfsim.stdout"

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

    # model time loop
    idx = 0
    while current_time < end_time:
        # get dt and prepare for non-linear iterations
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)

        # convergence loop
        kiter = 0
        mf6.prepare_solve()
        has_converged = mf6.solve()
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
@pytest.mark.skip(reason="storage term disagreement")
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()
