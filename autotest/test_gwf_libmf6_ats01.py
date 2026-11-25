"""
Run the test_gwf_ats03 model through the API with update_until, and check
that it reaches the requested end time after the first time step is retried.
"""

import os
from types import SimpleNamespace

import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from test_gwf_ats03 import build_models

cases = ["libgwf_ats01"]


def build(test):
    sim, _ = build_models(0, test)
    # comparison sim in a "libmf6" subdir: the framework runs this one via the API
    ws = SimpleNamespace(workspace=test.workspace / "libmf6")
    sim_api, _ = build_models(0, ws)
    return sim, sim_api


def api_func(exe, model_ws):
    from modflowapi import ModflowApi

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

    end_time = mf6.get_end_time()
    mf6.update_until(end_time)
    current_time = mf6.get_current_time()

    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    msg = f"update_until stopped at time {current_time}, not {end_time}"
    assert current_time == pytest.approx(end_time), msg

    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build,
        targets=targets,
        api_func=api_func,
    )
    test.run()
