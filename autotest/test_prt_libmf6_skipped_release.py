"""
Test that PRP drops unused release times: negative times,
and times after the end of the simulation unless
EXTEND_TRACKING is on. Run through the API and check the
committed particle store after each time step, and check
the warnings in the listing.

The flow system is the one on the FloPy readme. Particles
are released from the top left cell.

Cases:
  late:    a time after the end, dropped
  lateext: a time after the end, kept with EXTEND_TRACKING
  alllate: all times after the end, dropped, no default
           release at t=0
  neg:     a negative time, dropped
"""

from pathlib import Path

import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from prt_test_utils import FlopyReadmeCase, get_model_name

nstp = 3  # simulation ends at t=1
late_warning = "Skipping release time falling after the end"
neg_warning = "Skipping negative release time"

# released: release times in the store after each step
cases = {
    "prtlate": {
        "releasetimes": [0.0, 5.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {late_warning: 1, neg_warning: 0},
    },
    "prtlateext": {
        "releasetimes": [0.0, 5.0],
        "extend": True,
        "released": [[0.0], [0.0], [0.0, 5.0]],
        "warnings": {late_warning: 0, neg_warning: 0},
    },
    "prtalllate": {
        "releasetimes": [5.0, 6.0],
        "extend": False,
        "released": [[], [], []],
        "warnings": {late_warning: 2, neg_warning: 0},
    },
    "prtneg": {
        "releasetimes": [-1.0, 0.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {late_warning: 0, neg_warning: 1},
    },
}
names = list(cases)


def build_models(idx, test):
    case = cases[test.name]
    return FlopyReadmeCase.get_prt_sim(
        test.name,
        test.workspace,
        test.targets["libmf6"],
        nstp=nstp,
        releasetimes=case["releasetimes"],
        extend_tracking=case["extend"],
    )


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = names[idx]
    case = cases[name]
    model_ws = Path(model_ws or ".")
    prt_name = get_model_name(name, "prt").upper()
    nrpt = len(FlopyReadmeCase.releasepts_prt)

    mf6 = ModflowApi(exe, working_directory=model_ws)
    mf6.initialize()

    def get(var):
        return mf6.get_value(mf6.get_var_address(var, prt_name, "PRP1"))

    try:
        kstp = 0
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.update()
            kstp += 1
            released = case["released"][kstp - 1]
            irpt = get("PLIRPT")
            assert irpt.size == nrpt * len(released), (
                f"step {kstp}: expected {nrpt * len(released)} particles, "
                f"got {irpt.size}"
            )
            assert np.all((irpt >= 1) & (irpt <= nrpt))
            assert np.allclose(np.unique(get("PLTRELEASE")), released)
            assert get("NPARTICLES")[0] == irpt.size
        assert kstp == nstp
    finally:
        mf6.finalize()

    return True, open(model_ws / "mfsim.stdout").readlines()


def check_output(idx, test):
    stdout = (test.workspace / "mfsim.stdout").read_text()
    for warning, count in cases[test.name]["warnings"].items():
        assert stdout.count(warning) == count, warning


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(names))
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
