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
  default: no times or period blocks, default release at t=0
  blklate: only a period block numbered above NPER, no
           default release
  bothlate: all times after the end and a period block
           numbered above NPER, no default release

Cases that release nothing must warn once, naming why.
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
none_warning = "No particles will be released"
blk_cause = "all PERIOD blocks are numbered above NPER"
time_cause = "all release times are negative or fall after the end"
none_warnings = {none_warning: 0, blk_cause: 0, time_cause: 0}

# released: release times in the store after each step
cases = {
    "prtlate": {
        "releasetimes": [0.0, 5.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {late_warning: 1, neg_warning: 0, **none_warnings},
    },
    "prtlateext": {
        "releasetimes": [0.0, 5.0],
        "extend": True,
        "released": [[0.0], [0.0], [0.0, 5.0]],
        "warnings": {late_warning: 0, neg_warning: 0, **none_warnings},
    },
    "prtalllate": {
        "releasetimes": [5.0, 6.0],
        "extend": False,
        "released": [[], [], []],
        "warnings": {
            late_warning: 2,
            neg_warning: 0,
            none_warning: 1,
            blk_cause: 0,
            time_cause: 1,
        },
    },
    "prtneg": {
        "releasetimes": [-1.0, 0.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {late_warning: 0, neg_warning: 1, **none_warnings},
    },
    "prtdefault": {
        "releasetimes": [],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {late_warning: 0, neg_warning: 0, **none_warnings},
    },
    # nper is 1, so period block 2 is never loaded
    "prtblklate": {
        "releasetimes": [],
        "extend": False,
        "perioddata": {1: [("FIRST",)]},
        "released": [[], [], []],
        "warnings": {
            late_warning: 0,
            neg_warning: 0,
            none_warning: 1,
            blk_cause: 1,
            time_cause: 0,
        },
    },
    "prtbothlate": {
        "releasetimes": [5.0, 6.0],
        "extend": False,
        "perioddata": {1: [("FIRST",)]},
        "released": [[], [], []],
        "warnings": {
            late_warning: 2,
            neg_warning: 0,
            none_warning: 1,
            blk_cause: 1,
            time_cause: 1,
        },
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
        perioddata=case.get("perioddata"),
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
    # warnings are wrapped, so join lines
    stdout = " ".join((test.workspace / "mfsim.stdout").read_text().split())
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
