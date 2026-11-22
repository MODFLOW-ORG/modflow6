"""
Test the PRP release time coincidence tolerance.

The default RELEASE_TIME_TOLERANCE scales with the simulation
length T (the sum of all stress period lengths), so whether two
nearby release times are merged depends on their spacing relative
to T, not on the model's time units. Separately, release times
indistinguishable at floating-point precision are always merged,
even if the specified tolerance is below the ulp spacing.

Each case releases particles from one point at two explicit times
(or one explicit time and one period-block time), and checks how
many distinct releases occurred.

Cases:
    - floor: T = 2e10, period-block FIRST release in period 2 at
      t = 1e10, and an explicit release 10 ulps earlier. Tolerance
      is set far below the ulp spacing. Expect a single release.
    - sclsml: T = 1, explicit releases 1e-6 * T apart. Expect two.
      (The old fixed default, ~2.2e-5, merged these.)
    - sclbig: T = 1e6, explicit releases 1e-6 * T apart. Expect two.
    - short: T = 1, explicit releases 1e-7 apart. Expect two.
      (The old fixed default merged these.)
    - long: T = 1e6, explicit releases 1e-4 apart. Expect one.
      (The old fixed default kept these separate.)
"""

import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prttol"

# case name: (period lengths, explicit release times, period block, rttol, nrel)
FLOOR_T = 1e10
EPS = np.finfo(float).eps
DEFAULT_RTTOL = 100 * EPS * 1e5
CASES = {
    "floor": (
        [FLOOR_T, FLOOR_T],
        [FLOOR_T - 10 * np.spacing(FLOOR_T)],
        {0: [], 1: [("FIRST",)]},
        1e-12,
        1,
    ),
    "sclsml": ([1.0], [0.5, 0.5 + 1e-6], None, None, 2),
    "sclbig": ([1e6], [0.5e6, 0.5e6 + 1.0], None, None, 2),
    "short": ([1.0], [0.5, 0.5 + 1e-7], None, None, 2),
    "long": ([1e6], [0.5e6, 0.5e6 + 1e-4], None, None, 1),
}
cases = [f"{simname}{k}" for k in CASES]


def get_case(name):
    return CASES[name.removeprefix(simname)]


def build_models(test):
    perlens, releasetimes, perioddata, rttol, _ = get_case(test.name)
    sim = FlopyReadmeCase.get_prt_sim(
        test.name,
        test.workspace,
        test.targets["mf6"],
        nstp=1,
        releasetimes=[float(t) for t in releasetimes],
        perioddata=perioddata,
        perlens=perlens,
        releasepts=FlopyReadmeCase.releasepts_prt[:1],
        release_time_tolerance=rttol,
    )
    # write release times with enough digits to round-trip
    sim.simulation_data.float_precision = 16
    sim.simulation_data.float_characters = 24
    return sim


def check_output(test):
    name = test.name
    perlens, _, _, rttol, nrel = get_case(name)
    prt_ws = test.workspace
    prt_name = get_model_name(name, "prt")

    # check the number of distinct releases
    pls = pd.read_csv(prt_ws / f"{prt_name}.trk.csv")
    releases = pls[pls["ireason"] == 0]
    assert len(releases) == nrel, releases
    assert releases["trelease"].nunique() == nrel

    # check the effective tolerance is reported in the listing file
    lines = [l.strip() for l in open(prt_ws / f"{prt_name}.lst")]
    tol_lines = [l for l in lines if l.startswith("RELEASE TIME TOLERANCE =")]
    assert len(tol_lines) == 1
    tol_value = float(tol_lines[0].split("=")[1].split("(")[0])
    if rttol is None:
        assert "(DEFAULT)" in tol_lines[0]
        assert np.isclose(tol_value, DEFAULT_RTTOL * sum(perlens), rtol=1e-6)
    else:
        assert "(SPECIFIED)" in tol_lines[0]
        assert np.isclose(tol_value, rttol, rtol=1e-6)


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t),
        check=lambda t: check_output(t),
        targets=targets,
        compare=None,
    )
    test.run()
