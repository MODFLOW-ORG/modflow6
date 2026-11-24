"""
A failed ATS time step retried with a shorter step that ends within dtmin
of the end of the period should not end the period short.

The model is the test_gwf_ats03 model with a single 1-day period. The first
time step tries the whole period and fails to converge. With dtmin = 0.5 and
dtfailadj = 1.5, only one retry is possible (1/1.5 days; 1/2.25 is below
dtmin), and it ends 1/3 day short of the period end, which is within dtmin.

Currently the retried step ends the period and the simulation: the time
actually simulated (the sum of the time step lengths) is 1/1.5 days, the
period time written for it is 1/1.5, and the simulation time is snapped to
the period length. Whatever the fix, the simulated time, the period time
and the simulation time at the end of the run should all be the period
length. Marked as an expected failure until this is fixed.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_ats03 import build_models

cases = ["gwf_ats04"]
perlen = 1.0
dt0 = perlen
dtmin = 0.5
dtmax = perlen
dtadj = 2.0
dtfailadj = 1.5

# exactly one retry is possible, and it ends within dtmin of the period end
assert perlen - perlen / dtfailadj < dtmin <= perlen / dtfailadj
assert perlen / dtfailadj**2 < dtmin


def build(idx, test):
    sim, _ = build_models(0, test)
    sim.name = cases[idx]
    sim.tdis.perioddata.set_data([(perlen, 1, 1.0)])
    sim.tdis.ats.perioddata.set_data([(0, dt0, dtmin, dtmax, dtadj, dtfailadj)])
    gwf = sim.get_model()
    gwf.npf.save_flows.set_data(True)
    gwf.oc.saverecord.set_data({0: [("HEAD", "ALL"), ("BUDGET", "ALL")]})
    return sim, None


def check_output(idx, test):
    gwfname = "gwf_ats03a"

    # preconditions: pytest.fail() rather than assert, so that they are not
    # mistaken for the expected failure (xfail only accepts AssertionError)
    fpth = os.path.join(test.workspace, "mfsim.lst")
    with open(fpth) as f:
        if "will be retried" not in f.read():
            pytest.fail("first time step was not retried")

    fpth = os.path.join(test.workspace, gwfname + ".cbc")
    rec = flopy.utils.CellBudgetFile(fpth).recordarray
    rec = rec[rec["text"] == rec["text"][0]]  # one record per time step
    simulated = rec["delt"].sum()
    pertim = rec["pertim"][-1]
    totim = rec["totim"][-1]
    msg = (
        f"simulated time {simulated}, pertim {pertim} and totim {totim} "
        f"should all be the period length {perlen}"
    )
    assert np.allclose([simulated, pertim, totim], perlen), msg


@pytest.mark.xfail(
    reason="a retried ATS step ending within dtmin of the period end "
    "ends the period short",
    raises=AssertionError,
    strict=True,
)
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
