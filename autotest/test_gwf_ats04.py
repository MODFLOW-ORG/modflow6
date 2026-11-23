"""
Ensure that pertim is set to the period length when a failed ATS time step
is retried with a shorter step that ends within dtmin of the end of the
period, and so ends the period.

The model is the test_gwf_ats03 model with a single 1-day period. The first
time step tries the whole period and fails to converge. It is retried with
a step of 1/1.5 days, which converges. The retried step ends 1/3 day short
of the period end, which is within dtmin (0.5 days), so the step ends the
period and the simulation. The period time and the simulation time written
for the last step should both be the period length.
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


def build(idx, test):
    sim, _ = build_models(0, test)
    sim.name = cases[idx]
    sim.tdis.perioddata.set_data([(perlen, 1, 1.0)])
    sim.tdis.ats.perioddata.set_data([(0, dt0, dtmin, dtmax, dtadj, dtfailadj)])
    return sim, None


def check_output(idx, test):
    gwfname = "gwf_ats03a"

    # the first step must have failed and been retried
    fpth = os.path.join(test.workspace, "mfsim.lst")
    with open(fpth) as f:
        lines = f.read()
    assert "will be retried" in lines, "first time step was not retried"

    fpth = os.path.join(test.workspace, gwfname + ".hds")
    hds = flopy.utils.HeadFile(fpth)
    rec = hds.recordarray
    assert rec.size == 1, f"expected a single time step, found {rec.size}"
    pertim = rec["pertim"][-1]
    totim = rec["totim"][-1]
    assert np.isclose(pertim, perlen), f"pertim is {pertim}, not {perlen}"
    assert np.isclose(totim, perlen), f"totim is {totim}, not {perlen}"


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
