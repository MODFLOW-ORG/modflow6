"""
A failed ATS time step that would end the period is retried with a step
that either ends the period or leaves at least dtmin, so no time in the
period is skipped. If no such retry is possible without a step shorter
than dtmin, the simulation fails.

The model is the test_gwf_ats03 model with a single 1-day period. The first
time step tries the whole period and fails to converge (steps longer than
about 0.8 days fail for this model). Dividing by dtfailadj = 1.5 gives a
retry of 2/3 day, which would leave 1/3 day, less than dtmin.

With dtmin = 0.5, the retry is shortened to 0.5 days, leaving 0.5 for the
final step. The simulated time, the period time and the simulation time at
the end of the run are all the period length.

With dtmin = 0.6, a retry would have to be both at least 0.6 and at most
0.4 days long, so the simulation fails.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_ats03 import build_models

cases = ["gwf_ats04a", "gwf_ats04b"]
perlen = 1.0
dt0 = perlen
dtmin = [0.5, 0.6]
dtmax = perlen
dtadj = 2.0
dtfailadj = 1.5
xfail = [False, True]


def build(idx, test):
    sim, _ = build_models(0, test)
    sim.name = cases[idx]
    sim.tdis.perioddata.set_data([(perlen, 1, 1.0)])
    sim.tdis.ats.perioddata.set_data([(0, dt0, dtmin[idx], dtmax, dtadj, dtfailadj)])
    gwf = sim.get_model()
    gwf.npf.save_flows.set_data(True)
    gwf.oc.saverecord.set_data({0: [("HEAD", "ALL"), ("BUDGET", "ALL")]})
    return sim, None


def check_output(idx, test):
    gwfname = "gwf_ats03a"

    fpth = os.path.join(test.workspace, "mfsim.lst")
    with open(fpth) as f:
        lst = f.read()
    if xfail[idx]:
        assert "will be retried" not in lst, "first time step was retried"
        assert "convergence failure" in lst.lower(), "expected convergence failure"
        return
    assert "will be retried" in lst, "first time step was not retried"

    fpth = os.path.join(test.workspace, gwfname + ".cbc")
    rec = flopy.utils.CellBudgetFile(fpth).recordarray
    rec = rec[rec["text"] == rec["text"][0]]  # one record per time step
    assert np.all(rec["delt"] >= dtmin[idx]), f"step below dtmin: {rec['delt']}"
    simulated = rec["delt"].sum()
    pertim = rec["pertim"][-1]
    totim = rec["totim"][-1]
    msg = (
        f"simulated time {simulated}, pertim {pertim} and totim {totim} "
        f"should all be the period length {perlen}"
    )
    assert np.allclose([simulated, pertim, totim], perlen), msg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=xfail[idx],
    )
    test.run()
