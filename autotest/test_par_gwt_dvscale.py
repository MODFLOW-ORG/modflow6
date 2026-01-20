"""
This test reuses the simulation data and config in
test_gwt_dvscale.py and runs it in parallel mode.

The purpose of this test is to make sure that the dependent_variable_scaling
option works with parallel simulations.
"""

import numpy as np
import pytest
from flopy.mf6.utils import Mf6Splitter
from framework import TestFramework

cases = ["par_dvscale"]


def split_model(sim, ws):
    gwf = sim.get_model()
    nrow, ncol = gwf.dis.nrow.array, gwf.dis.ncol.array
    split_array = np.ones((nrow, ncol), dtype=int)
    split_array[:, int(ncol / 2) :] = 2

    mfsplit = Mf6Splitter(sim)
    new_sim = mfsplit.split_multi_model(split_array)
    new_sim.set_sim_path(ws)

    return new_sim


def build_models(idx, test):
    from test_gwt_dvscale import build_model as build

    ws = test.workspace
    sim = build(idx, ws, dvscale=True)
    sim = split_model(sim, ws)

    ws = test.workspace / "mf6"
    mc = build(idx, ws, dvscale=False)
    mc = split_model(mc, ws)

    return sim, mc


def check_results(test):
    from test_gwt_dvscale import check_results as base_check

    base_check(test)


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_results(t),
        compare="mf6",
        parallel=True,
        ncpus=2,
    )
    test.run()
