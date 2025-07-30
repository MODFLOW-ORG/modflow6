"""
This test reuses the simulation data and config in
test_gwt_dvscale.py and runs it in parallel mode.

The purpose of this test is to make sure that the dependent_variable_scaling 
option works with paralle simulations.
"""

import numpy as np
import pytest
from flopy.mf6.utils import Mf6Splitter
from framework import TestFramework

cases = ["par_dvscale"]
options = [True]


def build_models(idx, test):
    from test_gwt_dvscale import build_models as build

    sim, dummy = build(idx, test)

    gwf = sim.get_model()
    nrow, ncol = gwf.dis.nrow.array, gwf.dis.ncol.array
    split_array = np.ones((nrow, ncol), dtype=int)
    split_array[:, int(ncol / 2) :] = 2

    mfsplit = Mf6Splitter(sim)
    new_sim = mfsplit.split_multi_model(split_array)
    new_sim.set_sim_path(test.workspace)

    return new_sim, dummy


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
        compare=None,
        parallel=True,
        ncpus=2,
    )
    test.run()
