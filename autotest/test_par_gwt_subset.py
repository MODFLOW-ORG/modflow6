"""
This test reuses the simulation data and config in
test_gwt_subset.py and runs it in parallel mode.

The purpose of this test is to make sure that the
transport model ends up on the right process to be
coupled to flow by using the HPC subpackage.
"""

import flopy
import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_gwt_subset"]


def build_models(idx, test):
    from test_gwt_subset import build_models as build

    sim, dummy = build(idx, test)

    # create the load balance for a parallel run
    partitions = [("flow1", 0), ("flow2", 1), ("transport", 1)]
    hpc = flopy.mf6.ModflowUtlhpc(sim, partitions=partitions)

    return sim, dummy


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare", [None, pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4)]
)
def test_mf6model(idx, name, compare, function_tmpdir, targets):
    from test_gwt_subset import check_output

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
        parallel=True,
        ncpus=2,
    )
    test.run()
