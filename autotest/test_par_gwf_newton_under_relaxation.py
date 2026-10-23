"""
This tests reuses the simulation data in test_gwf_newton_under_relaxation
and runs it in parallel on one and two cpus with

so we can compare the parallel coupling of two models
with a serial model.

This test also checks that Newton under_relaxation works in parallel.
"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_nr_ur01", "par_nr_ur02"]


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare", [None, pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4)]
)
def test_mf6model(idx, name, compare, function_tmpdir, targets):
    from test_gwf_newton_under_relaxation import build_models, check_output

    ncpus = 2 if idx == 1 else 1
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
        parallel=True,
        ncpus=ncpus,
    )
    test.run()
