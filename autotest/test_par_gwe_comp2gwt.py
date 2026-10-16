"""
This tests reuses the simulation data in
test_gwegwe_exchng_with_comp2gwt and runs it
in parallel on two cpus
"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["gwegwe-gwtgwt"]


def build_models(idx, test):
    from test_gwegwe_exchng_with_comp2gwt import build_models as build

    sim, dummy = build(idx, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwegwe_exchng_with_comp2gwt import check_output as check

    check(idx, test)


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare", [None, pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4)]
)
def test_mf6model(idx, name, compare, function_tmpdir, targets):
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
