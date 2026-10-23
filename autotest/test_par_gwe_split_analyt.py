"""
This tests reuses the simulation data in test_gwe_split_analyt and runs it
in parallel on two cpus.  test_gwe_split_analyt runs 3 different conceptual
models documented in Carslaw and Jaeger (1947).
"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_eslcasei", "par_eslcaseii", "par_eslcaseiii"]


def build_models(idx, test):
    from test_gwe_split_analyt import build_models as build

    sim, dummy = build(idx, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwe_split_analyt import check_output as check

    check(idx, test)


@pytest.mark.parallel
@pytest.mark.developmode
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
