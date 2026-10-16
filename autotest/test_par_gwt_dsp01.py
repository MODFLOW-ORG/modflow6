"""
This test reuses the simulation data and config in
test_gwt_dsp01_gwtgwt.py and runs it in parallel mode.
"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_dsp01_gwtgwt"]


def build_models(idx, test):
    from test_gwt_dsp01_gwtgwt import build_models as build

    sim, dummy = build(idx, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwt_dsp01_gwtgwt import check_output as check

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
