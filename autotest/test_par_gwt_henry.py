"""
This test reuses the simulation data and config in
test_gwt_henry_gwtgwt.py and runs it in parallel mode.
"""

import pytest
from framework import TestFramework

cases = [
    pytest.param(0, "par-henry-ups"),
    pytest.param(1, "par-henry-cen"),
    pytest.param(2, "par-henry-tvd"),
    pytest.param(3, "par-henry-utvd", marks=pytest.mark.developmode),
]


def build_models(idx, name, test):
    from test_gwt_henry_gwtgwt import build_models as build

    sim, dummy = build(idx, name, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwt_henry_gwtgwt import check_output as check

    check(idx, test)


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", cases)
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, name, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=2,
    )
    test.run()
