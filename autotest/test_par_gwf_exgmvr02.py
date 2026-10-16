"""
This tests reuses the simulation data in test_gwf_exgmvr02.py
and runs it in parallel on two cpus with

cpu 1: 'left'
cpu 2: 'right'

so we can compare the parallel coupling of 'left' + 'right'
with a serial 'single'
"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_exgmvr02"]


def build_models(idx, test):
    from test_gwf_exgmvr02 import build_models as build

    sim, sim_ref = build(idx, test)
    return sim, sim_ref


def check_output(idx, test):
    from test_gwf_exgmvr02 import check_output as check

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
        parallel=True,
        ncpus=(2, 1),
        compare=compare,
    )
    test.run()
