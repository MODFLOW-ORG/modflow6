"""
This tests reuses the simulation data in test_gwf_ifmod_xt3d03.py
and runs it in parallel on 5 cpus with

cpu 0: 'ref'
cpu 1: 'tl'
cpu 2: 'bl'
cpu 3: 'tr'
cpu 4: 'br'

so we can compare the parallel coupling of the 4 models
with a serial reference case with XT3D enabled everywhere.
Particular aspect of this test is the capability to correctly
deal with the "Four Corners" issue, where for some connections
the flux has to be calculated using data from 4 subdomains.

"""

import pytest
from compare import Comparison
from framework import TestFramework

cases = ["par_xt3d03"]


def build_models(idx, test):
    from test_gwf_ifmod_xt3d03 import build_models as build

    sim, dummy = build(idx, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwf_ifmod_xt3d03 import check_output as check

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
        ncpus=5,
    )
    test.run()
