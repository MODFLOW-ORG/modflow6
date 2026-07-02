"""
This tests reuses the simulation setup in test_gwt_mvt and runs it
in parallel on two cpus.  test_gwt_mvt.py runs 1 simulation that is explained
in the header section of that test.  Basically there are two gwf and two gwt
models that are adjacent to one another.  To run in parallel, one of the gwf
models and its gwt counterpart will be run on one core while the other gwf-gwt
doublet will run on another core.
"""

import pytest
from framework import TestFramework

cases = ["par_mvt"]


def build_models(idx, test):
    from test_gwt_mvt import build_models as build

    sim, dummy = build(idx, test)
    return sim, dummy


def check_output(idx, test):
    from test_gwt_mvt import check_output as check

    check(idx, test)


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        xfail=True,
        ncpus=2,
    )
    test.run()
