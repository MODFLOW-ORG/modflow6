"""
Reuses the simulation data in test_gwe_idomain_subset_gwegwe and runs it in parallel
on two processes, so that the GWE-GWE exchange spans the process boundary while
the downgradient transport model uses a subset of the flow model cells.
"""

import pytest
from framework import TestFramework
from test_gwe_idomain_subset_gwegwe import build_models as build
from test_gwe_idomain_subset_gwegwe import cases as serial_cases
from test_gwe_idomain_subset_gwegwe import check_output as check

cases = ["par_" + name for name in serial_cases]


def build_models(idx, test):
    return build(idx, test)


def check_output(idx, test):
    check(idx, test)


@pytest.mark.parallel
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
        ncpus=2,
    )
    test.run()
