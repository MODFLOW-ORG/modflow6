"""
Parallel (PETSc) version of test_gwf_ims_period_output.

The same two coupled GWF models are split over two cpus. Each rank writes its
own mfsim.p*.lst and IMS csv files, and the shared check confirms every rank
recorded solver output consistent with the per-period linear settings.
"""

import pytest
from framework import TestFramework

cases = ["par_ims_period_out"]


def build_models(idx, test):
    from test_gwf_ims_period_output import build_models as build

    return build(idx, test)


def check_output(idx, test):
    from test_gwf_ims_period_output import check_output as check

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
        overwrite=False,  # keep the injected LINEAR_PERIODDATA + period file
    )
    test.run()
