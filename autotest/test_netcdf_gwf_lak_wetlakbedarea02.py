"""
NetCDF input version of test_gwf_lak_wetlakbedarea02: flopy4 re-writes the
simulation with NetCDF package arrays, and the run is compared to the ASCII run
and re-checked.
"""

import pytest
from compare import Comparison
from framework import TestFramework
from test_gwf_lak_wetlakbedarea02 import cases


@pytest.mark.netcdf
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
@pytest.mark.parametrize("compare", [Comparison.FP4_STRUCTURED, Comparison.FP4_LAYERED])
def test_mf6model_fp4(idx, name, function_tmpdir, targets, compare):
    """The base test re-written by flopy4 with NetCDF input matches the ASCII run."""
    from test_gwf_lak_wetlakbedarea02 import build_models as build
    from test_gwf_lak_wetlakbedarea02 import check_output as check

    def check_both(test):
        check(idx, test)
        test.workspace = test.workspace / compare.value
        check(idx, test)

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build(idx, t),
        check=check_both,
        targets=targets,
        compare=compare,
    )
    test.run()
