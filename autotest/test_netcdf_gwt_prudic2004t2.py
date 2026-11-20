"""
NetCDF input tests for test_gwt_prudic2004t2: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4.
"""

import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwt_prudic2004t2 import build_models, cases, check_output


@pytest.mark.netcdf
@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare",
    [
        Comparison.MF6_EXPORT_STRUCTURED,
        Comparison.MF6_EXPORT_LAYERED,
        Comparison.FP4_STRUCTURED,
        Comparison.FP4_LAYERED,
    ],
)
def test_netcdf_input(idx, name, compare, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
    )
    test.run()


def build_output_models(idx, test, export):
    """The base simulation, writing its NetCDF output file."""
    sims = build_models(idx, test)
    sims[0].tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    model = sims[0].gwt[0]
    add_netcdf_output(model, export)
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwt[0]
    check_model_netcdf_output(model)


@pytest.mark.netcdf
@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
def test_netcdf_output(idx, name, export, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_output_models(idx, t, export),
        check=lambda t: check_output_file(idx, t, export),
        compare=None,
    )
    test.run()
