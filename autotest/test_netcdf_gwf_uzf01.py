"""
NetCDF input version of test_gwf_uzf01: the simulations re-run with their
package arrays read from NetCDF, as exported by MODFLOW 6 and as written
by flopy4, and the run is compared to the ASCII run and
re-checked. Also verifies the list-based (GHB) and array-based (GHBG)
representations agree with each other, not just with their own ASCII
baseline, since the generic comparison only ever compares a
representation against its own rewrite.
"""

import flopy
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwf_uzf01 import build_models, cases, check_output, check_outputs


@pytest.mark.netcdf
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
@pytest.mark.parametrize(
    "compare",
    [
        Comparison.MF6_EXPORT_STRUCTURED,
        Comparison.MF6_EXPORT_LAYERED,
        Comparison.FP4_STRUCTURED,
        Comparison.FP4_LAYERED,
    ],
)
def test_netcdf_input(idx, name, function_tmpdir, targets, compare):
    """The base test re-run with NetCDF input matches the
    ASCII run, and the list-based and array-based representations agree
    with each other.
    """

    def check_rerun(test):
        # test.workspace is already retargeted to the comparison dir
        list_ws = test.workspace
        array_ws = list_ws / "mf6"
        check_output(list_ws, name)
        check_output(array_ws, name)

        # the generic comparison only ever checks a representation
        # against its own ASCII baseline; list-based (GHB) and array-based
        # (GHBG) must also agree with each other
        success = flopy.utils.compare.compare_heads(
            None,
            None,
            text="head",
            outfile=str(array_ws / f"{name}.hds.cmp.out"),
            files1=str(list_ws / f"{name}.hds"),
            files2=str(array_ws / f"{name}.hds"),
            difftol=True,
        )
        assert success, "list vs array head comparison failed"

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_outputs(idx, t),
        check_rerun=check_rerun,
        targets=targets,
        compare=compare,
    )
    test.run()


def build_output_models(idx, test, export):
    """The base simulation, writing its NetCDF output file."""
    sims = build_models(idx, test)
    model = sims[1].gwf[0]
    add_netcdf_output(model, export)
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_outputs(idx, test)
    sims = test.sims
    model = sims[1].gwf[0]
    check_model_netcdf_output(model)


@pytest.mark.netcdf
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
