"""
NetCDF input tests for test_gwf_rch01: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4,
and a NetCDF input file carrying only the legacy `mesh` global attribute.
"""

import mf6_export
import pytest
from compare import Comparison
from framework import TestFramework, _copy_for_rerun
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwf_rch01 import build_models, cases, check_output


def check_legacy_mesh_attr(idx, test):
    """Backward compatibility: a UGRID NetCDF input file carrying only the
    legacy unprefixed `mesh` global attribute (as written by MF6 versions
    before the modflow_mesh rename) must still be read correctly. The file
    is MF6's own validate-mode export, with the attribute renamed.
    """
    import flopy
    import netCDF4 as nc

    ws = test.workspace / "legacy_mesh"
    _copy_for_rerun(test.workspace, ws, ["."])
    mf6_export.enable_export(ws, "layered")
    success, _ = flopy.run_model(
        test.targets["mf6"], ws / "mfsim.nam", model_ws=ws, cargs=["--mode=validate"]
    )
    assert success, "NetCDF export failed"

    (path,) = ws.glob("*.export.nc")
    with nc.Dataset(path, "r+") as ds:
        value = ds.getncattr("modflow_mesh")
        ds.delncattr("modflow_mesh")
        ds.setncattr("mesh", value)

    mf6_export.use_netcdf_input(ws)
    success, _ = flopy.run_model(test.targets["mf6"], ws / "mfsim.nam", model_ws=ws)
    assert success, "legacy mesh attr re-run failed"
    test.workspace = ws
    check_output(idx, test)


@pytest.mark.netcdf
def test_netcdf_legacy_mesh_attr(function_tmpdir, targets):
    """An exported NetCDF input file carrying only the legacy unprefixed
    `mesh` global attribute is still read correctly.
    """
    idx = 0
    test = TestFramework(
        name=cases[idx],
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_legacy_mesh_attr(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
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
    model = sims[0].gwf[0]
    add_netcdf_output(model, export)
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwf[0]
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
