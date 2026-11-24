"""
NetCDF input tests for test_gwe_cnd: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4.
"""

from pathlib import Path

import flopy
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwe_cnd import build_models, cases, check_output

nc = pytest.importorskip("netCDF4")

WKT1 = (
    'PROJCS["NAD83 / UTM zone 18N",'
    'GEOGCS["NAD83",'
    'DATUM["North_American_Datum_1983",'
    'SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],'
    'AUTHORITY["EPSG","6269"]],'
    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
    'AUTHORITY["EPSG","4269"]],'
    'PROJECTION["Transverse_Mercator"],'
    'PARAMETER["latitude_of_origin",0],'
    'PARAMETER["central_meridian",-75],'
    'PARAMETER["scale_factor",0.9996],'
    'PARAMETER["false_easting",500000],'
    'PARAMETER["false_northing",0],'
    'UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
    'AXIS["Easting",EAST],'
    'AXIS["Northing",NORTH],'
    'AUTHORITY["EPSG","26918"]]'
)


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
    model = sims[0].gwe[0]
    add_netcdf_output(model, export)
    flopy.mf6.ModflowUtlncf(model.dis, wkt=WKT1, filename=f"{model.name}.dis.ncf")
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwe[0]
    check_model_netcdf_output(model)
    with nc.Dataset(Path(model.model_ws) / f"{model.name}.nc") as ds:
        if export == "structured":
            proj = ds.variables["projection"]
            assert "GeoTransform" not in proj.ncattrs()
            assert "spatial_ref" not in proj.ncattrs()
            assert proj.getncattr("wkt") == WKT1
            assert proj.getncattr("grid_mapping_name") == "transverse_mercator"
            assert proj.getncattr("longitude_of_central_meridian") == -75.0
            assert proj.getncattr("latitude_of_projection_origin") == 0.0
            assert proj.getncattr("scale_factor_at_central_meridian") == 0.9996
            assert proj.getncattr("false_easting") == 500000.0
            assert proj.getncattr("false_northing") == 0.0
            assert proj.getncattr("semi_major_axis") == 6378137.0
            assert proj.getncattr("inverse_flattening") == 298.257222101
            z = ds.variables["z"]
            coords = ds.variables["temperature"].getncattr("coordinates").split()
        else:
            # a single-layer model
            z = ds.variables["z_l1"]
            assert z.getncattr("layer") == 1
            coords = ds.variables["temperature_l1"].getncattr("coordinates").split()
        assert z.getncattr("standard_name") == "altitude"
        assert z.getncattr("positive") == "up"
        assert z.name in coords


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
