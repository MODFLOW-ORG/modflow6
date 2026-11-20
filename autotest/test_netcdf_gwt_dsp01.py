"""
NetCDF tests for test_gwt_dsp01: the NetCDF model output file, and
NetCDF package input, as exported by MODFLOW 6 and as written by flopy4.
"""

import flopy
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwt_dsp01 import build_models as build_base
from test_gwt_dsp01 import cases
from test_gwt_dsp01 import check_output as check_base

nc = pytest.importorskip("netCDF4")

# WKT1 for EPSG:26918 (NAD83 / UTM zone 18N)
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


def add_ncf(gwt, export):
    """NetCDF storage options and a CRS for the transport model's files."""
    if export == "ugrid":
        chunks = {"chunk_face": 5}
    else:
        chunks = {"chunk_z": 1, "chunk_y": 1, "chunk_x": 20, "wkt": WKT1}
    flopy.mf6.ModflowUtlncf(
        gwt.dis,
        deflate=3,
        shuffle=False,
        chunk_time=1,
        filename=f"{gwt.name}.dis.ncf",
        **chunks,
    )


def build_output_models(idx, test, export):
    """The base simulation, writing concentration every time step to a
    NetCDF output file with compression, chunking and a CRS."""
    sim, _ = build_base(idx, test)
    sim.tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    gwt = sim.gwt[0]
    gwtname = gwt.name
    add_netcdf_output(gwt, export)
    add_ncf(gwt, export)
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    )
    return sim, None


def build_input_models(idx, test, compare):
    """The base simulation; for MF6's own export, with the transport
    model's arrays exported compressed and chunked."""
    sim, _ = build_base(idx, test)
    if compare in EXPORT_FORMAT:
        add_ncf(sim.gwt[0], EXPORT_FORMAT[compare])
    return sim, None


def check_output_file(idx, test, export):
    """The NetCDF output file has the requested storage options and CRS,
    and holds the concentrations of the binary output file."""
    check_base(idx, test)
    gwtname = f"gwt_{cases[idx]}"
    fpth = test.workspace / f"{gwtname}.nc"

    with nc.Dataset(fpth) as ds:
        assert ds.data_model == "NETCDF4"
        var = ds.variables[
            "concentration" if export == "structured" else "concentration_l1"
        ]
        assert var.chunking() == ([1, 1, 1, 20] if export == "structured" else [1, 5])
        filters = var.filters()
        assert not filters["shuffle"]
        assert filters["complevel"] == 3
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

    check_model_netcdf_output(test.sims[0].gwt[0])


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


def check_exported_input(idx, test):
    """The re-run with MF6's exported input, read from a compressed file
    (the framework checks the exported arrays hold the flopy arrays)."""
    check_base(idx, test)
    with nc.Dataset(test.workspace / f"gwt_{cases[idx]}.export.nc") as ds:
        inputs = [v for v in ds.variables.values() if "modflow_input" in v.ncattrs()]
        assert inputs and all(v.filters()["complevel"] == 3 for v in inputs)


EXPORT_FORMAT = {
    Comparison.MF6_EXPORT_STRUCTURED: "structured",
    Comparison.MF6_EXPORT_LAYERED: "ugrid",
}


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
        build=lambda t: build_input_models(idx, t, compare),
        check=lambda t: check_base(idx, t),
        check_rerun=(
            (lambda t: check_exported_input(idx, t))
            if compare in EXPORT_FORMAT
            else None
        ),
        compare=compare,
    )
    test.run()
