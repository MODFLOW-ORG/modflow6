"""
NetCDF input tests for test_gwf_sto01: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4.
"""

from pathlib import Path

import flopy
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwf_sto01 import build_models, cases, check_output, htol

nc = pytest.importorskip("netCDF4")

wkt = (
    'PROJCS["NAD83 / UTM zone 18N", '
    'GEOGCS["NAD83", '
    'DATUM["North_American_Datum_1983", '
    'SPHEROID["GRS 1980",6378137,298.257222101], '
    "TOWGS84[0,0,0,0,0,0,0]], "
    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]], '
    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]], '
    'AUTHORITY["EPSG","4269"]], '
    'PROJECTION["Transverse_Mercator"], '
    'PARAMETER["latitude_of_origin",0], '
    'PARAMETER["central_meridian",-75], '
    'PARAMETER["scale_factor",0.9996], '
    'PARAMETER["false_easting",500000], '
    'PARAMETER["false_northing",0], '
    'UNIT["metre",1,AUTHORITY["EPSG","9001"]], '
    'AXIS["Easting",EAST], '
    'AXIS["Northing",NORTH], '
    'AUTHORITY["EPSG","26918"]]'
)


def build_input_models(idx, test, compare):
    """The base simulation; for MF6's own export, with its arrays exported
    compressed and chunked."""
    sims = build_models(idx, test)
    if compare in (Comparison.MF6_EXPORT_STRUCTURED, Comparison.MF6_EXPORT_LAYERED):
        model = sims[0].gwf[0]
        if compare == Comparison.MF6_EXPORT_LAYERED:
            chunks = {"chunk_face": 10}
        else:
            chunks = {"chunk_z": 1, "chunk_y": 5, "chunk_x": 5}
        flopy.mf6.ModflowUtlncf(
            model.dis,
            deflate=5,
            shuffle=True,
            chunk_time=1,
            wkt=wkt,
            filename=f"{model.name}.dis.ncf",
            **chunks,
        )
    return sims


def check_compressed_input(idx, test, complevel):
    """The base checks, on a re-run reading compressed input."""
    check_output(idx, test)
    for path in test.workspace.glob("*.export.nc"):
        with nc.Dataset(path) as ds:
            inputs = [
                v for v in ds.variables.values() if "modflow_input" in v.ncattrs()
            ]
            assert inputs and all(v.filters()["complevel"] == complevel for v in inputs)


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
        check=lambda t: check_output(idx, t),
        check_rerun=(
            (lambda t: check_compressed_input(idx, t, 5))
            if compare.value.startswith("mf6_export")
            else None
        ),
        compare=compare,
        htol=htol[idx],
    )
    test.run()


def build_output_models(idx, test, export):
    """The base simulation, writing its NetCDF output file."""
    sims = build_models(idx, test)
    model = sims[0].gwf[0]
    add_netcdf_output(model, export)
    model.dis.crs = wkt
    if export == "ugrid":
        chunks = {"chunk_face": 10}
    else:
        chunks = {"chunk_z": 1, "chunk_y": 5, "chunk_x": 5}
    flopy.mf6.ModflowUtlncf(
        model.dis,
        deflate=5,
        shuffle=True,
        chunk_time=1,
        wkt=wkt,
        filename=f"{model.name}.dis.ncf",
        **chunks,
    )
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwf[0]
    check_model_netcdf_output(model)
    grb = flopy.mf6.utils.MfGrdFile(Path(model.model_ws) / f"{model.name}.dis.grb")
    assert grb.version == 2
    assert grb.crs == wkt
    with nc.Dataset(Path(model.model_ws) / f"{model.name}.nc") as ds:
        head = ds.variables["head" if export == "structured" else "head_l1"]
        assert head.chunking() == ([1, 1, 5, 5] if export == "structured" else [1, 10])
        assert head.filters()["shuffle"]
        assert head.filters()["complevel"] == 5
        proj = ds.variables["projection"]
        assert proj.getncattr("wkt") == wkt
        assert proj.getncattr("crs_wkt") == wkt
        assert proj.getncattr("grid_mapping_name") == "transverse_mercator"


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
        htol=htol[idx],
    )
    test.run()
