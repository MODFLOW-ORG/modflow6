"""
NetCDF input tests for test_gwf_disv: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4.
"""

from pathlib import Path

import flopy
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwf_disv import build_models, cases, check_output

nc = pytest.importorskip("netCDF4")

wkt1 = (
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

wkt2 = (
    'PROJCRS["NAD83 / UTM zone 18N",'
    'BASEGEOGCRS["NAD83",'
    'DATUM["North American Datum 1983",'
    'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
    'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],'
    'ID["EPSG",4269]],'
    'CONVERSION["UTM zone 18N",'
    'METHOD["Transverse Mercator",ID["EPSG",9807]],'
    'PARAMETER["Latitude of natural origin",0,'
    'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],'
    'PARAMETER["Longitude of natural origin",-75,'
    'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
    'PARAMETER["Scale factor at natural origin",0.9996,'
    'SCALEUNIT["unity",1],ID["EPSG",8805]],'
    'PARAMETER["False easting",500000,LENGTHUNIT["metre",1],ID["EPSG",8806]],'
    'PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],'
    "CS[Cartesian,2],"
    'AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1]],'
    'AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1]],'
    'ID["EPSG",26918]]'
)


def build_input_models(idx, test, compare):
    """The base simulation, saving head and budget; for MF6's own export,
    with its arrays exported compressed and chunked."""
    sims = build_models(idx, test)
    model = sims[0].gwf[0]
    # the base model saves no output for the re-run to be compared on
    flopy.mf6.ModflowGwfoc(
        model,
        head_filerecord=f"{model.name}.hds",
        budget_filerecord=f"{model.name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    if compare == Comparison.MF6_EXPORT_LAYERED:
        flopy.mf6.ModflowUtlncf(
            model.disv,
            deflate=9,
            shuffle=True,
            chunk_time=1,
            chunk_face=3,
            crs_wkt=wkt2,
            filename=f"{model.name}.disv.ncf",
        )
    return sims


def check_compressed_input(idx, test, complevel):
    """The base checks, on a re-run reading compressed input, and the
    re-run's NetCDF output file written with the same NCF options."""
    check_output(idx, test)
    for path in test.workspace.glob("*.export.nc"):
        with nc.Dataset(path) as ds:
            inputs = [
                v for v in ds.variables.values() if "modflow_input" in v.ncattrs()
            ]
            assert inputs and all(v.filters()["complevel"] == complevel for v in inputs)
    model = test.sims[0].gwf[0]
    with nc.Dataset(test.workspace / f"{model.name}.nc") as ds:
        head = ds.variables["head_l1"]
        assert head.filters()["shuffle"]
        assert head.filters()["complevel"] == complevel
        assert head.chunking() == [1, 3]
        proj = ds.variables["projection"]
        assert proj.getncattr("crs_wkt") == wkt2
        assert proj.getncattr("grid_mapping_name") == "transverse_mercator"


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare",
    [
        Comparison.MF6_EXPORT_LAYERED,
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
            (lambda t: check_compressed_input(idx, t, 9))
            if compare.value.startswith("mf6_export")
            else None
        ),
        compare=compare,
    )
    test.run()


def build_output_models(idx, test, export, wkt_mode):
    """The base simulation, writing its NetCDF output file."""
    sims = build_models(idx, test)
    sims[0].tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    model = sims[0].gwf[0]
    add_netcdf_output(model, export)
    # the base model saves no head
    flopy.mf6.ModflowGwfoc(
        model, head_filerecord=f"{model.name}.hds", saverecord=[("HEAD", "ALL")]
    )
    if wkt_mode == "both":
        # distinct WKT1 and WKT2, so crs_wkt doesn't fall back to wkt
        model.disv.crs = wkt1
        crs = {"wkt": wkt1, "crs_wkt": wkt2}
    else:
        # grid_mapping_name derived from the WKT2 METHOD
        crs = {"crs_wkt": wkt2}
    flopy.mf6.ModflowUtlncf(
        model.disv,
        deflate=9,
        shuffle=True,
        chunk_time=1,
        chunk_face=3,
        filename=f"{model.name}.disv.ncf",
        **crs,
    )
    return sims


def check_output_file(idx, test, export, wkt_mode):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwf[0]
    check_model_netcdf_output(model)
    if wkt_mode == "both":
        grb = flopy.mf6.utils.MfGrdFile(Path(model.model_ws) / "disv.grb")
        assert grb.version == 2
        assert grb.crs == wkt1
    with nc.Dataset(Path(model.model_ws) / f"{model.name}.nc") as ds:
        head = ds.variables["head_l1"]
        assert head.filters()["shuffle"]
        assert head.filters()["complevel"] == 9
        assert head.chunking() == [1, 3]
        proj = ds.variables["projection"]
        if wkt_mode == "both":
            assert proj.getncattr("wkt") == wkt1
        else:
            assert "wkt" not in proj.ncattrs()
        assert proj.getncattr("crs_wkt") == wkt2
        assert proj.getncattr("grid_mapping_name") == "transverse_mercator"


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid"])
@pytest.mark.parametrize("wkt_mode", ["both", "wkt2_only"])
def test_netcdf_output(idx, name, export, wkt_mode, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_output_models(idx, t, export, wkt_mode),
        check=lambda t: check_output_file(idx, t, export, wkt_mode),
        compare=None,
    )
    test.run()
