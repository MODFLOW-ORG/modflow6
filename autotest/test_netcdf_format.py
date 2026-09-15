"""
NetCDF structured/mesh format, attribute, and CF/UGRID/MF6 convention
verification test.

Builds a minimal transient GWF model (DIS, IC, NPF, STO, CHD, WEL, OC) and
verifies that both OUTPUT and INPUT NetCDF files carry the expected convention
attributes for every combination of:

  fmt          : "structured" (DIS + NETCDF_STRUCTURED)
                 "ugrid"      (DIS + NETCDF_MESH2D)

  ncf_config   : "none"          - no utl-ncf package at all
                 "wkt_only"      - utl-ncf with wkt (WKT1) only
                 "crs_wkt_only"  - utl-ncf with crs_wkt (WKT2) only
                 "both"          - utl-ncf with distinct wkt (WKT1) and crs_wkt (WKT2)

  gridded_input: "ascii"  - only output (head) NC is checked
                 "netcdf" - additionally exports package arrays in validate mode
                            and checks the input NC attribute conventions

The test is focused on attribute presence and value correctness.  Data
round-trip correctness is covered by other tests.

Conventions
-----------
CF-1.13  (https://cfconventions.org/cf-conventions/cf-conventions.html)
  5.6  grid_mapping, crs_wkt, grid_mapping_name on the projection variable
  4.3  vertical (layer) coordinate: axis, positive, units
  4.4  time coordinate: standard_name, units, calendar
  7.3  cell_methods on the dependent variable ("time: point" -- an
       instantaneous state, not a time mean/accumulation)

UGRID-1.0  (https://ugrid-conventions.github.io/ugrid-conventions/)
  mesh topology variable; face-indexed data variables carrying mesh,
  location, and coordinates unconditionally (CRS-independent)

MF6 internal: modflow_model, modflow_grid, modflow_input
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

nc = pytest.importorskip("netCDF4")

# WKT strings — distinct WKT1 and WKT2 for EPSG:26918 (NAD83 / UTM zone 18N)
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

# Geographic WKT1 (EPSG:4269 NAD83) — no PROJECTION[ keyword;
# wkt_to_cf_gridmapping must return '' and grid_mapping_name must not be written.
WKT1_GEO = (
    'GEOGCS["NAD83",'
    'DATUM["North_American_Datum_1983",'
    'SPHEROID["GRS 1980",6378137,298.257222101]],'
    'PRIMEM["Greenwich",0],'
    'UNIT["degree",0.0174532925199433]]'
)

WKT2 = (
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

# small model appropriate for convention checks
NLAY, NROW, NCOL = 2, 3, 3
DELR = DELC = [100.0, 100.0, 100.0]
TOP = 20.0
BOTM = [10.0, 0.0]
STRT = 15.0
K = 1.0
SS = 1e-5
SY = 0.15
XORIGIN = 500_000.0
YORIGIN = 100_000_000.0
ANGROT_ROTATED = 15.0

# two stress periods: steady then transient
NPER = 2
PERIODDATA = [(1.0, 1, 1.0), (10.0, 3, 1.0)]  # (perlen, nstp, tsmult)

# CHD: fix head on first and last row of layer 1 — heads well within cell bounds
CHD_SP = {
    0: [((0, 0, j), 18.0) for j in range(NCOL)]
    + [((0, NROW - 1, j), 12.0) for j in range(NCOL)],
    1: [((0, 0, j), 18.0) for j in range(NCOL)]
    + [((0, NROW - 1, j), 12.0) for j in range(NCOL)],
}

# WEL: one well in layer 2, active only in transient period
WEL_SP = {1: [((1, 1, 1), -10.0)]}


def _has_crs(ncf_config):
    return ncf_config != "none"


def _input_nc_name(name, fmt):
    suffix = "structured" if fmt == "structured" else "ugrid"
    return f"{name}.{suffix}.nc"


def build_models(test, fmt, ncf_config):
    name = test.name
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name="mf6")

    flopy.mf6.ModflowTdis(
        sim,
        nper=NPER,
        perioddata=PERIODDATA,
        start_date_time="2000-01-01T00:00:00-00:00",
    )

    flopy.mf6.ModflowIms(sim, pname="ims", complexity="simple")

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        delr=DELR,
        delc=DELC,
        top=TOP,
        botm=BOTM,
        xorigin=XORIGIN,
        yorigin=YORIGIN,
        angrot=ANGROT_ROTATED if ncf_config == "rotated" else 0.0,
    )

    flopy.mf6.ModflowGwfic(gwf, strt=STRT)

    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=K)

    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=SS,
        sy=SY,
        steady_state={0: True},
        transient={1: True},
    )

    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=CHD_SP)

    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=WEL_SP)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )

    # configure NetCDF output format
    if fmt == "structured":
        gwf.name_file.nc_structured_filerecord = f"{name}.nc"
    elif fmt == "ugrid":
        gwf.name_file.nc_mesh2d_filerecord = f"{name}.nc"

    # configure utl-ncf when requested
    if ncf_config != "none":
        ncf_kwargs = dict(
            deflate=5,
            shuffle=True,
            chunk_time=1,
            filename=f"{name}.dis.ncf",
        )
        if fmt == "structured":
            ncf_kwargs.update(chunk_z=1, chunk_y=3, chunk_x=3)
        else:
            ncf_kwargs.update(chunk_face=9)

        if ncf_config == "wkt_only":
            ncf_kwargs["wkt"] = WKT1
        elif ncf_config == "crs_wkt_only":
            ncf_kwargs["crs_wkt"] = WKT2
        elif ncf_config == "both":
            ncf_kwargs["wkt"] = WKT1
            ncf_kwargs["crs_wkt"] = WKT2
        elif ncf_config == "geographic_wkt":
            ncf_kwargs["wkt"] = WKT1_GEO
        elif ncf_config == "rotated":
            ncf_kwargs["wkt"] = WKT1
        elif ncf_config == "latlon":
            # explicit ncpl works around a flopy bug: without it, flopy
            # infers NCPL as NROW instead of NROW*NCOL for structured
            # latitude/longitude griddata, corrupting the written .ncf file
            ncf_kwargs["ncpl"] = NROW * NCOL
            ncf_kwargs["latitude"] = np.linspace(40.0, 40.02, NROW * NCOL)
            ncf_kwargs["longitude"] = np.linspace(-75.0, -74.98, NROW * NCOL)

        flopy.mf6.ModflowUtlncf(dis, **ncf_kwargs)

    return sim, None


def _check_global_attrs(ds, name, fmt, ncf_config, label=""):
    ctx = f" [{label}]" if label else ""

    modflow_model = ds.getncattr("modflow_model")
    assert modflow_model == modflow_model.lower(), (
        f"modflow_model must be lowercase{ctx}: {modflow_model!r}"
    )
    assert "6:" not in modflow_model, (
        f"modflow_model must not contain '6:'{ctx}: {modflow_model!r}"
    )
    parts = modflow_model.split(":")
    assert len(parts) == 2, (
        f"modflow_model must be 'type: name'{ctx}: {modflow_model!r}"
    )
    assert parts[1].strip() == name.lower(), (
        f"modflow_model name mismatch{ctx}: {modflow_model!r}"
    )

    modflow_grid = ds.getncattr("modflow_grid")
    assert modflow_grid == modflow_grid.lower(), (
        f"modflow_grid must be lowercase{ctx}: {modflow_grid!r}"
    )
    assert len(modflow_grid) > 0, f"modflow_grid must be non-empty{ctx}"
    # For DIS structured format, expect "structured"; for DIS ugrid export the
    # value may differ — leave an explicit note for review
    if fmt == "structured":
        assert modflow_grid == "structured", (
            f"DIS structured modflow_grid must be 'structured'{ctx}: {modflow_grid!r}"
        )

    conventions = ds.getncattr("Conventions")
    assert "CF-1.13" in conventions, (
        f"Conventions must include CF-1.13{ctx}: {conventions!r}"
    )
    if fmt == "ugrid":
        assert "UGRID-1.0" in conventions, (
            f"UGRID Conventions must include UGRID-1.0{ctx}: {conventions!r}"
        )
    else:
        assert "UGRID" not in conventions, (
            f"Structured Conventions must not include UGRID{ctx}: {conventions!r}"
        )

    if fmt == "ugrid":
        assert "mesh" in ds.ncattrs(), f"mesh global attr missing{ctx}"
        assert ds.getncattr("mesh") == "layered", (
            f"mesh global attr must be 'layered'{ctx}"
        )
    else:
        assert "mesh" not in ds.ncattrs(), (
            f"Structured format must not carry mesh global attr{ctx}"
        )

    assert "source" in ds.ncattrs(), f"source global attr missing{ctx}"
    assert "history" in ds.ncattrs(), f"history global attr missing{ctx}"


def _check_layer_coord(ds, label=""):
    ctx = f" [{label}]" if label else ""
    assert "layer" in ds.variables, f"layer coordinate variable missing{ctx}"
    layer = ds.variables["layer"]
    assert layer.dtype == "int32", (
        f"layer must be stored as a 32-bit integer (NC_INT){ctx}, got {layer.dtype}"
    )
    assert layer.getncattr("axis") == "Z", f"layer axis must be 'Z'{ctx}"
    assert layer.getncattr("positive") == "down", f"layer positive must be 'down'{ctx}"
    assert layer.getncattr("units") == "1", f"layer units must be '1'{ctx}"
    assert layer.getncattr("long_name") == "model layer", (
        f"layer long_name must be 'model layer'{ctx}"
    )
    vals = layer[:]
    assert list(vals) == list(range(1, NLAY + 1)), (
        f"layer values must be 1..nlay{ctx}, got {list(vals)}"
    )


def _check_z_coord(ds, fmt, label=""):
    """z (structured) / z_l1, z_l2, ... (mesh/UGRID): unconditional
    auxiliary vertical position coordinate -- distinct from the discrete
    layer index. Split per layer for mesh since face-indexed variables
    carry no layer dimension to legally reference a
    combined z(layer, nmesh_face) coordinate (CF-1.13 5.2)."""
    ctx = f" [{label}]" if label else ""
    expected = [
        (TOP + BOTM[0]) / 2.0,
        (BOTM[0] + BOTM[1]) / 2.0,
    ]

    if fmt == "structured":
        assert "z" in ds.variables, f"z coordinate variable missing{ctx}"
        z = ds.variables["z"]
        assert z.dimensions == ("layer", "y", "x"), (
            f"z must be dimensioned (layer, y, x){ctx}, got {z.dimensions}"
        )
        assert "axis" not in z.ncattrs(), (
            f"z must not carry an axis attribute (layer remains the sole Z axis){ctx}"
        )
        assert z.getncattr("standard_name") == "altitude", (
            f"z standard_name must be 'altitude'{ctx}"
        )
        assert z.getncattr("positive") == "up", f"z positive must be 'up'{ctx}"
        assert z.getncattr("long_name") == "cell center elevation", (
            f"z long_name must be 'cell center elevation'{ctx}"
        )
        assert "units" in z.ncattrs(), f"z units missing{ctx}"

        vals = z[:]
        for k in range(NLAY):
            assert np.allclose(vals[k, :, :], expected[k]), (
                f"z layer {k + 1} value mismatch{ctx}: "
                f"got {vals[k, :, :]}, expected {expected[k]}"
            )
        return

    # mesh/UGRID: z_l1, z_l2, ... one per layer, dimensioned (nmesh_face,)
    assert "z" not in ds.variables, (
        f"a bare 'z' variable must not appear in mesh/UGRID output -- "
        f"mesh format uses per-layer z_lN{ctx}"
    )
    for k in range(1, NLAY + 1):
        vname = f"z_l{k}"
        assert vname in ds.variables, f"{vname} variable missing{ctx}"
        z = ds.variables[vname]
        assert z.dimensions == ("nmesh_face",), (
            f"{vname} must be dimensioned (nmesh_face,){ctx}, got {z.dimensions}"
        )
        assert "axis" not in z.ncattrs(), (
            f"{vname} must not carry an axis attribute{ctx}"
        )
        assert z.getncattr("standard_name") == "altitude", (
            f"{vname} standard_name must be 'altitude'{ctx}"
        )
        assert z.getncattr("positive") == "up", f"{vname} positive must be 'up'{ctx}"
        assert z.getncattr("long_name") == f"cell center elevation (layer {k})", (
            f"{vname} long_name must be 'cell center elevation (layer {k})'{ctx}"
        )
        assert "units" in z.ncattrs(), f"{vname} units missing{ctx}"
        assert "layer" in z.ncattrs(), f"layer attr missing on {vname}{ctx}"
        assert z.getncattr("layer") == k, (
            f"layer attr must be {k} on {vname}{ctx}: {z.getncattr('layer')}"
        )

        vals = z[:]
        assert np.allclose(vals, expected[k - 1]), (
            f"{vname} value mismatch{ctx}: got {vals}, expected {expected[k - 1]}"
        )


def _check_time_coord(ds, label=""):
    ctx = f" [{label}]" if label else ""
    assert "time" in ds.variables, f"time coordinate variable missing{ctx}"
    time = ds.variables["time"]
    assert "calendar" in time.ncattrs(), f"time calendar missing{ctx}"
    assert "units" in time.ncattrs(), f"time units missing{ctx}"
    assert time.getncattr("units").startswith("days since"), (
        f"time units must start with 'days since'{ctx}: {time.getncattr('units')!r}"
    )
    assert time.getncattr("axis") == "T", f"time axis must be 'T'{ctx}"
    assert time.getncattr("standard_name") == "time", (
        f"time standard_name must be 'time'{ctx}"
    )


def _check_projection(ds, fmt, ncf_config, label=""):
    ctx = f" [{label}]" if label else ""

    if not _has_crs(ncf_config):
        assert "projection" not in ds.variables, (
            f"projection variable must not be written without utl-ncf{ctx}"
        )
        return

    assert "projection" in ds.variables, (
        f"projection variable missing despite CRS being configured{ctx}"
    )
    proj = ds.variables["projection"]

    if ncf_config in ("wkt_only", "both"):
        assert "wkt" in proj.ncattrs(), f"wkt missing for {ncf_config!r}{ctx}"
        assert proj.getncattr("wkt") == WKT1, f"wkt value mismatch{ctx}"
    else:
        assert "wkt" not in proj.ncattrs(), (
            f"wkt must not be written for {ncf_config!r}{ctx}"
        )

    assert "crs_wkt" in proj.ncattrs(), f"crs_wkt missing{ctx}"
    if ncf_config == "crs_wkt_only":
        assert proj.getncattr("crs_wkt") == WKT2, f"crs_wkt must be WKT2{ctx}"
    elif ncf_config == "both":
        assert proj.getncattr("crs_wkt") == WKT2, (
            f"crs_wkt must be explicit WKT2 when both provided{ctx}"
        )
    elif ncf_config == "wkt_only":
        assert proj.getncattr("crs_wkt") == WKT1, (
            f"crs_wkt must fall back to wkt when only wkt provided{ctx}"
        )

    assert "grid_mapping_name" in proj.ncattrs(), f"grid_mapping_name missing{ctx}"
    assert proj.getncattr("grid_mapping_name") == "transverse_mercator", (
        f"expected transverse_mercator for UTM18N{ctx}: "
        f"{proj.getncattr('grid_mapping_name')!r}"
    )

    # GeoTransform / spatial_ref: structured-only, requires WKT1
    # (this%wkt) since MF6 has no CRS library to derive WKT1 from WKT2, and
    # only for unrotated grids (this model has ANGROT=0).
    expect_geotransform = fmt == "structured" and ncf_config in ("wkt_only", "both")
    if expect_geotransform:
        assert "GeoTransform" in proj.ncattrs(), (
            f"GeoTransform missing for {ncf_config!r}{ctx}"
        )
        gt = [float(v) for v in proj.getncattr("GeoTransform").split()]
        assert len(gt) == 6, f"GeoTransform must have 6 values{ctx}: {gt}"
        expected_gt = [
            XORIGIN,
            sum(DELR) / NCOL,
            0.0,
            YORIGIN + sum(DELC),
            0.0,
            -sum(DELC) / NROW,
        ]
        assert np.allclose(gt, expected_gt), (
            f"GeoTransform value mismatch{ctx}: {gt} != {expected_gt}"
        )
        assert "spatial_ref" in proj.ncattrs(), f"spatial_ref missing{ctx}"
        assert proj.getncattr("spatial_ref") == WKT1, f"spatial_ref must be WKT1{ctx}"
    else:
        assert "GeoTransform" not in proj.ncattrs(), (
            f"GeoTransform must not be written for fmt={fmt!r}, "
            f"ncf_config={ncf_config!r}{ctx}"
        )
        assert "spatial_ref" not in proj.ncattrs(), (
            f"spatial_ref must not be written for fmt={fmt!r}, "
            f"ncf_config={ncf_config!r}{ctx}"
        )


def _check_coord_gridmapping(ds, fmt, ncf_config, label=""):
    """Verify that coordinate variables carry grid_mapping when CRS is set."""
    ctx = f" [{label}]" if label else ""

    if not _has_crs(ncf_config):
        return

    if fmt == "structured":
        for coord in ("x", "y"):
            assert "grid_mapping" in ds.variables[coord].ncattrs(), (
                f"{coord} coordinate missing grid_mapping "
                f"(ncf_config={ncf_config!r}){ctx}"
            )
            assert ds.variables[coord].getncattr("grid_mapping") == "projection", (
                f"{coord} grid_mapping must be 'projection'{ctx}"
            )
    else:
        for coord in (
            "mesh_node_x",
            "mesh_node_y",
            "mesh_face_x",
            "mesh_face_y",
        ):
            assert "grid_mapping" in ds.variables[coord].ncattrs(), (
                f"{coord} missing grid_mapping (ncf_config={ncf_config!r}){ctx}"
            )
            assert ds.variables[coord].getncattr("grid_mapping") == "projection", (
                f"{coord} grid_mapping must be 'projection'{ctx}"
            )


def _check_data_var(var, vname, fmt, ncf_config, label=""):
    """Common assertions for any data variable (output or input)."""
    ctx = f" [{label}]" if label else ""

    assert "long_name" in var.ncattrs(), f"long_name missing on {vname}{ctx}"
    assert len(var.getncattr("long_name")) > 0, (
        f"long_name must be non-empty on {vname}{ctx}"
    )

    if _has_crs(ncf_config):
        assert "grid_mapping" in var.ncattrs(), (
            f"grid_mapping missing on {vname} (ncf_config={ncf_config!r}){ctx}"
        )
        assert var.getncattr("grid_mapping") == "projection", (
            f"grid_mapping must be 'projection' on {vname}{ctx}"
        )
    else:
        assert "grid_mapping" not in var.ncattrs(), (
            f"grid_mapping must not be written on {vname} without CRS{ctx}"
        )

    # z: structured only, unconditional -- present in
    # coordinates iff the variable actually carries the layer dimension
    # (CF-1.13 5.2 subset rule; z's own dims are (layer, y, x)).
    if fmt == "structured":
        has_layer = "layer" in var.dimensions
        coords = (
            var.getncattr("coordinates").split()
            if "coordinates" in var.ncattrs()
            else []
        )
        if has_layer:
            assert "z" in coords, (
                f"z missing from coordinates on layered var {vname}{ctx}: {coords}"
            )
        else:
            assert "z" not in coords, (
                f"z must not appear in coordinates on non-layered var "
                f"{vname}{ctx}: {coords}"
            )
        # x/y are never listed here, CRS or not -- they are true CF
        # dimension coordinates, already discoverable by dimension-name
        # matching alone (CF-1.13 Ch.5 preamble), unlike lon/lat (real 2D
        # auxiliary coordinates, which do need to be listed).
        assert "x" not in coords and "y" not in coords, (
            f"x/y must never appear in coordinates on {vname}{ctx}: {coords}"
        )
    else:
        # z_lN: the layer attribute (not the variable name) identifies
        # whether vname is a per-layer-split variable, and which layer it
        # belongs to
        is_layered = "layer" in var.ncattrs()
        coords = (
            var.getncattr("coordinates").split()
            if "coordinates" in var.ncattrs()
            else []
        )
        if is_layered:
            expected_z = f"z_l{var.getncattr('layer')}"
            assert expected_z in coords, (
                f"{expected_z} missing from coordinates on {vname}{ctx}: {coords}"
            )
        else:
            assert not any(c.startswith("z_l") for c in coords), (
                f"no z_lN may appear in coordinates on non-layered var "
                f"{vname}{ctx}: {coords}"
            )


def _check_output_nc(ds, name, fmt, ncf_config):
    _check_global_attrs(ds, name, fmt, ncf_config, label="output")
    _check_layer_coord(ds, label="output")
    _check_z_coord(ds, fmt, label="output")
    _check_time_coord(ds, label="output")
    _check_projection(ds, fmt, ncf_config, label="output")
    _check_coord_gridmapping(ds, fmt, ncf_config, label="output")

    if fmt == "structured":
        assert "head" in ds.variables, "head variable missing from structured output"
        _check_data_var(ds.variables["head"], "head", fmt, ncf_config, label="output")
        # cell_methods (item 8): dependent variable is an instantaneous
        # state at each output time, not a time mean/accumulation
        assert ds.variables["head"].getncattr("cell_methods") == "time: point", (
            "head cell_methods must be 'time: point'"
        )
    else:
        for k in range(1, NLAY + 1):
            vname = f"head_l{k}"
            assert vname in ds.variables, f"{vname} missing from UGRID output"
            head = ds.variables[vname]
            _check_data_var(head, vname, fmt, ncf_config, label="output")

            assert head.getncattr("cell_methods") == "time: point", (
                f"{vname} cell_methods must be 'time: point'"
            )

            # layer attribute must match this variable's actual layer, not
            # just be present -- independent check against the known k,
            # not derived from the attribute itself
            assert "layer" in head.ncattrs(), f"layer attr missing on {vname}"
            assert head.getncattr("layer") == k, (
                f"layer attr must be {k} on {vname}: {head.getncattr('layer')}"
            )

            # UGRID topology attrs always present on face vars regardless of CRS
            assert head.getncattr("mesh") == "mesh", f"{vname} mesh attr must be 'mesh'"
            assert head.getncattr("location") == "face", (
                f"{vname} location attr must be 'face'"
            )
            assert "coordinates" in head.ncattrs(), f"{vname} coordinates attr missing"

    if ncf_config != "none":
        sample = "head" if fmt == "structured" else "head_l1"
        cmpr = ds.variables[sample].filters()
        assert cmpr["complevel"] == 5, "expected deflate level 5"
        assert cmpr["shuffle"], "expected shuffle enabled"


# Structured input vars: split by whether they have full 2D spatial extent.
# grid_mapping is only meaningful on variables whose dimensions include the full
# horizontal spatial extent (nrow x ncol or nlay x nrow x ncol). 1D dimension arrays
# (delr=ncol, delc=nrow) define grid geometry but are not georeferenced fields,
# so grid_mapping does not apply — consistent with CF-1.13 5.6.
_STRUCTURED_VARS_NONSPATIAL = [
    "dis_delr",  # 1D (ncol) — grid spacing, not a georeferenced field
    "dis_delc",  # 1D (nrow) — grid spacing, not a georeferenced field
]
_STRUCTURED_VARS_SPATIAL = [
    "dis_top",  # 2D (nrow, ncol)
    "dis_botm",  # 3D (nlay, nrow, ncol)
    "ic_strt",  # 3D
    "npf_icelltype",  # 3D
    "npf_k",  # 3D
    "sto_iconvert",  # 3D
    "sto_ss",  # 3D
    "sto_sy",  # 3D
]

# UGRID: all face variables are indexed by ncpl — inherently 2D in geographic
# space — so grid_mapping applies to all of them, including single-layer vars.
_UGRID_VARS_1D = [
    "dis_top",  # 1D (ncpl)
]
_UGRID_VARS_LAYERED = [
    "dis_botm",
    "ic_strt",
    "npf_icelltype",
    "npf_k",
    "sto_iconvert",
    "sto_ss",
    "sto_sy",
]


def _check_input_var(
    var,
    vname,
    name,
    fmt,
    ncf_config,
    is_layered,
    expected_layer,
    spatial=True,
    label="",
):
    """Check convention attrs on a single input NC variable."""
    ctx = f" [{label}]" if label else ""

    # long_name
    assert "long_name" in var.ncattrs(), f"long_name missing on {vname}{ctx}"
    assert len(var.getncattr("long_name")) > 0, (
        f"long_name must be non-empty on {vname}{ctx}"
    )

    # modflow_input: modelname/pkgtype_or_instname/paramtag (all lowercase, "/" sep)
    assert "modflow_input" in var.ncattrs(), f"modflow_input missing on {vname}{ctx}"
    mi = var.getncattr("modflow_input")
    assert mi == mi.lower(), f"modflow_input must be lowercase on {vname}{ctx}: {mi!r}"
    parts = mi.split("/")
    assert len(parts) == 3, (
        f"modflow_input must have 3 slash-separated parts on {vname}{ctx}: {mi!r}"
    )
    assert parts[0] == name.lower(), (
        f"modflow_input model name must match model name on {vname}{ctx}: {mi!r}"
    )
    assert all(len(p) > 0 for p in parts), (
        f"modflow_input parts must all be non-empty on {vname}{ctx}: {mi!r}"
    )

    # grid_mapping: only check on 2D/3D spatial fields (ncpl or nrow x ncol extent),
    # not on 1D dimension arrays like dis_delr/dis_delc
    if spatial:
        _check_data_var(var, vname, fmt, ncf_config, label=label)

    # UGRID-specific: per-layer face variables
    if fmt == "ugrid" and is_layered:
        # layer integer attr (wiki: input layer vars carry layer = k)
        assert "layer" in var.ncattrs(), (
            f"layer attr missing on UGRID layered input var {vname}{ctx}"
        )
        assert var.getncattr("layer") == expected_layer, (
            f"layer attr must be {expected_layer} on {vname}{ctx}: "
            f"{var.getncattr('layer')}"
        )

        # UGRID topology linkage always present on face vars regardless of CRS
        assert "mesh" in var.ncattrs(), f"mesh attr missing on {vname}{ctx}"
        assert var.getncattr("mesh") == "mesh", (
            f"mesh attr must be 'mesh' on {vname}{ctx}"
        )
        assert "location" in var.ncattrs(), f"location attr missing on {vname}{ctx}"
        assert var.getncattr("location") == "face", (
            f"location attr must be 'face' on {vname}{ctx}"
        )
        assert "coordinates" in var.ncattrs(), (
            f"coordinates attr missing on {vname}{ctx}"
        )


def _check_input_nc(ds, name, fmt, ncf_config):
    _check_global_attrs(ds, name, fmt, ncf_config, label="input")
    _check_layer_coord(ds, label="input")
    _check_z_coord(ds, fmt, label="input")
    _check_projection(ds, fmt, ncf_config, label="input")
    _check_coord_gridmapping(ds, fmt, ncf_config, label="input")
    # time coord only present on period-data variables; don't require it on input

    if fmt == "structured":
        for vname in _STRUCTURED_VARS_NONSPATIAL:
            assert vname in ds.variables, (
                f"Expected input variable '{vname}' missing from structured input NC"
            )
            _check_input_var(
                ds.variables[vname],
                vname,
                name,
                fmt,
                ncf_config,
                is_layered=False,
                expected_layer=None,
                spatial=False,
                label="input",
            )
        for vname in _STRUCTURED_VARS_SPATIAL:
            assert vname in ds.variables, (
                f"Expected input variable '{vname}' missing from structured input NC"
            )
            _check_input_var(
                ds.variables[vname],
                vname,
                name,
                fmt,
                ncf_config,
                is_layered=False,
                expected_layer=None,
                spatial=True,
                label="input",
            )
    else:
        for vname in _UGRID_VARS_1D:
            assert vname in ds.variables, (
                f"Expected input variable '{vname}' missing from UGRID input NC"
            )
            _check_input_var(
                ds.variables[vname],
                vname,
                name,
                fmt,
                ncf_config,
                is_layered=False,
                expected_layer=None,
                label="input",
            )
        for base in _UGRID_VARS_LAYERED:
            for k in range(1, NLAY + 1):
                vname = f"{base}_l{k}"
                assert vname in ds.variables, (
                    f"Expected input variable '{vname}' missing from UGRID input NC"
                )
                _check_input_var(
                    ds.variables[vname],
                    vname,
                    name,
                    fmt,
                    ncf_config,
                    is_layered=True,
                    expected_layer=k,
                    label="input",
                )


def check_output(test, fmt, ncf_config, gridded_input):
    name = test.name
    ws = test.workspace

    # --- output NC ---
    with nc.Dataset(ws / f"{name}.nc") as ds:
        assert ds.data_model == "NETCDF4"
        _check_output_nc(ds, name, fmt, ncf_config)

    # --- input NC (validate mode) ---
    if gridded_input == "ascii":
        return

    gwf = test.sims[0].gwf[0]
    for pkg_name in ("DIS", "IC", "NPF", "STO"):
        gwf.get_package(pkg_name).export_array_netcdf = True

    input_nc = _input_nc_name(name, fmt)
    if fmt == "structured":
        gwf.name_file.nc_structured_filerecord = input_nc
    else:
        gwf.name_file.nc_mesh2d_filerecord = input_nc

    test.sims[0].write_simulation()
    success, _ = flopy.run_model(
        test.targets["mf6"],
        ws / "mfsim.nam",
        model_ws=ws,
        report=True,
        cargs=["--mode=validate"],
    )
    assert success, "validate-mode run failed"

    with nc.Dataset(ws / input_nc) as ds:
        assert ds.data_model == "NETCDF4"
        _check_input_nc(ds, name, fmt, ncf_config)


def _check_geo_crs_output(test, fmt):
    """projection variable present, crs_wkt written, grid_mapping_name absent."""
    name = test.name
    ws = test.workspace
    with nc.Dataset(ws / f"{name}.nc") as ds:
        _check_global_attrs(ds, name, fmt, "geographic_wkt", label="geo_output")
        _check_layer_coord(ds, label="geo_output")
        _check_z_coord(ds, fmt, label="geo_output")
        if fmt == "structured":
            # grid_mapping name (the variable "projection" is still pointed
            # to) is independent of grid_mapping_name (the CF attribute
            # value on that variable, empty here since geographic CRS has no
            # PROJECTION[ keyword) -- head still gets grid_mapping + "z"
            # (x/y are never listed in coordinates, CRS or not).
            _check_data_var(
                ds.variables["head"], "head", fmt, "geographic_wkt", label="geo_output"
            )
        assert "projection" in ds.variables, (
            "projection variable missing for geographic CRS"
        )
        proj = ds.variables["projection"]
        # crs_wkt falls back to wkt when no explicit crs_wkt provided
        assert "crs_wkt" in proj.ncattrs(), "crs_wkt missing for geographic CRS"
        assert proj.getncattr("crs_wkt") == WKT1_GEO, "crs_wkt must equal the input wkt"
        # geographic CRS has no PROJECTION[ keyword → wkt_to_cf_gridmapping returns ''
        assert "grid_mapping_name" not in proj.ncattrs(), (
            "grid_mapping_name must not be written for geographic (unprojected) CRS"
        )


def _check_latlon_output(test):
    """latitude/longitude griddata: coordinates = "lon lat z", no grid_mapping."""
    name = test.name
    ws = test.workspace
    with nc.Dataset(ws / f"{name}.nc") as ds:
        _check_global_attrs(ds, name, "structured", "latlon", label="latlon_output")
        _check_layer_coord(ds, label="latlon_output")
        _check_z_coord(ds, "structured", label="latlon_output")

        assert "lat" in ds.variables, "lat variable missing"
        assert "lon" in ds.variables, "lon variable missing"
        lat = ds.variables["lat"]
        lon = ds.variables["lon"]
        assert lat.getncattr("standard_name") == "latitude", (
            "lat standard_name must be 'latitude'"
        )
        assert lon.getncattr("standard_name") == "longitude", (
            "lon standard_name must be 'longitude'"
        )
        assert lat.getncattr("units") == "degrees_north", (
            "lat units must be 'degrees_north'"
        )
        assert lon.getncattr("units") == "degrees_east", (
            "lon units must be 'degrees_east'"
        )

        assert "projection" not in ds.variables, (
            "projection variable must not be written for latlon config"
        )
        head = ds.variables["head"]
        assert "grid_mapping" not in head.ncattrs(), (
            "grid_mapping must not be written for latlon config"
        )
        coords = head.getncattr("coordinates").split()
        assert "lon" in coords and "lat" in coords, (
            f"lon/lat missing from head coordinates: {coords}"
        )
        assert "z" in coords, f"z missing from head coordinates: {coords}"
        assert "x" not in coords and "y" not in coords, (
            f"x/y must not appear in coordinates for latlon config: {coords}"
        )


def _rotated_geotransform_corners(gt):
    """Real-world (x, y) at the four pixel-grid corners implied by a
    6-element GDAL GeoTransform."""
    corners = []
    for xpixel in (0, NCOL):
        for yline in (0, NROW):
            xg = gt[0] + xpixel * gt[1] + yline * gt[2]
            yg = gt[3] + xpixel * gt[4] + yline * gt[5]
            corners.append((xg, yg))
    return corners


def _check_rotated_output(test, fmt):
    """Rotated structured grid: GeoTransform must be present and correct
    (rotation-aware formula, DELR/DELC uniform in this model so the
    effective-pixel-size approximation is exact, not just close). x/y
    dimension coordinates remain grid-local regardless (CF dimension
    coordinates cannot represent a rotated position) -- unaffected by this
    change, not re-checked here."""
    name = test.name
    ws = test.workspace
    with nc.Dataset(ws / f"{name}.nc") as ds:
        proj = ds.variables["projection"]
        if fmt == "structured":
            assert "GeoTransform" in proj.ncattrs(), (
                "GeoTransform missing for rotated structured grid"
            )
            gt = [float(v) for v in proj.getncattr("GeoTransform").split()]
            assert len(gt) == 6, f"GeoTransform must have 6 values: {gt}"

            ang = np.radians(ANGROT_ROTATED)
            dx_eff = sum(DELR) / NCOL
            dy_eff = -sum(DELC) / NROW
            expected_gt = [
                XORIGIN - sum(DELC) * np.sin(ang),
                dx_eff * np.cos(ang),
                -dy_eff * np.sin(ang),
                YORIGIN + sum(DELC) * np.cos(ang),
                dx_eff * np.sin(ang),
                dy_eff * np.cos(ang),
            ]
            assert np.allclose(gt, expected_gt), (
                f"rotated GeoTransform value mismatch: {gt} != {expected_gt}"
            )
            assert "spatial_ref" in proj.ncattrs(), "spatial_ref missing"
            assert proj.getncattr("spatial_ref") == WKT1, "spatial_ref must be WKT1"
        else:
            # mesh format already handled rotation correctly before this
            # change (via dis_transform_xy) -- confirm unaffected/unchanged
            assert "GeoTransform" not in proj.ncattrs(), (
                "GeoTransform must never be written for mesh format"
            )
            assert "spatial_ref" not in proj.ncattrs(), (
                "spatial_ref must never be written for mesh format"
            )


def _check_rotated_cross_format(struct_test, ugrid_test):
    """Strong cross-check: verify structured's new GeoTransform-implied grid
    corners (a newly-Fortran-implemented rotation formula) against mesh
    format's independently-implemented, already-shipped, real-world node
    coordinates (dis_transform_xy) -- both generated by real MF6 runs of the
    identical rotated grid, not just checked against a Python derivation."""
    with nc.Dataset(struct_test.workspace / f"{struct_test.name}.nc") as ds:
        gt = [
            float(v)
            for v in ds.variables["projection"].getncattr("GeoTransform").split()
        ]
    expected_corners = _rotated_geotransform_corners(gt)

    with nc.Dataset(ugrid_test.workspace / f"{ugrid_test.name}.nc") as ds:
        node_x = ds.variables["mesh_node_x"][:]
        node_y = ds.variables["mesh_node_y"][:]

    # tolerance reflects GeoTransform's ES16.8 text-formatting precision
    # (9 significant digits) at YORIGIN's ~1e8 magnitude, i.e. ~0.1 unit --
    # not the underlying formula's precision, which matches exactly (see
    # test_rotated_geotransform's direct np.allclose check on the raw GT
    # values, unaffected by this compounding-corner-arithmetic concern)
    for xg, yg in expected_corners:
        dist = np.sqrt((node_x - xg) ** 2 + (node_y - yg) ** 2)
        assert dist.min() < 1.0, (
            f"structured GeoTransform corner ({xg}, {yg}) has no matching "
            f"mesh node (closest distance {dist.min()!r}) -- structured "
            f"and mesh rotation implementations disagree"
        )


cases = ["gwf_cf_conv"]


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("fmt", ["structured", "ugrid"])
@pytest.mark.parametrize("ncf_config", ["none", "wkt_only", "crs_wkt_only", "both"])
@pytest.mark.parametrize("gridded_input", ["ascii", "netcdf"])
def test_mf6model(idx, name, function_tmpdir, targets, fmt, ncf_config, gridded_input):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, fmt, ncf_config),
        check=lambda t: check_output(t, fmt, ncf_config, gridded_input),
        targets=targets,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
@pytest.mark.parametrize("fmt", ["structured", "ugrid"])
def test_geographic_crs_no_grid_mapping_name(fmt, function_tmpdir, targets):
    """Geographic CRS: projection variable written but grid_mapping_name absent."""
    test = TestFramework(
        name="gwf_geo_crs",
        workspace=function_tmpdir,
        build=lambda t: build_models(t, fmt, "geographic_wkt"),
        check=lambda t: _check_geo_crs_output(t, fmt),
        targets=targets,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
def test_latlon_griddata(function_tmpdir, targets):
    """latitude/longitude griddata (structured only): coordinates = "lon lat z",
    no grid_mapping/projection."""
    test = TestFramework(
        name="gwf_latlon",
        workspace=function_tmpdir,
        build=lambda t: build_models(t, "structured", "latlon"),
        check=lambda t: _check_latlon_output(t),
        targets=targets,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
@pytest.mark.parametrize("fmt", ["structured", "ugrid"])
def test_rotated_geotransform(fmt, function_tmpdir, targets):
    """Rotated grid (ANGROT != 0) with a CRS configured: structured now gets
    a rotation-aware GeoTransform (reopened 2026-07-09 -- previously
    withheld entirely for rotated grids); mesh format is unaffected (it
    already handled rotation correctly via dis_transform_xy)."""
    test = TestFramework(
        name="gwf_rotated",
        workspace=function_tmpdir,
        build=lambda t: build_models(t, fmt, "rotated"),
        check=lambda t: _check_rotated_output(t, fmt),
        targets=targets,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
def test_rotated_geotransform_matches_mesh(function_tmpdir, targets):
    """Strong cross-check: build the identical rotated grid in both
    structured and mesh format, and confirm structured's new
    GeoTransform-implied corner coordinates match mesh's independently
    computed, already-shipped real-world node coordinates -- validates the
    new Fortran implementation against a second, trusted code path, not
    just a Python-side derivation."""
    (function_tmpdir / "structured").mkdir()
    (function_tmpdir / "ugrid").mkdir()

    struct_test = TestFramework(
        name="gwf_rotated",
        workspace=function_tmpdir / "structured",
        build=lambda t: build_models(t, "structured", "rotated"),
        check=lambda t: _check_rotated_output(t, "structured"),
        targets=targets,
        compare=None,
    )
    struct_test.run()

    ugrid_test = TestFramework(
        name="gwf_rotated",
        workspace=function_tmpdir / "ugrid",
        build=lambda t: build_models(t, "ugrid", "rotated"),
        check=lambda t: _check_rotated_output(t, "ugrid"),
        targets=targets,
        compare=None,
    )
    ugrid_test.run()

    _check_rotated_cross_format(struct_test, ugrid_test)
