"""
NetCDF export test version of test_gwf_disv.  This test compares
the heads and input arrays in the the NetCDF file to those
in the FloPy binary output head file and package data objects.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_disv import cases

xa = pytest.importorskip("xarray")
xu = pytest.importorskip("xugrid")
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


def build_models(idx, test, gridded_input):
    from test_gwf_disv import build_models as build

    sim, dummy = build(idx, test)
    sim.tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    gwf = sim.gwf[0]
    gwf.disv.export_array_netcdf = True
    gwf.disv.crs = wkt
    gwf.ic.export_array_netcdf = True
    gwf.npf.export_array_netcdf = True

    name = cases[idx]

    gwf.name_file.nc_mesh2d_filerecord = (
        f"{name}.ugrid.nc" if gridded_input == "netcdf" else f"{name}.nc"
    )

    # netcdf config
    ncf = flopy.mf6.ModflowUtlncf(
        gwf.disv,
        deflate=9,
        shuffle=True,
        chunk_time=1,
        chunk_face=3,
        wkt=wkt,
        filename=f"{name}.disv.ncf",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[
            ("HEAD", "ALL"),
        ],
    )

    return sim, dummy


def check_output(idx, test, gridded_input):
    from test_gwf_disv import check_output as check

    name = test.name

    # verify crs data string in grb version 2 file
    fname = os.path.join(test.workspace, "disv.grb")
    grbobj = flopy.mf6.utils.MfGrdFile(fname)
    crs = grbobj._datadict["CRS"]
    assert crs == wkt

    if gridded_input == "netcdf":
        # re-run the simulation with model netcdf input
        with open(test.workspace / f"{name}.nam", "w") as f:
            f.write("BEGIN options\n")
            f.write(f"  NETCDF_MESH2D  FILEOUT  {name}.nc\n")
            f.write(f"  NETCDF  FILEIN {name}.ugrid.nc\n")
            f.write("END options\n\n")
            f.write("BEGIN packages\n")
            f.write(f"  DISV6  {name}.disv  disv\n")
            f.write(f"  IC6  {name}.ic  ic\n")
            f.write(f"  NPF6  {name}.npf  npf\n")
            f.write(f"  CHD6  {name}.chd  chd_0\n")
            f.write(f"  OC6  {name}.oc  oc\n")
            f.write("END packages\n")

        with open(test.workspace / f"{name}.disv", "w") as f:
            f.write("BEGIN options\n")
            f.write(f"  NCF6  FILEIN  {name}.disv.ncf\n")
            f.write("END options\n\n")
            f.write("BEGIN dimensions\n")
            f.write("  NLAY  3\n")
            f.write("  NCPL  9\n")
            f.write("  NVERT  16\n")
            f.write("END dimensions\n\n")
            f.write("BEGIN griddata\n")
            f.write("  top NETCDF\n")
            f.write("  botm NETCDF\n")
            if name == "disv01b":
                f.write("  idomain NETCDF\n")
            f.write("END griddata\n\n")
            f.write("BEGIN vertices\n")
            f.write("  1  1.00000000E+08  1.00000030E+08\n")
            f.write("  2  1.00000010E+08  1.00000030E+08\n")
            f.write("  3  1.00000020E+08  1.00000030E+08\n")
            f.write("  4  1.00000030E+08  1.00000030E+08\n")
            f.write("  5  1.00000000E+08  1.00000020E+08\n")
            f.write("  6  1.00000010E+08  1.00000020E+08\n")
            f.write("  7  1.00000020E+08  1.00000020E+08\n")
            f.write("  8  1.00000030E+08  1.00000020E+08\n")
            f.write("  9  1.00000000E+08  1.00000010E+08\n")
            f.write("  10  1.00000010E+08  1.00000010E+08\n")
            f.write("  11  1.00000020E+08  1.00000010E+08\n")
            f.write("  12  1.00000030E+08  1.00000010E+08\n")
            f.write("  13  1.00000000E+08  1.00000000E+08\n")
            f.write("  14  1.00000010E+08  1.00000000E+08\n")
            f.write("  15  1.00000020E+08  1.00000000E+08\n")
            f.write("  16  1.00000030E+08  1.00000000E+08\n")
            f.write("END vertices\n\n")
            f.write("BEGIN cell2d\n")
            f.write("  1  1.00000005E+08  1.00000025E+08  4  1  2  6  5\n")
            f.write("  2  1.00000015E+08  1.00000025E+08  4  2  3  7  6\n")
            f.write("  3  1.00000025E+08  1.00000025E+08  4  3  4  8  7\n")
            f.write("  4  1.00000005E+08  1.00000015E+08  4  5  6  10  9\n")
            f.write("  5  1.00000015E+08  1.00000015E+08  4  6  7  11  10\n")
            f.write("  6  1.00000025E+08  1.00000015E+08  4  7  8  12  11\n")
            f.write("  7  1.00000005E+08  1.00000005E+08  4  9  10  14  13\n")
            f.write("  8  1.00000015E+08  1.00000005E+08  4  10  11  15  14\n")
            f.write("  9  1.00000025E+08  1.00000005E+08  4  11  12  16  15\n")
            f.write("END cell2d\n\n")

        with open(test.workspace / f"{name}.ic", "w") as f:
            f.write("BEGIN options\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  strt NETCDF\n")
            f.write("END griddata\n")

        with open(test.workspace / f"{name}.npf", "w") as f:
            f.write("BEGIN options\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  icelltype  NETCDF\n")
            f.write("  k  NETCDF\n")
            f.write("END griddata\n")

        success, buff = flopy.run_model(
            test.targets["mf6"],
            test.workspace / "mfsim.nam",
            model_ws=test.workspace,
            report=True,
        )

        assert success
        test.success = success

    check(idx, test)

    # verify format of generated netcdf file
    with nc.Dataset(test.workspace / f"{name}.nc") as ds:
        assert ds.data_model == "NETCDF4"
        cmpr = ds.variables["head_l1"].filters()
        chnk = ds.variables["head_l1"].chunking()
        assert cmpr["shuffle"]
        assert cmpr["complevel"] == 9
        assert chnk == [1, 3]
        assert ds.variables["projection"].getncattr("wkt").lower() == wkt.lower()

    # Check NetCDF output
    nc_fpth = os.path.join(test.workspace, name + ".nc")
    ds = xu.open_dataset(nc_fpth)
    xds = ds.ugrid.to_dataset()

    hds_fpth = os.path.join(test.workspace, name + ".hds")
    hds = flopy.utils.HeadFile(hds_fpth, precision="double")

    # Compare NetCDF head arrays with binary headfile
    gwf = test.sims[0].gwf[0]
    disv = getattr(gwf, "disv")
    tdis = getattr(test.sims[0], "tdis")
    nper = getattr(tdis, "nper").data
    nlay = getattr(disv, "nlay").data
    pd = getattr(tdis, "perioddata").array
    timestep = 0
    for i in range(nper):
        for j in range(pd[i][1]):
            rec = hds.get_data(kstpkper=(j, i))
            for l in range(nlay):
                assert np.allclose(
                    np.array(rec[l]).flatten(),
                    xds[f"head_l{l + 1}"][timestep, :].fillna(1.00000000e30).data,
                ), f"NetCDF-Headfile comparison failure in timestep {timestep + 1}"
            timestep += 1

    for v in xds.data_vars.keys():
        assert not v.startswith("disv")
        assert not v.startswith("ic")
        assert not v.startswith("npf")

    xds.close()

    if gridded_input == "ascii":
        return

    # Check NetCDF input
    nc_fpth = os.path.join(test.workspace, f"{name}.ugrid.nc")
    ds = xu.open_dataset(nc_fpth)
    xds = ds.ugrid.to_dataset()

    # NetCDF variables, layered variables end with "_l"
    vlist = ["disv_top", "disv_botm_l", "npf_icelltype_l", "npf_k_l", "ic_strt_l"]

    # Compare NetCDF package input arrays with FloPy arrays
    for i, var in enumerate(vlist):
        tokens = var.split("_", 1)
        package_name = tokens[0]
        array_name = tokens[1].split("_")[0]
        package = getattr(gwf, package_name)
        b = getattr(package, array_name).array
        if var.endswith("_l"):
            for l in range(nlay):
                assert np.allclose(
                    np.array(b[l]).flatten(), xds[f"{var}{l + 1}"].data
                ), f"NetCDF input array comparison failure, variable={var}{l + 1}"
        else:
            assert np.allclose(np.array(b).flatten(), xds[var].data), (
                f"NetCDF input array comparison failure, variable={var}"
            )


@pytest.mark.netcdf
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("gridded_input", ["ascii", "netcdf"])
def test_mf6model(idx, name, function_tmpdir, targets, gridded_input):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t, gridded_input),
        check=lambda t: check_output(idx, t, gridded_input),
        cargs=["--mode=validate"] if gridded_input == "netcdf" else None,
        compare=None,
    )
    test.run()
