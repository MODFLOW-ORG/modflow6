"""
NetCDF export test version of test_gwf_lak_wetlakbedarea02.  This test compares
the temperature and input arrays in the the NetCDF file to those
in the FloPy binary output head file and package data objects.
"""

# Imports

import os

import numpy as np
import pytest

try:
    import flopy
except:
    msg = "Error. FloPy package is not available.\n"
    msg += "Try installing using the following command:\n"
    msg += " pip install flopy"
    raise Exception(msg)

from framework import TestFramework
from test_gwf_lak_wetlakbedarea02 import cases

xa = pytest.importorskip("xarray")
xu = pytest.importorskip("xugrid")
nc = pytest.importorskip("netCDF4")


def build_models(idx, test, export, gridded_input):
    from test_gwf_lak_wetlakbedarea02 import build_models as build

    sim, dummy = build(idx, test)
    sim.tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    gwf = sim.gwf[0]
    gwf.dis.export_array_netcdf = True
    gwf.ic.export_array_netcdf = True
    gwf.npf.export_array_netcdf = True
    gwf.rch.export_array_netcdf = True
    gwf.evt.export_array_netcdf = True

    name = "gwf-" + cases[idx]

    fname = f"{name}.{export}.nc" if gridded_input == "netcdf" else f"{name}.nc"
    if export == "ugrid":
        gwf.name_file.nc_mesh2d_filerecord = fname
    elif export == "structured":
        gwf.name_file.nc_structured_filerecord = fname

    # netcdf config
    ncf = flopy.mf6.ModflowUtlncf(
        gwf.dis,
        filename=f"{name}.dis.ncf",
    )

    return sim, dummy


def check_output(idx, test, export, gridded_input):
    from test_gwf_lak_wetlakbedarea02 import check_output as check

    name = "gwf-" + cases[idx]

    # verify format of generated netcdf file
    fname = f"{name}.{export}.nc" if gridded_input == "netcdf" else f"{name}.nc"
    with nc.Dataset(test.workspace / fname) as ds:
        assert ds.data_model == "NETCDF4"

    if gridded_input == "netcdf":
        # re-run the simulation with model netcdf input

        if export == "ugrid":
            fileout_tag = "NETCDF_MESH2D"
        elif export == "structured":
            fileout_tag = "NETCDF_STRUCTURED"

        with open(test.workspace / f"{name}.nam", "w") as f:
            f.write("BEGIN options\n")
            f.write("  SAVE_FLOWS\n")
            f.write("  NEWTON\n")
            f.write(f"  {fileout_tag}  FILEOUT  {name}.nc\n")
            f.write(f"  NETCDF  FILEIN {name}.{export}.nc\n")
            f.write("END options\n\n")
            f.write("BEGIN packages\n")
            f.write(f"  DIS6  {name}.dis  dis\n")
            f.write(f"  NPF6  {name}.npf  npf\n")
            f.write(f"  STO6  {name}.sto  sto\n")
            f.write(f"  IC6  {name}.ic  ic\n")
            f.write(f"  CHD6  {name}.chd  chd_0\n")
            f.write(f"  RCH6  {name}.rcha  rcha_0\n")
            f.write(f"  EVT6  {name}.evta  evta_0\n")
            f.write(f"  LAK6  {name}.lak  lak-1\n")
            f.write(f"  OC6  {name}.oc  oc\n")
            f.write("END packages\n")

        with open(test.workspace / f"{name}.dis", "w") as f:
            f.write("BEGIN options\n")
            f.write("  LENGTH_UNITS feet\n")
            f.write(f"  NCF6  FILEIN  {name}.dis.ncf\n")
            f.write("END options\n\n")
            f.write("BEGIN dimensions\n")
            f.write("  NLAY  6\n")
            f.write("  NROW  17\n")
            f.write("  NCOL  17\n")
            f.write("END dimensions\n\n")
            f.write("BEGIN griddata\n")
            f.write("  delr NETCDF\n")
            f.write("  delc NETCDF\n")
            f.write("  top NETCDF\n")
            f.write("  botm NETCDF\n")
            f.write("  idomain NETCDF\n")
            f.write("END griddata\n\n")

        with open(test.workspace / f"{name}.ic", "w") as f:
            f.write("BEGIN options\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  strt NETCDF\n")
            f.write("END griddata\n")

        with open(test.workspace / f"{name}.npf", "w") as f:
            f.write("BEGIN options\n")
            f.write("  SAVE_SPECIFIC_DISCHARGE\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  icelltype  NETCDF\n")
            f.write("  k  NETCDF\n")
            f.write("  k33  NETCDF\n")
            f.write("END griddata\n")

        with open(test.workspace / f"{name}.rcha", "w") as f:
            f.write("BEGIN options\n")
            f.write("  READASARRAYS\n")
            f.write("END options\n\n")
            f.write("BEGIN period 1\n")
            f.write("  recharge NETCDF\n")
            f.write("END period 1\n")

        with open(test.workspace / f"{name}.evta", "w") as f:
            f.write("BEGIN options\n")
            f.write("  READASARRAYS\n")
            f.write("END options\n\n")
            f.write("BEGIN period 1\n")
            f.write("  surface NETCDF\n")
            f.write("  rate NETCDF\n")
            f.write("  depth NETCDF\n")
            f.write("END period 1\n")

        success, buff = flopy.run_model(
            test.targets["mf6"],
            test.workspace / "mfsim.nam",
            model_ws=test.workspace,
            report=True,
        )

        assert success
        test.success = success

    check(idx, test)

    try:
        # load heads
        fname = name + ".hds"
        fpth = os.path.join(test.workspace, fname)
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        heads = hobj.get_alldata()
    except:
        assert False, f'could not load headfile data from "{fpth}"'

    # Check NetCDF output
    nc_fpth = os.path.join(test.workspace, f"{name}.nc")
    if export == "ugrid":
        ds = xu.open_dataset(nc_fpth)
        xds = ds.ugrid.to_dataset()
    elif export == "structured":
        xds = xa.open_dataset(nc_fpth)

    # Compare NetCDF head arrays with binary headfile temperatures
    gwf = test.sims[0].gwf[0]
    dis = getattr(gwf, "dis")
    tdis = getattr(test.sims[0], "tdis")
    nper = getattr(tdis, "nper").data
    nlay = getattr(dis, "nlay").data
    pd = getattr(tdis, "perioddata").array
    kstp = 0
    for i in range(nper):
        for j in range(int(pd[i][1])):
            rec = hobj.get_data(kstpkper=(j, i))
            if export == "ugrid":
                for l in range(nlay):
                    assert np.allclose(
                        np.array(rec[l]).flatten(),
                        # xds[f"head_l{l+1}"][timestep, :].data,
                        xds[f"head_l{l + 1}"][kstp, :].fillna(1.00000000e30).data,
                    ), f"NetCDF-temperature comparison failure in timestep {kstp + 1}"
                kstp += 1
            elif export == "structured":
                assert np.allclose(
                    # np.array(rec).flatten(),
                    np.array(rec),
                    xds["head"][kstp, :].fillna(1.00000000e30).data,
                ), f"NetCDF-head comparison failure in timestep {kstp + 1}"
                kstp += 1

    for v in xds.data_vars.keys():
        assert not v.startswith("dis")
        assert not v.startswith("ic")
        assert not v.startswith("npf")

    xds.close()

    if gridded_input == "ascii":
        return

    # Check NetCDF input
    nc_fpth = os.path.join(test.workspace, f"{name}.{export}.nc")
    if export == "ugrid":
        ds = xu.open_dataset(nc_fpth)
        xds = ds.ugrid.to_dataset()
    elif export == "structured":
        xds = xa.open_dataset(nc_fpth)

    vlist = [
        "dis_delr",
        "dis_delc",
        "dis_top",
        "dis_botm_l",
        "dis_idomain_l",
        "ic_strt_l",
        "npf_icelltype_l",
        "npf_k_l",
    ]

    # Compare NetCDF package input arrays with FloPy arrays
    gwf = test.sims[0].gwf[0]
    for i, var in enumerate(vlist):
        tokens = var.split("_", 1)
        package_name = tokens[0]
        array_name = tokens[1].split("_")[0]
        package = getattr(gwf, package_name)
        b = getattr(package, array_name).array
        if export == "ugrid":
            if var.endswith("_l"):
                for l in range(nlay):
                    assert np.allclose(
                        np.array(b[l]).flatten(), xds[f"{var}{l + 1}"].data
                    ), f"NetCDF input array comparison failure, variable={var}{l + 1}"
            else:
                assert np.allclose(np.array(b).flatten(), xds[var].data), (
                    f"NetCDF input array comparison failure, variable={var}"
                )
        elif export == "structured":
            var = var.replace("_l", "")
            assert np.allclose(
                # np.array(b).flatten(), xds[var].data
                np.array(b),
                xds[var].data,
            ), f"NetCDF input array comparison failure, variable={var}"


@pytest.mark.netcdf
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
@pytest.mark.parametrize("export", ["ugrid", "structured"])
@pytest.mark.parametrize("gridded_input", ["ascii", "netcdf"])
def test_mf6model(idx, name, function_tmpdir, targets, export, gridded_input):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t, export, gridded_input),
        check=lambda t: check_output(idx, t, export, gridded_input),
        targets=targets,
        cargs=["--mode=validate"] if gridded_input == "netcdf" else None,
    )
    test.run()
