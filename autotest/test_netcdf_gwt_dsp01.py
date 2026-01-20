"""
NetCDF export test version of test_gwt_dsp01.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwt_dsp01 import cases, xt3d

xa = pytest.importorskip("xarray")
xu = pytest.importorskip("xugrid")
nc = pytest.importorskip("netCDF4")


def build_models(idx, test, export, gridded_input):
    from test_gwt_dsp01 import build_models as build

    sim, dummy = build(idx, test)
    sim.tdis.start_date_time = "2041-01-01T00:00:00-05:00"
    gwt = sim.gwt[0]
    gwt.dis.export_array_netcdf = True
    gwt.ic.export_array_netcdf = True
    gwt.dsp.export_array_netcdf = True

    # output control
    gwtname = "gwt_" + cases[idx]

    fname = f"{gwtname}.{export}.nc" if gridded_input == "netcdf" else f"{gwtname}.nc"
    if export == "ugrid":
        gwt.name_file.nc_mesh2d_filerecord = fname
        ncf = flopy.mf6.ModflowUtlncf(
            gwt.dis,
            deflate=3,
            shuffle=False,
            chunk_time=1,
            chunk_face=5,
            filename=f"{gwtname}.dis.ncf",
        )
    elif export == "structured":
        gwt.name_file.nc_structured_filerecord = fname
        ncf = flopy.mf6.ModflowUtlncf(
            gwt.dis,
            deflate=3,
            shuffle=False,
            chunk_time=1,
            chunk_z=1,
            chunk_y=1,
            chunk_x=20,
            filename=f"{gwtname}.dis.ncf",
        )

    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    )

    return sim, dummy


def check_output(idx, test, export, gridded_input):
    from test_gwt_dsp01 import check_output as check

    name = cases[idx]
    gwtname = "gwt_" + name

    if gridded_input == "netcdf":
        # re-run the simulation with model netcdf input
        if export == "ugrid":
            fileout_tag = "NETCDF_MESH2D"
        elif export == "structured":
            fileout_tag = "NETCDF_STRUCTURED"

        with open(test.workspace / f"{gwtname}.nam", "w") as f:
            f.write("BEGIN options\n")
            f.write("  SAVE_FLOWS\n")
            f.write(f"  {fileout_tag}  FILEOUT  {gwtname}.nc\n")
            f.write(f"  NETCDF  FILEIN {gwtname}.{export}.nc\n")
            f.write("END options\n\n")
            f.write("BEGIN packages\n")
            f.write(f"  OBS6  {gwtname}.obs  gwt_obs\n")
            f.write(f"  SSM6  {gwtname}.ssm  ssm\n")
            f.write(f"  MST6  {gwtname}.mst  mst\n")
            f.write(f"  CNC6  {gwtname}.cnc  cnc-1\n")
            f.write(f"  DSP6  {gwtname}.dsp  dsp\n")
            f.write(f"  ADV6  {gwtname}.adv  adv\n")
            f.write(f"  IC6  {gwtname}.ic  ic\n")
            f.write(f"  DIS6  {gwtname}.dis  dis\n")
            f.write(f"  OC6  {gwtname}.oc  oc\n")
            f.write("END packages\n")

        with open(test.workspace / f"{gwtname}.dis", "w") as f:
            f.write("BEGIN options\n")
            f.write(f"  NCF6  FILEIN  {gwtname}.dis.ncf\n")
            f.write("END options\n\n")
            f.write("BEGIN dimensions\n")
            f.write("  NLAY  1\n")
            f.write("  NROW  1\n")
            f.write("  NCOL  100\n")
            f.write("END dimensions\n\n")
            f.write("BEGIN griddata\n")
            f.write("  delr NETCDF\n")
            f.write("  delc NETCDF\n")
            f.write("  top NETCDF\n")
            f.write("  botm NETCDF\n")
            f.write("  idomain NETCDF\n")
            f.write("END griddata\n\n")

        with open(test.workspace / f"{gwtname}.ic", "w") as f:
            f.write("BEGIN options\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  strt NETCDF\n")
            f.write("END griddata\n")

        with open(test.workspace / f"{gwtname}.dsp", "w") as f:
            f.write("BEGIN options\n")
            if not xt3d[idx]:
                f.write("  XT3D_OFF\n")
            f.write("END options\n\n")
            f.write("BEGIN griddata\n")
            f.write("  diffc  NETCDF\n")
            f.write("  alh  NETCDF\n")
            f.write("  alv  NETCDF\n")
            f.write("  ath1  NETCDF\n")
            f.write("  atv  NETCDF\n")
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
    with nc.Dataset(test.workspace / f"{gwtname}.nc") as ds:
        assert ds.data_model == "NETCDF4"
        if export == "structured":
            cmpr = ds.variables["concentration"].filters()
            chnk = ds.variables["concentration"].chunking()
            assert chnk == [1, 1, 1, 20]
        elif export == "ugrid":
            cmpr = ds.variables["concentration_l1"].filters()
            chnk = ds.variables["concentration_l1"].chunking()
            assert chnk == [1, 5]
        assert not cmpr["shuffle"]
        assert cmpr["complevel"] == 3

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # netcdf
    nc_fpth = os.path.join(test.workspace, f"{gwtname}.nc")
    if export == "ugrid":
        ds = xu.open_dataset(nc_fpth)
        xds = ds.ugrid.to_dataset()
    elif export == "structured":
        xds = xa.open_dataset(nc_fpth)

    # Compare NetCDF head arrays with binary headfile concentrations
    gwt = test.sims[0].gwt[0]
    dis = getattr(gwt, "dis")
    tdis = getattr(test.sims[0], "tdis")
    nper = getattr(tdis, "nper").data
    nlay = getattr(dis, "nlay").data
    pd = getattr(tdis, "perioddata").array
    timestep = 0
    for i in range(nper):
        for j in range(pd[i][1]):
            rec = cobj.get_data(kstpkper=(j, i))
            if export == "ugrid":
                for l in range(nlay):
                    assert np.allclose(
                        np.array(rec[l]).flatten(),
                        xds[f"concentration_l{l + 1}"][timestep, :].data,
                    ), (
                        "NetCDF-concentration comparison failure "
                        f"in timestep {timestep + 1}"
                    )
                timestep += 1
            elif export == "structured":
                assert np.allclose(
                    # np.array(rec).flatten(),
                    np.array(rec),
                    xds["concentration"][timestep, :].data,
                ), f"NetCDF-concentration comparison failure in timestep {timestep + 1}"
                timestep += 1

    for v in xds.data_vars.keys():
        assert not v.startswith("dis")
        assert not v.startswith("ic")
        assert not v.startswith("dsp")

    xds.close()

    if gridded_input == "ascii":
        return

    # Check NetCDF input
    nc_fpth = os.path.join(test.workspace, f"{gwtname}.{export}.nc")
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
        "dsp_diffc_l",
        "dsp_alh_l",
        "dsp_alv_l",
        "dsp_ath1_l",
        "dsp_atv_l",
    ]

    # Compare NetCDF package input arrays with FloPy arrays
    gwt = test.sims[0].gwt[0]
    for i, var in enumerate(vlist):
        tokens = var.split("_", 1)
        package_name = tokens[0]
        array_name = tokens[1].split("_")[0]
        package = getattr(gwt, package_name)
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
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
@pytest.mark.parametrize("gridded_input", ["ascii", "netcdf"])
def test_mf6model(idx, name, function_tmpdir, targets, export, gridded_input):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t, export, gridded_input),
        check=lambda t: check_output(idx, t, export, gridded_input),
        cargs=["--mode=validate"] if gridded_input == "netcdf" else None,
    )
    test.run()
