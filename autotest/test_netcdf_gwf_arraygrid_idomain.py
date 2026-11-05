"""
Grid array (WELG, CHDG) input with values at cells with IDOMAIN less than
one, exported to NetCDF in validate mode and read back with NETCDF period
input, with a value written at the inactive cell in the exported file.
Results must match the ASCII input run.
"""

import flopy
import pytest
from framework import TestFramework
from test_gwf_arraygrid_idomain import build_models, cases, cell_inactive, ncol
from test_gwf_arraygrid_idomain import check_output as check

pytest.importorskip("xarray")
pytest.importorskip("xugrid")
pytest.importorskip("netCDF4")

# grid array package and param per case
pkg_param = [("welg", "q"), ("chdg", "head"), ("welg", "q")]


def set_inactive_value(ncfile, param, value=99.0):
    """Write a value at the inactive cell in every period."""
    import netCDF4

    k, i, j = cell_inactive
    with netCDF4.Dataset(ncfile, "r+") as ds:
        for v in ds.variables.values():
            tag = v.modflow_input if "modflow_input" in v.ncattrs() else ""
            if not tag.lower().endswith(f"/{param}"):
                continue
            if v.ndim == 4:
                v[:, k, i, j] = value
            elif v.layer == k + 1:
                v[:, i * ncol + j] = value


def check_output(idx, test, export):
    name = cases[idx]
    ws = test.workspace
    sim = test.sims[0]
    gwf = sim.get_model(name)
    pkgtype, param = pkg_param[idx]
    pkg = gwf.get_package(pkgtype)
    nper = sim.tdis.nper.get_data()

    # check ascii input run
    check(idx, test)

    # re-run in validate mode to generate netcdf input
    pkg.export_array_netcdf = True
    if export == "ugrid":
        gwf.name_file.nc_mesh2d_filerecord = f"{name}.{export}.nc"
    elif export == "structured":
        gwf.name_file.nc_structured_filerecord = f"{name}.{export}.nc"
    sim.write_simulation()
    success, _ = flopy.run_model(
        test.targets["mf6"],
        ws / "mfsim.nam",
        model_ws=ws,
        report=True,
        cargs=["--mode=validate"],
    )
    assert success

    # use netcdf input file
    gwf.name_file.nc_mesh2d_filerecord = None
    gwf.name_file.nc_structured_filerecord = None
    gwf.name_file.nc_filerecord = f"{name}.{export}.nc"
    gwf.name_file.write()
    set_inactive_value(ws / f"{name}.{export}.nc", param)

    # rewrite grid array package with netcdf period input
    with open(ws / pkg.filename, "w") as f:
        f.write("BEGIN options\n")
        f.write("  READARRAYGRID\n")
        f.write("END options\n\n")
        f.write("BEGIN dimensions\n")
        f.write("  MAXBOUND  3\n")
        f.write("END dimensions\n\n")
        for kper in range(1, nper + 1):
            f.write(f"BEGIN period {kper}\n")
            f.write(f"  {param} NETCDF\n")
            f.write(f"END period {kper}\n\n")

    success, _ = flopy.run_model(
        test.targets["mf6"],
        ws / "mfsim.nam",
        model_ws=ws,
        report=True,
    )
    assert success

    # check netcdf input run
    check(idx, test)


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases[: len(pkg_param)]))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
def test_mf6model(idx, name, function_tmpdir, targets, export):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, export),
        targets=targets,
    )
    test.run()
