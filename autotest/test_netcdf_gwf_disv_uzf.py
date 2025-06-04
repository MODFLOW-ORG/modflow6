"""
NetCDF test version of test_gwf_disv_uzf.  The primary aim is to test
that GHBA package NetCDF array input (bhead and cond) gives the same
results as test_gwf_disv_uzf list based (GHB) and array based (GHBA)
ascii input runs.  This test compares heads in the the NetCDF file to
those in the FloPy binary output head file.
"""

# Imports

import os
from pathlib import Path

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
from test_gwf_disv_uzf import cases

xa = pytest.importorskip("xarray")
xu = pytest.importorskip("xugrid")
nc = pytest.importorskip("netCDF4")


def build_models(idx, test):
    from test_gwf_disv_uzf import build_models as build

    sim, mc = build(idx, test)
    gwf = mc.gwf[0]
    gwf.get_package("GHBA_0").export_array_netcdf = True

    name = cases[idx]

    gwf.name_file.nc_mesh2d_filerecord = f"{name}.nc"

    return sim, mc


def check_output(idx, test):
    from test_gwf_disv_uzf import check_output as check

    name = test.name
    ws = Path(test.workspace / "mf6")

    # check outputs of GHB / GHBA ascii input runs
    check(test.workspace, name)
    check(ws, name)

    # verify format of generated netcdf file
    with nc.Dataset(ws / f"{name}.nc") as ds:
        assert ds.data_model == "NETCDF4"

    # re-run the simulation with model netcdf input
    input_fname = f"{name}.nc"
    nc_fname = f"{name}.ugrid.nc"
    os.rename(ws / input_fname, ws / nc_fname)

    fileout_tag = "NETCDF_MESH2D"

    with open(ws / f"{name}.nam", "w") as f:
        f.write("BEGIN options\n")
        f.write("  SAVE_FLOWS\n")
        f.write("  NEWTON\n")
        f.write(f"  {fileout_tag}  FILEOUT  {name}.nc\n")
        f.write(f"  NETCDF  FILEIN {name}.ugrid.nc\n")
        f.write("END options\n\n")
        f.write("BEGIN packages\n")
        f.write(f"  DISV6  {name}.disv  disv\n")
        f.write(f"  IC6  {name}.ic  ic\n")
        f.write(f"  NPF6  {name}.npf  npf\n")
        f.write(f"  STO6  {name}.sto  sto\n")
        f.write(f"  GHBA6  {name}.ghba  ghba_0\n")
        f.write(f"  UZF6  {name}.uzf  uzf_0\n")
        f.write(f"  OC6  {name}.oc  oc\n")
        f.write(f"  OBS6  {name}.obs  head_obs\n")
        f.write("END packages\n")

    with open(ws / f"{name}.ghba", "w") as f:
        f.write("BEGIN options\n")
        f.write("  PRINT_INPUT\n")
        f.write("  PRINT_FLOWS\n")
        f.write("END options\n\n")
        f.write("BEGIN period 1\n")
        f.write("  bhead NETCDF\n")
        f.write("  cond NETCDF\n")
        f.write("END period 1\n")

    success, buff = flopy.run_model(
        test.targets["mf6"],
        ws / "mfsim.nam",
        model_ws=ws,
        report=True,
    )

    assert success
    test.success = success

    # check netcdf input based run
    check(ws, test.name)

    # compare head files for original
    # list based and netcdf input runs
    ext = ["hds"]
    text = ["head"]
    names = [test.name]
    for i, e in enumerate(ext):
        fpth1 = os.path.join(
            test.workspace,
            f"{names[i]}.{e}",
        )
        fpth2 = os.path.join(ws, f"{names[i]}.{e}")
        fout = os.path.join(
            ws,
            f"{names[i]}.{e}.cmp.out",
        )
        success_tst = flopy.utils.compare.compare_heads(
            None,
            None,
            text=f"{text[i]}",
            outfile=fout,
            files1=fpth1,
            files2=fpth2,
            difftol=True,
        )
        msg = f"initial {text[i]} comparison success = {success_tst}"
        if success_tst:
            test.success = True
            print(msg)
        else:
            test.success = False
            assert success_tst, msg

    # now compare heads in head file and
    # netcdf export for netcdf input run
    try:
        # load heads
        fpth = os.path.join(ws, f"{name}.hds")
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        heads = hobj.get_alldata()
    except:
        assert False, f'could not load headfile data from "{fpth}"'

    # open dataset
    nc_fpth = os.path.join(ws, f"{name}.nc")
    ds = xu.open_dataset(nc_fpth)
    xds = ds.ugrid.to_dataset()

    # Compare NetCDF head arrays with binary headfile
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
            for l in range(nlay):
                assert np.allclose(
                    np.array(rec[l]).ravel(),
                    xds[f"head_l{l + 1}"][kstp, :].data,
                ), f"NetCDF-head comparison failure in timestep {kstp + 1}"
            kstp += 1


@pytest.mark.netcdf
@pytest.mark.developmode
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
