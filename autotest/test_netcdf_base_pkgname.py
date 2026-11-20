"""
NetCDF export names a package that a model can have only one of (DIS, NPF,
IC) by its package type, also when the model name file gives it a name,
so the exported file can be read back as input.
"""

import shutil

import flopy
import numpy as np
import pytest
from framework import TestFramework

nc = pytest.importorskip("netCDF4")

cases = ["nc_pkgname"]
NLAY, NROW, NCOL = 2, 3, 4


def build_models(idx, test, export):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        top=10.0,
        botm=[5.0, 0.0],
        pname="DIS-1",
        export_array_netcdf=True,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        k=np.arange(1.0, NLAY * NROW * NCOL + 1).reshape(NLAY, NROW, NCOL),
        pname="NPF-1",
        export_array_netcdf=True,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0, pname="IC-1", export_array_netcdf=True)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 10.0], [(0, NROW - 1, NCOL - 1), 6.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    if export == "ugrid":
        gwf.name_file.nc_mesh2d_filerecord = f"{name}.export.nc"
    else:
        gwf.name_file.nc_structured_filerecord = f"{name}.export.nc"
    return sim, None


def read_netcdf_input(ws, name):
    """Make the model read its GRIDDATA arrays from the exported file."""
    nam = (ws / f"{name}.nam").read_text().splitlines()
    nam = [
        f"  NETCDF FILEIN {name}.export.nc" if "FILEOUT" in line.upper() else line
        for line in nam
    ]
    (ws / f"{name}.nam").write_text("\n".join(nam) + "\n")
    griddata = {
        "dis": f"BEGIN dimensions\n  NLAY {NLAY}\n  NROW {NROW}\n  NCOL {NCOL}\n"
        "END dimensions\n\nBEGIN griddata\n  delr NETCDF\n  delc NETCDF\n"
        "  top NETCDF\n  botm NETCDF\nEND griddata\n",
        "npf": "BEGIN griddata\n  icelltype NETCDF\n  k NETCDF\nEND griddata\n",
        "ic": "BEGIN griddata\n  strt NETCDF\nEND griddata\n",
    }
    for ext, text in griddata.items():
        (ws / f"{name}.{ext}").write_text(f"BEGIN options\nEND options\n\n{text}")


def check_output(idx, test, export):
    name = cases[idx]
    ws = test.workspace
    suffix = "_l1" if export == "ugrid" else ""

    # each array is named by package type, not by the name file's name
    with nc.Dataset(ws / f"{name}.export.nc") as ds:
        for pkg, tag in (("dis", "botm"), ("npf", "k"), ("ic", "strt")):
            var = ds.variables[f"{pkg}_{tag}{suffix}"]
            assert var.modflow_input == f"{name}/{pkg}/{tag}", var.modflow_input
        assert not any("-1_" in v for v in ds.variables)

    # the ASCII run, then the run reading the exported arrays as input
    ascii_ws = ws / "ascii"
    shutil.copytree(ws, ascii_ws)
    success, _ = flopy.run_model(
        test.targets["mf6"], ascii_ws / "mfsim.nam", model_ws=ascii_ws
    )
    assert success, "ASCII run failed"

    read_netcdf_input(ws, name)
    success, buff = flopy.run_model(test.targets["mf6"], ws / "mfsim.nam", model_ws=ws)
    assert success, "run with NetCDF input failed:\n" + "\n".join(buff[-20:])

    heads = [flopy.utils.HeadFile(w / f"{name}.hds").get_data() for w in (ascii_ws, ws)]
    assert np.allclose(*heads)


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
def test_mf6model(idx, name, function_tmpdir, targets, export):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t, export),
        check=lambda t: check_output(idx, t, export),
        compare=None,
        cargs=["--mode=validate"],
    )
    test.run()
