"""
Grid array version of test_gwf_csub_db01_nr. The CHD time series (TS6,
linear) is replaced by a CHDG head array given as a full grid time-array
series (TAS6, linear) with the same records at the CHD cell and DNODATA
elsewhere. Heads must match the list based CHD run.

The netcdf variant also reads the DIS, IC, NPF, and STO grid data from a
NetCDF input file exported in validate mode, while CHDG keeps the ASCII
TAS6 file; this combines NetCDF array input with grid time-array series
input.
"""

import shutil

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_csub_db01_nr import build_models as build_list
from test_gwf_csub_db01_nr import c6, cases, hmin, nlay, perlen, strt, tsp0

modes = ["ascii", pytest.param("netcdf", marks=pytest.mark.netcdf)]
tas_name = "chdts"
DNODATA = 3.0e30
chd_cell = tuple(c6[0][0:3])

# CHD time series records of test_gwf_csub_db01_nr
records = [
    (0.0, strt),
    (tsp0, strt),
    (365.25 + tsp0, hmin),
    (np.array(perlen).sum(), strt),
]

# packages with grid data read from NetCDF in the netcdf variant
nc_packages = ["dis", "ic", "npf", "sto"]


class Workspace:
    """Stand-in for the test object used by build_list."""

    def __init__(self, workspace):
        self.workspace = workspace


def build_models(idx, test):
    sim, _ = build_list(idx, test)
    gwf = sim.gwf[0]
    gwf.remove_package(gwf.get_package("chd"))
    shape = (nlay, 1, 2)
    tas = {}
    for t, v in records:
        a = np.full(shape, DNODATA)
        a[chd_cell] = v
        tas[t] = a
    chdg = flopy.mf6.ModflowGwfchdg(
        gwf, maxbound=1, head={0: f"TIMEARRAYSERIES {tas_name}"}
    )
    chdg.tas.initialize(
        filename=f"{cases[idx]}.chdg.tas",
        tas_array=tas,
        time_series_namerecord=tas_name,
        interpolation_methodrecord="linear",
    )
    return sim, None


def griddata_to_netcdf(text):
    """Replace each GRIDDATA array in a package file with NETCDF input."""
    out, ingrid, skip = [], False, False
    for line in text.splitlines(keepends=True):
        words = line.split()
        if words[:2] == ["BEGIN", "griddata"]:
            ingrid = True
        elif words[:2] == ["END", "griddata"]:
            ingrid, skip = False, False
        elif ingrid and words:
            indent = len(line) - len(line.lstrip())
            if indent <= 2:
                out.append(f"  {words[0]}  NETCDF\n")
                skip = True
                continue
            if skip:
                continue
        out.append(line)
    return "".join(out)


def to_netcdf(ws, name, exe):
    """Export grid data in validate mode and rerun from NetCDF input."""
    ncfile = f"{name}.structured.nc"
    nam = (ws / f"{name}.nam").read_text()
    files = {p: (ws / f"{name}.{p}").read_text() for p in nc_packages}
    (ws / f"{name}.nam").write_text(
        nam.replace(
            "END options", f"  NETCDF_STRUCTURED  FILEOUT  {ncfile}\nEND options", 1
        )
    )
    for p, text in files.items():
        (ws / f"{name}.{p}").write_text(
            text.replace("END options", "  EXPORT_ARRAY_NETCDF\nEND options", 1)
        )
    success, _ = flopy.run_model(
        exe, ws / "mfsim.nam", model_ws=ws, cargs=["--mode=validate"], silent=True
    )
    assert success, "validate mode NetCDF export failed"

    (ws / f"{name}.nam").write_text(
        nam.replace("END options", f"  NETCDF  FILEIN  {ncfile}\nEND options", 1)
    )
    for p, text in files.items():
        converted = griddata_to_netcdf(text)
        assert "NETCDF" in converted, f"{p} has no grid data"
        (ws / f"{name}.{p}").write_text(converted)
    success, _ = flopy.run_model(exe, ws / "mfsim.nam", model_ws=ws, silent=True)
    assert success, "NetCDF input run failed"


def check_output(idx, test, mode):
    name = cases[idx]
    exe = test.targets["mf6"]
    ws = test.workspace
    if mode == "netcdf":
        ws = test.workspace / "netcdf"
        shutil.copytree(test.workspace, ws, ignore=shutil.ignore_patterns("list"))
        to_netcdf(ws, name, exe)

    # reference: list based CHD with a TS6 time series
    ref = test.workspace / "list"
    sim, _ = build_list(idx, Workspace(ref))
    sim.write_simulation(silent=True)
    success, _ = flopy.run_model(exe, ref / "mfsim.nam", model_ws=ref, silent=True)
    assert success, "reference list CHD simulation failed"

    h = flopy.utils.HeadFile(ws / f"{name}.hds").get_alldata()
    href = flopy.utils.HeadFile(ref / f"{name}.hds").get_alldata()
    assert h.shape == href.shape
    assert np.allclose(h, href, atol=1e-8), (
        f"CHDG TAS heads differ from CHD TS: max diff {np.abs(h - href).max()}"
    )


@pytest.mark.parametrize("mode", modes)
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, mode, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, mode),
        targets=targets,
    )
    test.run()
