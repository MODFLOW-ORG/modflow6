"""
Framework helpers for the re-run comparisons (flopy4 re-writes and MF6's
own NetCDF export), tested without flopy4 or a MODFLOW 6 run.
"""

from pathlib import Path

import pytest
from flopy4_input import check_netcdf_input
from rerun import copy_for_rerun

nc = pytest.importorskip("netCDF4")


def _write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_copy_for_rerun(tmp_path):
    """Outputs and other comparisons are left out of the copy; a file an
    FMI package reads is kept, even with an output file's extension."""
    src = tmp_path / "src"
    _write(src / "mfsim.nam", "BEGIN models\n  gwt6 gwt.nam gwt\nEND models\n")
    _write(
        src / "gwt.nam",
        "BEGIN packages\n  DIS6 gwt.dis dis\n  FMI6 gwt.fmi fmi\nEND packages\n",
    )
    _write(src / "gwt.dis", "BEGIN dimensions\n  NCPL 1\nEND dimensions\n")
    _write(
        src / "gwt.fmi",
        "BEGIN packagedata\n  GWFBUDGET FILEIN flow.bud\nEND packagedata\n",
    )
    for name in ("flow.bud", "gwt.lst", "gwt.ucn", "mf6_regression/gwt.nam"):
        _write(src / name, "")

    dst = tmp_path / "dst"
    assert copy_for_rerun(src, dst, ["."]) == [Path(".")]
    copied = {p.relative_to(dst).as_posix() for p in dst.rglob("*") if p.is_file()}
    assert copied == {"mfsim.nam", "gwt.nam", "gwt.dis", "gwt.fmi", "flow.bud"}


def _workspace(ws: Path, npf_k: str, variables: dict[str, str]):
    """A GWF model reading DIS and NPF arrays from `m.nc`, which holds
    `variables` (name to modflow_input)."""
    _write(ws / "mfsim.nam", "BEGIN models\n  gwf6 m.nam m\nEND models\n")
    _write(
        ws / "m.nam",
        "BEGIN options\n  NETCDF FILEIN m.nc\nEND options\n"
        "BEGIN packages\n  DIS6 m.dis dis\n  NPF6 m.npf npf\nEND packages\n",
    )
    _write(
        ws / "m.dis",
        "BEGIN griddata\n  delr NETCDF\n  delc NETCDF\n  top NETCDF\n"
        "  botm NETCDF\nEND griddata\n",
    )
    _write(ws / "m.npf", f"BEGIN griddata\n  icelltype NETCDF\n{npf_k}END griddata\n")
    with nc.Dataset(ws / "m.nc", "w", format="NETCDF4") as ds:
        ds.modflow_model = "gwf: m"
        ds.createDimension("x", 1)
        for name, attr in variables.items():
            ds.createVariable(name, "f8", ("x",)).modflow_input = attr


ARRAYS = {
    f"{p}_{t}": f"m/{p}/{t}"
    for p, t in [
        ("dis", "delr"),
        ("dis", "delc"),
        ("dis", "top"),
        ("dis", "botm"),
        ("npf", "icelltype"),
        ("npf", "k"),
    ]
}


def test_check_netcdf_input(tmp_path):
    """Every array MF6 can read from NetCDF is, from a variable in the
    model's file."""
    _workspace(tmp_path, "  k NETCDF\n", ARRAYS)
    check_netcdf_input(tmp_path)


def test_check_netcdf_input_ascii_array(tmp_path):
    """An array MF6 could read from NetCDF, given as text, is an error."""
    _workspace(tmp_path, "  k\n    CONSTANT 1.0\n", ARRAYS)
    with pytest.raises(AssertionError, match="k not read from NetCDF"):
        check_netcdf_input(tmp_path)


def test_check_netcdf_input_missing_variable(tmp_path):
    """An array read with NETCDF needs a variable in the model's file."""
    arrays = {k: v for k, v in ARRAYS.items() if k != "npf_k"}
    _workspace(tmp_path, "  k NETCDF\n", arrays)
    with pytest.raises(AssertionError, match="k has no variable"):
        check_netcdf_input(tmp_path)
