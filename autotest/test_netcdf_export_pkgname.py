"""
NetCDF export names each package's arrays the way MODFLOW 6 reads them back
as input: a package a model can have only one of (DIS, NPF, IC) by its
package type, whatever name the model name file gives it, and a package a
model can have several of (RCHA, EVTA, WELG, GHBG) by its name file name,
or MODFLOW 6's default `<type>-<n>` when the name file gives none. Each
case is exported in validate mode, read back as input with NETCDF in the
GRIDDATA and PERIOD blocks, and run, and its heads compared to the run
with text input.
"""

import re
import shutil

import flopy
import numpy as np
import pytest
from framework import TestFramework

nc = pytest.importorskip("netCDF4")

cases = ["nc_pkgname"]
NLAY, NROW, NCOL = 2, 3, 4
NO = 3.0e30

# naming mode: the package names expected in the exported file. In "named"
# mode DIS and NPF are given names other than their types, IC its own type,
# RCHA and GHBG custom names, and EVTA and WELG flopy's names. In "unnamed"
# mode the name file gives no names.
EXPECTED = {
    "named": {
        "dis": "dis",
        "npf": "npf",
        "ic": "ic",
        "rch": "rch-x",
        "evt": "evta_0",
        "wel": ["welg_0", "welg_1"],
        "ghb": "ghb-x",
    },
    "unnamed": {
        "dis": "dis",
        "npf": "npf",
        "ic": "ic",
        "rch": "rch-1",
        "evt": "evt-1",
        "wel": ["wel-1", "wel-2"],
        "ghb": "ghb-1",
    },
}

# arrays given in each package file, read back with NETCDF
ARRAYS = {
    "dis": {"delr", "delc", "top", "botm"},
    "npf": {"icelltype", "k"},
    "ic": {"strt"},
    "rcha": {"recharge"},
    "evta": {"surface", "rate", "depth"},
    "welg": {"q"},
    "ghbg": {"bhead", "cond", "conc"},
}


def grid_array(cells):
    """A grid array holding values at (layer, row, column) cells only."""
    a = np.full((NLAY, NROW, NCOL), NO)
    for cell, value in cells.items():
        a[cell] = value
    return a


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    export = {"export_array_netcdf": True}
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        top=10.0,
        botm=[5.0, 0.0],
        pname="DIS-X",
        **export,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        k=np.arange(1.0, NLAY * NROW * NCOL + 1).reshape(NLAY, NROW, NCOL),
        pname="NPF-1",
        **export,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0, pname="IC", **export)
    flopy.mf6.ModflowGwfrcha(
        gwf, recharge={0: 1.0e-3, 1: 2.0e-3}, pname="RCH-X", **export
    )
    flopy.mf6.ModflowGwfevta(
        gwf, surface={0: 9.0}, rate={0: 1.0e-4}, depth={0: 2.0}, **export
    )
    flopy.mf6.ModflowGwfwelg(
        gwf, maxbound=1, q={0: grid_array({(1, 1, 1): -0.5})}, **export
    )
    flopy.mf6.ModflowGwfwelg(
        gwf, maxbound=1, q={1: grid_array({(1, 1, 2): -0.25})}, **export
    )
    bnd = {(0, 0, 0): 10.0, (0, NROW - 1, NCOL - 1): 6.0}
    flopy.mf6.ModflowGwfghbg(
        gwf,
        maxbound=2,
        auxiliary=["conc"],
        bhead={0: grid_array(bnd)},
        cond={0: grid_array(dict.fromkeys(bnd, 100.0))},
        aux={0: np.array([grid_array(dict.fromkeys(bnd, 1.0))])},
        pname="GHB-X",
        **export,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    return sim, None


def strip_names(path):
    """Remove the package names from a name file's PACKAGES block."""
    lines, inside = [], False
    for line in path.read_text().splitlines():
        tokens = line.split()
        upper = [t.upper() for t in tokens[:2]]
        if upper in (["BEGIN", "PACKAGES"], ["END", "PACKAGES"]):
            inside = upper[0] == "BEGIN"
        elif inside and len(tokens) > 2:
            line = f"  {tokens[0]} {tokens[1]}"
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")


def set_netcdf_file(path, option):
    """Add a NetCDF file option to a name file's OPTIONS block, replacing
    any NetCDF file option already there."""
    lines = [
        line
        for line in path.read_text().splitlines()
        if not re.match(r"\s*NETCDF", line, re.IGNORECASE)
    ]
    i = next(i for i, l in enumerate(lines) if l.strip().upper() == "END OPTIONS")
    lines.insert(i, f"  {option}")
    path.write_text("\n".join(lines) + "\n")


def read_netcdf_input(path, names):
    """Replace each array in `names`, given with values in the GRIDDATA
    or a PERIOD block, with its NETCDF entry."""
    lines, skip = [], False
    for line in path.read_text().splitlines():
        tokens = line.split()
        if tokens and tokens[0].upper() in ("BEGIN", "END"):
            skip = False
        elif tokens and tokens[0].lower() in names:
            skip = True
            lines.append(f"  {tokens[0]} NETCDF")
            continue
        if not skip:
            lines.append(line)
    path.write_text("\n".join(lines) + "\n")


def check_output(idx, test, export, naming):
    name = cases[idx]
    ws = test.workspace / "export"
    shutil.copytree(test.workspace, ws)
    namefile = ws / f"{name}.nam"
    if naming == "unnamed":
        strip_names(namefile)
    record = "NETCDF_MESH2D" if export == "ugrid" else "NETCDF_STRUCTURED"
    set_netcdf_file(namefile, f"{record} FILEOUT {name}.export.nc")
    success, buff = flopy.run_model(
        test.targets["mf6"], ws / "mfsim.nam", model_ws=ws, cargs=["--mode=validate"]
    )
    assert success, "validate run failed:\n" + "\n".join(buff[-20:])

    # each array is named for its package as MF6 reads it back
    expected = EXPECTED[naming]
    packages = {
        p for v in expected.values() for p in (v if isinstance(v, list) else [v])
    }
    with nc.Dataset(ws / f"{name}.export.nc") as ds:
        found = set()
        for var in ds.variables.values():
            if "modflow_input" not in var.ncattrs():
                continue
            model, package, tag = var.modflow_input.split("/")
            assert model == name, var.modflow_input
            assert var.name.startswith(f"{package}_"), (var.name, var.modflow_input)
            found.add(package)
            if tag == "aux":
                assert var.modflow_iaux == 1
                assert var.name.startswith(f"{package}_conc"), var.name
        assert found == packages, f"exported {sorted(found)}, not {sorted(packages)}"
        # the bug fixed: DIS and NPF exported by the names the file gives
        assert not {"dis-x", "npf-1"} & found
        assert f"{expected['ghb']}_conc" in "".join(ds.variables)

    # read every array back with NETCDF, then run and compare heads
    set_netcdf_file(namefile, f"NETCDF FILEIN {name}.export.nc")
    for ext, names in ARRAYS.items():
        for path in ws.glob(f"{name}*.{ext}"):
            read_netcdf_input(path, names)
    success, buff = flopy.run_model(test.targets["mf6"], ws / "mfsim.nam", model_ws=ws)
    assert success, "run with NetCDF input failed:\n" + "\n".join(buff[-20:])
    heads = [
        flopy.utils.HeadFile(w / f"{name}.hds").get_alldata()
        for w in (test.workspace, ws)
    ]
    assert np.allclose(*heads)


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
@pytest.mark.parametrize("naming", ["named", "unnamed"])
def test_mf6model(idx, name, function_tmpdir, targets, export, naming):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, export, naming),
        compare=None,
    )
    test.run()
