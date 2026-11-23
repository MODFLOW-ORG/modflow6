"""
Package input re-written with flopy4.

For `Comparison.FP4_STRUCTURED` and `Comparison.FP4_LAYERED`: a copy of a
test simulation is loaded and written again by flopy4, with each model's
package arrays in a NetCDF input file (`rewrite_netcdf`), and run. The
re-written input is checked to read every array MODFLOW 6 can read from
NetCDF that way (`check_netcdf_input`). flopy4 is imported only when a
simulation is re-written.
"""

from pathlib import Path

import mf6_input


def rewrite_netcdf(sim_ws: Path, netcdf_format: str):
    """
    Re-write the simulation in `sim_ws` in place with flopy4, with each
    model's package arrays in a NetCDF file, `<model>.input.nc`.

    Parameters
    ----------
    sim_ws : path-like
        Simulation workspace, a copy of the test's.
    netcdf_format : str
        "structured" (DIS grids only) or "layered" (UGRID layered mesh).
    """
    from flopy4.mf6 import NetCDFFormat, Simulation
    from flopy4.mf6.netcdf import NetCDFModel
    from flopy4.mf6.write_context import WriteContext

    sim_ws = Path(sim_ws)
    fmt = NetCDFFormat(netcdf_format)
    sim = Simulation.load(sim_ws / "mfsim.nam")
    for model in sim.models:
        name = model.name
        if fmt == NetCDFFormat.STRUCTURED and type(model.dis).__name__ != "Dis":
            raise ValueError(f"Structured NetCDF input requires a DIS grid: {name}")
        # the name file refers to it relative to the workspace
        nc_name = f"{name}.input.nc"
        model.netcdf_input_file = Path(nc_name)
        NetCDFModel.from_model(model, netcdf_format=fmt).to_netcdf(sim_ws / nc_name)
    with WriteContext(use_netcdf=True):
        sim.write()


def rewrite_ascii(sim_ws: Path):
    """Re-write the simulation in `sim_ws` in place with flopy4, as text
    input."""
    from flopy4.mf6 import Simulation

    Simulation.load(Path(sim_ws) / "mfsim.nam").write()


def package_arrays(path: Path) -> list[tuple[str, bool]]:
    """(name, read from NetCDF) of each array given in a package file."""
    rows = [t for t in map(mf6_input.tokens, mf6_input.read_lines(path)) if t]
    arrays = []
    for i, row in enumerate(rows):
        upper = [t.upper() for t in row]
        if (
            len(upper) > 1
            and upper[-1] == "NETCDF"
            and upper[0] not in ("BEGIN", "END")
        ):
            arrays.append((row[0].lower(), True))
        elif (
            upper[1:] in ([], ["LAYERED"])
            and i + 1 < len(rows)
            and rows[i + 1][0].upper() in ("CONSTANT", "INTERNAL", "OPEN/CLOSE")
        ):
            arrays.append((row[0].lower(), False))
    return arrays


def check_netcdf_input(sim_ws: Path):
    """
    Check a simulation re-written by `rewrite_netcdf` reads every array
    MODFLOW 6 can read from NetCDF (per the DFNs) with NETCDF, from a
    variable in its model's NetCDF4 input file.
    """
    import netCDF4

    sim_ws = Path(sim_ws)
    models = mf6_input.models(sim_ws)
    assert models, f"No models in {sim_ws / 'mfsim.nam'}"
    for mtype, namefile, mname in models:
        nc = [
            t[2]
            for t in map(mf6_input.tokens, mf6_input.read_lines(sim_ws / namefile))
            if [w.upper() for w in t[:2]] == ["NETCDF", "FILEIN"]
        ]
        assert nc, f"No NETCDF FILEIN in {namefile}"
        with netCDF4.Dataset(sim_ws / nc[0]) as ds:
            assert ds.data_model == "NETCDF4", f"{nc[0]}: {ds.data_model}"
            assert "modflow_model" in ds.ncattrs(), f"{nc[0]}: no modflow_model"
            inputs = {
                ds[v].modflow_input.lower()
                for v in ds.variables
                if "modflow_input" in ds[v].ncattrs()
            }
        n_netcdf = 0
        for ptype, fname, pname in mf6_input.packages(sim_ws, namefile, mtype):
            lines = mf6_input.read_lines(sim_ws / fname)
            component = mf6_input.component(mtype, ptype, lines)
            if component is None:
                continue
            capable = mf6_input.netcdf_arrays(component)
            # variables are named by package type, or by package name for
            # types a model can have several of
            package = pname if mf6_input.multi_package(component) else ptype
            aux = set(mf6_input.aux_names(lines))
            for name, from_netcdf in package_arrays(sim_ws / fname):
                tag = "aux" if name in aux else name
                if from_netcdf:
                    n_netcdf += 1
                    key = f"{mname}/{package}/{tag}".lower()
                    assert key in inputs, f"{fname}: {name} has no variable {key}"
                else:
                    assert tag not in capable, f"{fname}: {name} not read from NetCDF"
        assert n_netcdf, f"No package array of {namefile} read from NetCDF"
