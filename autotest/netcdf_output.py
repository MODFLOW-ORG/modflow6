"""
Checks of a model's NetCDF output file against its binary output file.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest


def check_netcdf_output(nc_path: Path, binary_path: Path, text: str = "head"):
    """
    Check a NetCDF output file is NetCDF-4, holds no package input, and
    holds the dependent variable (`text`: head, concentration or
    temperature) of the binary output file at each time that file saved.
    The variable is one array in a structured file and one per layer
    (`<text>_l<k>`) in a mesh file. Also check the file decodes with
    xarray (CF times) and, for a mesh file, with xugrid (UGRID topology).
    """
    import netCDF4

    hds = flopy.utils.HeadFile(binary_path, precision="double", text=text.upper())
    binary = hds.get_alldata()
    ntimes, nlay = binary.shape[:2]
    binary = binary.reshape(ntimes, nlay, -1)

    with netCDF4.Dataset(nc_path) as ds:
        ds.set_auto_mask(False)
        assert ds.data_model == "NETCDF4", f"{nc_path.name}: {ds.data_model}"
        inputs = [v for v in ds.variables if "modflow_input" in ds[v].ncattrs()]
        assert not inputs, f"{nc_path.name} holds package input: {inputs}"
        mesh = text not in ds.variables
        if mesh:
            values = np.stack([ds[f"{text}_l{k + 1}"][:] for k in range(nlay)], axis=1)
        else:
            values = ds[text][:].reshape(ds[text].shape[0], nlay, -1)
        # the NetCDF file holds every time step, the binary file those saved
        times = ds["time"][:]
    steps = [int(np.argmin(np.abs(times - t))) for t in hds.get_times()]
    assert np.allclose(times[steps], hds.get_times()), (
        f"{nc_path.name}: binary output times not in the NetCDF file"
    )
    values = values[steps]
    assert values.shape == binary.shape, (
        f"{nc_path.name}: {text} shape {values.shape}, binary {binary.shape}"
    )
    assert np.allclose(values, binary), f"{nc_path.name}: {text} differs from binary"

    xr = pytest.importorskip("xarray")
    with xr.open_dataset(nc_path) as xds:
        assert np.issubdtype(xds["time"].dtype, np.datetime64), xds["time"].dtype
        if mesh:
            xu = pytest.importorskip("xugrid")
            assert xu.UgridDataset(xds).ugrid.grid.n_face == values.shape[-1]


# model type: (dependent variable, its output control file record)
DEPENDENT = {
    "gwf6": ("head", "head_filerecord"),
    "gwt6": ("concentration", "concentration_filerecord"),
    "gwe6": ("temperature", "temperature_filerecord"),
}


def add_netcdf_output(model, export: str) -> str:
    """Write a model's dependent variable to `<model>.nc` as well, in the
    structured or (`ugrid`) layered mesh format. Returns the file name."""
    fname = f"{model.name}.nc"
    if export == "ugrid":
        model.name_file.nc_mesh2d_filerecord = fname
    else:
        model.name_file.nc_structured_filerecord = fname
    return fname


def check_model_netcdf_output(model):
    """
    Check the NetCDF output file `add_netcdf_output` had a model write:
    its global and layer attributes, and its dependent variable against
    the model's binary output file (see `check_netcdf_output`).
    """
    import netCDF4

    ws = Path(model.model_ws)
    text, record = DEPENDENT[model.model_type.lower()]
    binary = getattr(model.get_package("oc"), record).get_data()[0][0]
    nc_path = ws / f"{model.name}.nc"
    with netCDF4.Dataset(nc_path) as ds:
        mtype = model.model_type.lower().removesuffix("6")
        expected = f"{mtype}: {model.name.lower()}"
        assert ds.getncattr("modflow_model") == expected, ds.getncattr("modflow_model")
        layer = ds.variables["layer"]
        assert layer.getncattr("units") == "1"
        assert layer.getncattr("axis") == "Z"
        assert layer.getncattr("positive") == "down"
        assert layer.getncattr("long_name") == "model layer"
        dependent = ds.variables[text if text in ds.variables else f"{text}_l1"]
        assert "standard_name" not in dependent.ncattrs()
        # head is a length; other dependent variables have no units known
        # from the input
        assert ("units" in dependent.ncattrs()) == (text == "head"), dependent.name
    check_netcdf_output(nc_path, ws / binary, text)
