"""
NetCDF input tests for test_gwf_vsc03_sfr: the simulation re-run with its package
arrays read from NetCDF, as exported by MODFLOW 6 and as written by flopy4.
"""

from pathlib import Path

import numpy as np
import pytest
from compare import Comparison
from framework import TestFramework
from netcdf_output import add_netcdf_output, check_model_netcdf_output
from test_gwf_vsc03_sfr import build_models, cases, check_output

nc = pytest.importorskip("netCDF4")


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare",
    [
        Comparison.MF6_EXPORT_STRUCTURED,
        Comparison.MF6_EXPORT_LAYERED,
        Comparison.FP4_STRUCTURED,
        Comparison.FP4_LAYERED,
    ],
)
def test_netcdf_input(idx, name, compare, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
    )
    test.run()


def build_output_models(idx, test, export):
    """The base simulation, writing its NetCDF output file."""
    sims = build_models(idx, test)
    model = sims[0].gwf[0]
    add_netcdf_output(model, export)
    return sims


def check_output_file(idx, test, export):
    """The base checks, and the NetCDF output file."""
    check_output(idx, test)
    sims = test.sims
    model = sims[0].gwf[0]
    check_model_netcdf_output(model)
    # z cell-center elevations, from a spatially varying top
    top = np.asarray(model.dis.top.array)
    botm = np.asarray(model.dis.botm.array)
    above = np.concatenate([top[np.newaxis], botm[:-1]])
    expected = (above + botm) / 2.0
    with nc.Dataset(Path(model.model_ws) / f"{model.name}.nc") as ds:
        if export == "structured":
            assert np.allclose(ds.variables["z"][:], expected)
        else:
            for k in range(botm.shape[0]):
                assert np.allclose(ds.variables[f"z_l{k + 1}"][:], expected[k].ravel())


@pytest.mark.netcdf
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("export", ["ugrid", "structured"])
def test_netcdf_output(idx, name, export, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_output_models(idx, t, export),
        check=lambda t: check_output_file(idx, t, export),
        compare=None,
    )
    test.run()
