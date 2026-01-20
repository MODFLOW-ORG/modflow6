"""
Test of GWF DISU Package.  Use the flopy disu tool to create
a simple regular grid example, but using DISU instead of DIS.
The first case is just a simple test.  For the second case, set
one of the cells inactive and test to make sure connectivity
in binary grid file is correct.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework

cases = ["disu01a", "disu01b"]
grb_filename = "disu.grb"


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nlay = 3
    nrow = 3
    ncol = 3
    delr = 10.0 * np.ones(ncol)
    delc = 10.0 * np.ones(nrow)
    top = 0
    botm = [-10, -20, -30]
    disukwargs = get_disu_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
    if idx == 1:
        # for the second test, set one cell to idomain = 0
        idomain = np.ones((nlay, nrow * ncol), dtype=int)
        idomain[0, 1] = 0
        disukwargs["idomain"] = idomain

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    tdis = flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ims = flopy.mf6.ModflowIms(sim, print_option="SUMMARY")
    disu = flopy.mf6.ModflowGwfdisu(gwf, grb_filerecord=grb_filename, **disukwargs)
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf)
    spd = {0: [[(0,), 1.0], [(nrow * ncol - 1), 0.0]]}
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    return sim, None


def check_output(idx, test):
    fname = Path(test.workspace) / grb_filename
    grbobj = flopy.mf6.utils.MfGrdFile(fname)
    nodes = grbobj._datadict["NODES"]
    ia = grbobj._datadict["IA"]
    ja = grbobj._datadict["JA"]
    idomain = grbobj._datadict["IDOMAIN"]

    if idx == 0:
        assert np.array_equal(idomain, np.array(27 * [1]), int)

    if idx == 1:
        assert np.array_equal(ia[0:4], np.array([1, 4, 4, 7]))
        assert np.array_equal(ja[:6], np.array([1, 4, 10, 3, 6, 12]))
        assert ia[-1] == 127
        assert ia.shape[0] == 28, "ia should have size of 28"
        assert ja.shape[0] == 126, "ja should have size of 126"
        assert np.array_equal(idomain, np.array([1, 0] + 25 * [1]), int)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()
