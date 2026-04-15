"""
Regression test for issue #2772.

Build a 3x3x3 DISU model, deactivate one cell via IDOMAIN, then add an
HFB barrier whose CELLIDs use user (full-grid) node numbers that exceed
the active-cell count. Both endpoints of the barrier are active and
adjacent. The HFB IO guide documents CELLID for DISU as a user node
number, so this should run cleanly.

Prior to the fix in `Disu.f90:get_nodenumber_idx1`, the range check used
`this%nodes` (reduced/active count) rather than `this%nodesuser`,
inconsistent with the equivalent checks in DISV/DISV2D. Reading the HFB
file then failed with:

    Node number (N) is less than 1 or greater than nodes (n_active).
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework

cases = ["disuhfb"]


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nlay, nrow, ncol = 3, 3, 3
    delr = 10.0 * np.ones(ncol)
    delc = 10.0 * np.ones(nrow)
    top = 0.0
    botm = [-10.0, -20.0, -30.0]
    disukwargs = get_disu_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
    # Deactivate user node 2 so n_active = 26 < n_user = 27
    idomain = np.ones((nlay, nrow * ncol), dtype=int)
    idomain[0, 1] = 0
    disukwargs["idomain"] = idomain

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY")
    flopy.mf6.ModflowGwfdisu(gwf, **disukwargs)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    # CHD on first and last user cells (both active)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0,), 1.0], [(nlay * nrow * ncol - 1), 0.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )

    # HFB pair (26, 27): user IDs that span the active/inactive threshold.
    # Both cells are active and laterally adjacent in the bottom layer.
    hfb_data = [[(25,), (26,), 1.0e-6]]
    flopy.mf6.ModflowGwfhfb(gwf, stress_period_data={0: hfb_data})

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwf = sim.get_model(name)
    head = gwf.output.head().get_data()
    assert np.all(np.isfinite(head)), "head should be finite everywhere active"


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
