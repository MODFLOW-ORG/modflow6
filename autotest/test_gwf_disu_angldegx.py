"""
Test the ANGLDEGX checks of the DISU Package.  The models are built with
the flopy disu tool, which writes consistent ANGLDEGX values, and then the
values are corrupted in one of two ways so that the run must issue an
informative warning message:

  1. the value for the connection from cell m to cell n is set equal to
     the value for the connection from cell n to cell m (the two outward
     normals must differ by 180 degrees), and
  2. the value for the connection from cell n to cell m is rotated by
     180 degrees so that the face normal points away from cell m (only
     detectable when VERTICES and CELL2D are specified).

The values are only warnings because ANGLDEGX is not used unless XT3D,
horizontal anisotropy, or SAVE_SPECIFIC_DISCHARGE is active.
"""

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework

cases = ["disu-anglsym", "disu-angldir"]


def build_models(idx, test):
    name = cases[idx]
    nlay, nrow, ncol = 1, 3, 3
    delr = 10.0 * np.ones(ncol)
    delc = 10.0 * np.ones(nrow)
    top = 0.0
    botm = [-10.0]
    disukwargs = get_disu_kwargs(
        nlay, nrow, ncol, delr, delc, top, botm, return_vertices=True
    )
    ja = np.array(disukwargs["ja"])
    angldegx = np.array(disukwargs["angldegx"], dtype=float)
    iac = np.array(disukwargs["iac"])
    ia = np.zeros(iac.shape[0] + 1, dtype=int)
    ia[1:] = np.cumsum(iac)

    def conn_index(n, m):
        for ipos in range(ia[n] + 1, ia[n + 1]):
            if ja[ipos] == m:
                return ipos
        raise ValueError(f"cells {n} and {m} are not connected")

    # corrupt the connection between cell 1 (row 0, col 1) and cell 4
    # (row 1, col 1); cell 4 is in the negative y direction from cell 1, so
    # the correct values are 270 degrees from cell 1 to cell 4 and 90
    # degrees from cell 4 to cell 1
    n, m = 1, 4
    if idx == 0:
        # forward connection gets the same normal as the reverse connection
        angldegx[conn_index(n, m)] = angldegx[conn_index(m, n)]
    else:
        # forward and reverse normals are still opposite, but both point
        # away from the connected cell
        angldegx[conn_index(n, m)] = 90.0
        angldegx[conn_index(m, n)] = 270.0
    disukwargs["angldegx"] = angldegx

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    tdis = flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ims = flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", linear_acceleration="BICGSTAB"
    )
    disu = flopy.mf6.ModflowGwfdisu(gwf, **disukwargs)
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf, xt3doptions=True)
    spd = {0: [[(0,), 1.0], [(nrow * ncol - 1,), 0.0]]}
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    return sim, None


def check_output(idx, test):
    if idx == 0:
        tag = (
            "ANGLDEGX values for 1 horizontal connections in the DISU Package "
            "are inconsistent with the values for the reverse connections "
            "(for example, ANGLDEGX = 90.000 for the connection from cell 2 "
            "to cell 5 and ANGLDEGX = 90.000 for the connection from cell 5 "
            "to cell 2)"
        )
    else:
        tag = (
            "ANGLDEGX values for 2 horizontal connections in the DISU Package "
            "point away from the connected cell (for example, ANGLDEGX = "
            "90.000 for the connection from cell 2 to cell 5, but the "
            "direction from the center of cell 2 to the center of cell 5 is "
            "270.000 degrees)"
        )
    # the warning report wraps long messages; join the lines before searching
    with open(test.workspace / "mfsim.lst", "r") as f:
        text = " ".join(line.strip() for line in f.readlines())
    assert "WARNING REPORT" in text, "no warning report in mfsim.lst"
    assert tag in text, f"expected warning message not found: {tag}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
