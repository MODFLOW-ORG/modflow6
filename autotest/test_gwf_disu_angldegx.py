"""
Test the ANGLDEGX checks of the DISU Package.  The models are built with
the flopy disu tool, which writes consistent ANGLDEGX values, and then the
values are corrupted so that exactly one of the checks is triggered:

  anglsym   the value for the reverse connection is rotated by 30 degrees,
            so the two values no longer differ by 180 degrees; vertices are
            omitted so that the vertex-based checks cannot fire
  angldir   both values are rotated by 180 degrees, so they still differ by
            180 degrees but both point away from the connected cell
  anglvrt   both values are rotated by 20 degrees, so they still differ by
            180 degrees and still point into the connected cell, but differ
            by more than 1 degree from the normal computed from the shared
            vertices
  angl45    both values are rotated by 60 degrees, which is within 90
            degrees of the direction between the cell centers but beyond
            the 45-degree warning threshold

Each of these runs to completion with warnings because the NPF Package does
not use ANGLDEGX.  The last case repeats anglvrt with the XT3D option, in
which case the NPF Package terminates with an error.
"""

import re
from pathlib import Path

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework

cases = ["anglsym", "angldir", "anglvrt", "angl45", "anglerr"]

nlay, nrow, ncol = 1, 3, 3
# cell 1 (row 0, col 1) and cell 4 (row 1, col 1) in zero-based numbering;
# cell 4 is in the negative y direction from cell 1, so the correct values
# are 270 degrees from cell 1 to cell 4 and 90 degrees from cell 4 to cell 1
n_cell, m_cell = 1, 4

# fragments of the four warnings and the error, with the cell pair
tag_sym = (
    "ANGLDEGX values for 1 cell faces associated with horizontal connections "
    "in the DISU Package are inconsistent with the values for the reverse "
    "connections (for example, ANGLDEGX = 270.000 for the connection from "
    "cell 2 to cell 5 and ANGLDEGX = 120.000 for the connection from cell 5 "
    "to cell 2)"
)
tag_dir = (
    "ANGLDEGX values for 2 cell faces associated with horizontal connections "
    "in the DISU Package point away from the connected cell (for example, "
    "ANGLDEGX = 90.000 for the connection from cell 2 to cell 5, but the "
    "direction from the center of cell 2 to the center of cell 5 is 270.000 "
    "degrees)"
)
# the vertex-normal warning quotes the value given for the example
# connection, which differs from case to case
vrt_example = {"angldir": "90.000", "anglvrt": "290.000", "angl45": "330.000"}
vrt_example["anglerr"] = vrt_example["anglvrt"]


def tag_vrt(name):
    return (
        "ANGLDEGX values for 1 cell faces associated with horizontal "
        "connections in the DISU Package differ by more than 1 degree from the "
        "outward normal of the face computed from the two VERTICES shared by "
        "the connected cells (for example, ANGLDEGX = "
        f"{vrt_example.get(name, '290.000')} "
        "for the connection from cell 2 to cell 5, but the normal computed "
        "from the vertices is 270.000 degrees)"
    )


tag_45 = (
    "ANGLDEGX values for 2 cell faces associated with horizontal connections "
    "in the DISU Package deviate by more than 45 degrees from the direction "
    "between the centers of the connected cells (for example, ANGLDEGX = "
    "330.000 for the connection from cell 2 to cell 5, but the direction "
    "between the cell centers is 270.000 degrees)"
)
tag_err = "ANGLDEGX values in the DISU Package are inconsistent for 1 cell faces"

# which warning each case must and must not produce
expected = {
    "anglsym": {"sym"},
    "angldir": {"dir", "vrt"},
    "anglvrt": {"vrt"},
    "angl45": {"45", "vrt"},
    "anglerr": {"vrt"},
}


def build_models(idx, test):
    name = cases[idx]
    delr = 10.0 * np.ones(ncol)
    delc = 10.0 * np.ones(nrow)
    disukwargs = get_disu_kwargs(
        nlay, nrow, ncol, delr, delc, 0.0, [-10.0], return_vertices=True
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

    inm = conn_index(n_cell, m_cell)
    imn = conn_index(m_cell, n_cell)
    assert angldegx[inm] == 270.0 and angldegx[imn] == 90.0
    if name == "anglsym":
        angldegx[imn] = 120.0
        # no vertices, so only the reciprocity check applies
        for key in ("nvert", "vertices", "cell2d"):
            disukwargs[key] = None
    elif name == "angldir":
        angldegx[inm] = 90.0
        angldegx[imn] = 270.0
    elif name in ("anglvrt", "anglerr"):
        angldegx[inm] = 290.0
        angldegx[imn] = 110.0
    elif name == "angl45":
        angldegx[inm] = 330.0
        angldegx[imn] = 150.0
    disukwargs["angldegx"] = angldegx

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY", linear_acceleration="BICGSTAB")
    flopy.mf6.ModflowGwfdisu(gwf, **disukwargs)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, xt3doptions=(name == "anglerr"))
    spd = {0: [[(0,), 1.0], [(nrow * ncol - 1,), 0.0]]}
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    return sim, None


def listing_text(path):
    # the reports wrap long messages; join the lines before searching
    with open(path, "r") as f:
        return " ".join(line.strip() for line in f.readlines())


def check_output(idx, test):
    name = cases[idx]
    ws = Path(test.workspace)
    text = listing_text(ws / "mfsim.lst")
    tags = {
        "sym": tag_sym,
        "dir": tag_dir,
        "vrt": tag_vrt(name),
        "45": tag_45,
    }
    for key, tag in tags.items():
        if key in expected[name]:
            assert tag in text, f"{name}: expected warning not found: {tag}"
        else:
            assert tag not in text, f"{name}: unexpected warning: {tag}"
    # the warnings are also written to the model listing file immediately
    ntag = len(
        re.findall(r"WARNING: ANGLDEGX values for", listing_text(ws / f"{name}.lst"))
    )
    assert ntag == len(expected[name]), f"{name}: {ntag} warnings in model listing"
    if name == "anglerr":
        assert "ERROR REPORT" in text
        assert tag_err in text, f"{name}: expected error not found"
    else:
        assert "ERROR REPORT" not in text


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=(name == "anglerr"),
    )
    test.run()
