"""
Test local x and y release point coordinates in the PRT PRP package.

With LOCAL_XY on a DIS grid, x and y are scaled to the unit interval
within the cell, 0 at the west and south faces and 1 at the east and
north faces. On a DISV grid the cell is divided into triangles joining
the cell's area centroid to each face. The integer part of x selects
the face from vertex i to vertex i + 1 and the fractional part is the
fraction along it. y is the fraction from the centroid (0) to the face (1).

With LOCAL_XY_OFFSET, x and y are offsets from the cell center, which
for DISV is the user-specified cell center.

Local coordinates are relative to the cell, so the grid's origin and
rotation don't affect where particles are released. Release positions
are reported in model coordinates.

Error cases: both options at once, and local coordinates out of range.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = [
    "prtlxy_dis",
    "prtlxy_disoff",
    "prtlxy_disv",
    "prtlxy_disvoff",
    "prtlxy_both",
    "prtlxy_disrng",
    "prtlxy_disvrng",
    "prtlxy_disvtri",
]
xfail = {
    "prtlxy_both": "may not both be specified",
    "prtlxy_disrng": "Local y coordinate must fall in the interval [0, 1]",
    "prtlxy_disvrng": "Local x coordinate must fall in the interval [1, 6)",
    "prtlxy_disvtri": "Local x coordinate must fall in the interval [1, 4)",
}

top, botm = 1.0, [0.0]
z = 0.5
georef = {"xorigin": 100.0, "yorigin": 200.0, "angrot": 30.0}

# DIS, non-uniform spacing. Rows run north to south.
delr = [1.0, 2.0, 4.0]
delc = [3.0, 5.0]
ymax = sum(delc)

# DISV, a 10x10 square cut into a pentagon and a triangle.
vertices = [
    (0, 0.0, 0.0),
    (1, 0.0, 10.0),
    (2, 6.0, 10.0),
    (3, 10.0, 4.0),
    (4, 10.0, 0.0),
    (5, 10.0, 10.0),
]
cells = [[0, 1, 2, 3, 4], [2, 5, 3]]  # clockwise


def cell_center(iverts):
    """Vertex average, deliberately not the area centroid."""
    return np.mean([vertices[i][1:] for i in iverts], axis=0)


def centroid(xy):
    x, y = xy[:, 0], xy[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    c = x * yn - xn * y
    a = c.sum() / 2
    return np.array([((x + xn) * c).sum(), ((y + yn) * c).sum()]) / (6 * a)


def fan(icell, x, y):
    xy = np.array([vertices[i][1:] for i in cells[icell]])
    i = int(x)
    t = x - i
    p = (1 - t) * xy[i - 1] + t * xy[i % len(xy)]
    c = centroid(xy)
    return c + y * (p - c)


def releases(name):
    """Release points (cellid, local x, local y) and expected model x, y."""
    if name in ["prtlxy_dis", "prtlxy_both"]:
        # cell (row 1, col 2) spans x [3, 7], y [0, 5]
        return [
            ((0, 1, 2), 0.25, 0.6, (4.0, 3.0)),
            ((0, 1, 2), 0.0, 0.0, (3.0, 0.0)),
            ((0, 1, 2), 1.0, 1.0, (7.0, 5.0)),
            ((0, 0, 0), 0.5, 0.5, (0.5, 6.5)),
        ]
    if name == "prtlxy_disoff":
        # cell (row 1, col 2) center (5, 2.5)
        return [
            ((0, 1, 2), -1.0, 2.0, (4.0, 4.5)),
            ((0, 0, 1), 0.0, 0.0, (2.0, 6.5)),
        ]
    if name == "prtlxy_disrng":
        return [((0, 1, 2), 0.5, 1.5, None)]
    if name == "prtlxy_disv":
        pts = [
            (0, 1.0, 0.0),  # centroid
            (0, 1.0, 1.0),  # vertex 1
            (0, 3.5, 1.0),  # middle of face 3
            (0, 5.5, 0.5),  # halfway to middle of the last face
            (0, 2.25, 0.8),
            (1, 1.0, 1.0),
            (1, 3.75, 0.3),
        ]
        return [((0, ic), x, y, tuple(fan(ic, x, y))) for ic, x, y in pts]
    if name == "prtlxy_disvoff":
        pts = [(0, 1.0, -2.0), (1, 0.0, 0.0)]
        return [
            ((0, ic), x, y, tuple(cell_center(cells[ic]) + (x, y))) for ic, x, y in pts
        ]
    if name == "prtlxy_disvrng":
        return [((0, 0), 6.0, 0.5, None)]
    if name == "prtlxy_disvtri":
        return [((0, 1), 4.0, 0.5, None)]


def build_models(test):
    name = test.name
    disv = "disv" in name
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=test.targets["mf6"], sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])

    def add_dis(model):
        if disv:
            cell2d = [
                [ic, *cell_center(iv), len(iv), *iv] for ic, iv in enumerate(cells)
            ]
            flopy.mf6.ModflowGwfdisv(
                model,
                nlay=1,
                ncpl=len(cells),
                nvert=len(vertices),
                top=top,
                botm=botm,
                vertices=vertices,
                cell2d=cell2d,
                **georef,
            )
        else:
            flopy.mf6.ModflowGwfdis(
                model,
                nlay=1,
                nrow=len(delc),
                ncol=len(delr),
                delr=delr,
                delc=delc,
                top=top,
                botm=botm,
                **georef,
            )

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    add_dis(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.bud", saverecord=[])
    ims = flopy.mf6.ModflowIms(sim, filename="gwf.ims")
    sim.register_solution_package(ims, [gwf.name])

    prt = flopy.mf6.ModflowPrt(sim, modelname="prt")
    add_dis(prt)
    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)
    options = {}
    if "off" in name:
        options["local_xy_offset"] = True
    else:
        options["local_xy"] = True
    if name == "prtlxy_both":
        options["local_xy_offset"] = True
    pts = releases(name)
    flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=len(pts),
        packagedata=[[i, *cid, x, y, z] for i, (cid, x, y, _) in enumerate(pts)],
        perioddata={0: [("FIRST",)]},
        **options,
    )
    flopy.mf6.ModflowPrtoc(prt, trackcsv_filerecord="prt.trk.csv")
    flopy.mf6.ModflowPrtfmi(prt)
    flopy.mf6.ModflowGwfprt(sim, exgmnamea=gwf.name, exgmnameb=prt.name)
    ems = flopy.mf6.ModflowEms(sim, filename="prt.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def check_output(test):
    name = test.name
    if name in xfail:
        assert any(xfail[name] in line for line in test.buffs[0])
        return
    pls = pd.read_csv(test.workspace / "prt.trk.csv")
    rel = pls[pls.ireason == 0].sort_values("irpt")
    expected = np.array([e for *_, e in releases(name)])
    assert np.allclose(rel[["x", "y"]].to_numpy(), expected)
    assert np.allclose(rel.z, z)


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build_models,
        check=check_output,
        targets=targets,
        xfail=name in xfail,
    )
    test.run()
