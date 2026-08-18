"""
Reproduce https://github.com/MODFLOW-ORG/modflow6/issues/2825.

A PRT particle release point on the edge of a cell could be
erroneously rejected as being outside the cell.

- "near" case: a distance safely inside PRT's tolerance but above the
  floating point noise floor. The simulation should complete normally.
- "far" case: a distance far outside PRT's tolerance. The simulation
  should fail with a "not in cell" error.
"""

import math

import flopy
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtedge"
cases = [f"{simname}_near", f"{simname}_far"]

nlay, nrow, ncol = 1, 5, 5
delr = delc = 10.0
top = 10.0
botm = [0.0]
angle = 41.3  # degrees; arbitrary non-axis-aligned rotation

xorigin, yorigin = 512345.678901, 4567890.123456

icell = 12  # an interior cell of the 5x5 grid

DIST_NEAR = 1.0e-6
DIST_FAR = 1.0e-3


def rotate(x, y, deg):
    r = math.radians(deg)
    c, s = math.cos(r), math.sin(r)
    return x * c - y * s, x * s + y * c


def get_gridprops():
    gridprops = get_disv_kwargs(
        nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, tp=top, botm=botm
    )

    vxy = {}
    vertices = []
    for iv, x, y in gridprops["vertices"]:
        rx, ry = rotate(x, y, angle)
        rx, ry = rx + xorigin, ry + yorigin
        vertices.append([iv, rx, ry])
        vxy[iv] = (rx, ry)
    gridprops["vertices"] = vertices

    cell2d = []
    for cell in gridprops["cell2d"]:
        ic, xc, yc = cell[0], cell[1], cell[2]
        rxc, ryc = rotate(xc, yc, angle)
        cell2d.append([ic, rxc + xorigin, ryc + yorigin] + cell[3:])
    gridprops["cell2d"] = cell2d

    return gridprops, vxy


def release_point(gridprops, vxy, distance):
    """A point offset outward from the midpoint of one of icell's
    edges by the given perpendicular distance."""
    ivs = gridprops["cell2d"][icell][4:]
    xc, yc = gridprops["cell2d"][icell][1], gridprops["cell2d"][icell][2]
    xa, ya = vxy[ivs[0]]
    xb, yb = vxy[ivs[1]]
    xmid, ymid = (xa + xb) / 2.0, (ya + yb) / 2.0
    edx, edy = xb - xa, yb - ya
    elen = math.hypot(edx, edy)
    nx, ny = edy / elen, -edx / elen
    if nx * (xc - xmid) + ny * (yc - ymid) > 0:
        nx, ny = -nx, -ny  # point outward, away from the cell center
    return xmid + distance * nx, ymid + distance * ny


def build_models(idx, test):
    name = cases[idx]
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    ws = test.workspace

    gridprops, vxy = get_gridprops()
    distance = DIST_NEAR if "near" in name else DIST_FAR
    x, y = release_point(gridprops, vxy, distance)

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=test.targets["mf6"], sim_ws=ws)
    # write coordinates at full double precision so the deliberately
    # tiny release point offset survives the round trip through the
    # input files
    sim.simulation_data.float_precision = 17
    sim.simulation_data.float_characters = 24
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **gridprops)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0), top], [(0, gridprops["ncpl"] - 1), 0.0]],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{gwfname}.cbc")
    ims = flopy.mf6.ModflowIms(sim)
    sim.register_solution_package(ims, [gwf.name])

    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    flopy.mf6.ModflowPrtdisv(prt, **gridprops)
    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)

    releasepts = [(0, (0, icell), x, y, 5.0)]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[f"{prtname}.trk"],
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
    )
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwfname,
        exgmnameb=prtname,
        filename=f"{gwfname}.gwfprt",
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prtname}.ems")
    sim.register_solution_package(ems, [prt.name])

    return sim


def build(idx, test):
    return build_models(idx, test)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build(idx, t),
        targets=targets,
        compare=None,
        xfail="far" in name,
    )
    test.run()
