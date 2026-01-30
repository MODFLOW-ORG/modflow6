"""
Minimal test for rotated DISV cells in PRT.

This test verifies that particle tracking through rotated DISV cells
correctly transforms coordinates. The bug being tested occurred when
composing/inverting transforms: a pure translation inversion was
incorrectly applying rotation to the stored origin.

Two cases are tested:

1. A simple (unrefined) rectangular DISV grid with rotated vertices.
2. A quad-refined DISV grid (via gridgen) with rotated vertices.

Both grids have their vertex coordinates rotated manually; the bug
causes coordinate corruption in rectilinear cells whose vertices are
rotated away from the axes. Rotation via the grid's `angrot` parameter
is not sufficient to trigger the bug, since that rotation is not used
for particle tracking, just metadata written to the grb file.

"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.binaryfile import HeadFile
from flopy.utils.gridgen import Gridgen
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtrotrect"
cases = [simname, f"{simname}q"]  # simple, quad-refined

# Use large coordinates to amplify the compose bug effect
# The bug causes larger errors when coordinates are far from origin
nlay = 1
nrow = 10
ncol = 10
Lx = 100000.0
Ly = 100000.0
delr = Lx / ncol
delc = Ly / nrow
top = 10.0
botm = [0.0]

# Non-zero rotation triggers the compose bug in quad cells
angle = 30.0  # degrees

# Expected grid bounds (after 30 degree rotation)
# Grid starts at origin (0,0) before rotation. After rotation:
# - Corner (0,0) -> (0,0)
# - Corner (Lx,0) -> (Lx*cos30, Lx*sin30) = (86602, 50000)
# - Corner (0,Ly) -> (-Ly*sin30, Ly*cos30) = (-50000, 86602)
# - Corner (Lx,Ly) -> (Lx*cos30-Ly*sin30, Lx*sin30+Ly*cos30) = (36602, 136602)
# Add small margin for numerical tolerance
margin = 1000.0
xmin = -Ly * np.sin(np.radians(angle)) - margin  # -50000 - margin
xmax = Lx * np.cos(np.radians(angle)) + margin  # 86602 + margin
ymin = -margin  # 0 - margin
ymax = (
    Lx * np.sin(np.radians(angle)) + Ly * np.cos(np.radians(angle)) + margin
)  # 136602 + margin


def rotate(x, y, angle_deg):
    angle_rad = np.radians(angle_deg)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    return x * cos_a - y * sin_a, x * sin_a + y * cos_a


def rotate_gridprops(gridprops):
    """Rotate all vertices and cell centers by the global angle."""
    vertices = gridprops["vertices"]
    rotated_vertices = []
    for v in vertices:
        iv, x, y = v[0], v[1], v[2]
        rx, ry = rotate(x, y, angle)
        rotated_vertices.append([iv, rx, ry])
    gridprops["vertices"] = rotated_vertices

    cell2d = gridprops["cell2d"]
    rotated_cell2d = []
    for cell in cell2d:
        icell = cell[0]
        xc, yc = cell[1], cell[2]
        rxc, ryc = rotate(xc, yc, angle)
        rotated_cell2d.append([icell, rxc, ryc] + cell[3:])
    gridprops["cell2d"] = rotated_cell2d

    return gridprops


def get_gridprops_quad(test):
    """Quad-refined DISV grid via gridgen, with rotated vertices."""
    workspace = test.workspace
    targets = test.targets

    ms = flopy.modflow.Modflow()
    flopy.modflow.ModflowDis(
        ms,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    gridgen_ws = workspace / "gridgen"
    gridgen_ws.mkdir(parents=True, exist_ok=True)

    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=targets["gridgen"],
    )

    polygon = [
        [(30000, 30000), (30000, 70000), (70000, 70000), (70000, 30000), (30000, 30000)]
    ]
    refinement_levels = 2
    g.add_refinement_features([polygon], "polygon", refinement_levels, range(nlay))
    g.build(verbose=False)

    return rotate_gridprops(g.get_gridprops_disv())


def get_gridprops_simple():
    """Simple (unrefined) rectangular DISV grid with rotated vertices."""
    gridprops = get_disv_kwargs(
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        tp=top,
        botm=botm,
    )
    return rotate_gridprops(gridprops)


def build_sim(idx, test):
    name = cases[idx]
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    ws = test.workspace

    gridprops = get_gridprops_simple() if idx == 0 else get_gridprops_quad(test)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=test.targets["mf6"], sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units="FEET",
        **gridprops,
    )

    flopy.mf6.ModflowGwfnpf(
        gwf,
        k=1.0,
        save_specific_discharge=True,
        save_saturation=True,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)

    # CHD: high head on one corner, low on opposite corner
    ncpl = gridprops["ncpl"]
    chd_data = [
        [(0, 0), top],
        [(0, ncpl - 1), 0.0],
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_data)

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    ims = flopy.mf6.ModflowIms(sim, print_option="SUMMARY")
    sim.register_solution_package(ims, [gwf.name])

    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    flopy.mf6.ModflowPrtdisv(
        prt,
        length_units="FEET",
        **gridprops,
    )

    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)

    # release particle from center of grid
    x_local = 50000.0
    y_local = 50000.0
    z_local = 5.0
    x_model, y_model = rotate(x_local, y_local, angle)
    icell = gwf.modelgrid.intersect(x_model, y_model)
    releasepts = [
        # particle index, (k, icell), x, y, z (model coordinates)
        (0, (0, icell), x_model, y_model, z_local),
    ]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )

    track_file = f"{prtname}.trk"
    track_csv = f"{prtname}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[track_file],
        trackcsv_filerecord=[track_csv],
    )

    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwfname,
        exgmnameb=prtname,
        filename=f"{gwfname}.gwfprt",
    )

    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test):
    return build_sim(idx, test)


def check_output(idx, test, snapshot):
    name = cases[idx]
    prtname = get_model_name(name, "prt")
    ws = test.workspace

    mf6_pls = pd.read_csv(ws / f"{prtname}.trk.csv")

    # check coordinates are within reasonable bounds
    x_ok = (mf6_pls["x"] >= xmin) & (mf6_pls["x"] <= xmax)
    y_ok = (mf6_pls["y"] >= ymin) & (mf6_pls["y"] <= ymax)
    z_ok = (mf6_pls["z"] >= -0.1) & (mf6_pls["z"] <= top + 0.1)
    bad_points = mf6_pls[~(x_ok & y_ok & z_ok)]
    assert len(bad_points) == 0, (
        f"Found {len(bad_points)} track points outside grid bounds. "
        f"This indicates the coordinate transform composition bug."
    )

    # check for the compose bug: with the bug, coordinates become erratic.
    # the particle should move generally in a consistent direction (toward low head).
    # check that x coordinates don't jump backwards by more than half a cell.
    x_vals = mf6_pls["x"].values
    for i in range(1, len(x_vals)):
        dx = x_vals[i] - x_vals[i - 1]
        assert dx >= -delr, (
            f"backward x jump at point {i}:"
            f"x[{i - 1}] = {x_vals[i - 1]:.2f}"
            f"x[{i}] = {x_vals[i]:.2f}"
            f"dx = {dx:.2f}"
        )

    # compare pathlines with snapshot
    actual_data = mf6_pls.drop("name", axis=1).round(3)
    actual_records = actual_data.to_records(index=False)
    assert snapshot == actual_records


def plot_output(idx, test):
    name = cases[idx]
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    ws = Path(test.workspace)
    sim = test.sims[0]
    gwf = sim.get_model(gwfname)
    mg = gwf.modelgrid

    # Load pathline results
    prt_track_csv = ws / f"{prtname}.trk.csv"
    mf6_pls = pd.read_csv(prt_track_csv).replace(r"^\s*$", np.nan, regex=True)

    # Load head and specific discharge
    gwf_head_file = ws / f"{gwfname}.hds"
    hds = HeadFile(gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # Set up plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax.set_aspect("equal")
    ax.set_title(f"Rotated DISV Grid - {name}")

    # Plot grid, head, and velocity vectors
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    pmv.plot_grid(alpha=0.5)
    pmv.plot_array(hds[0], alpha=0.2)
    pmv.plot_vector(qx, qy, normalize=True, color="gray", alpha=0.5)

    # Plot pathlines
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / max(len(mf6_plines), 1)),
            linewidth=2,
        )
        # Mark start and end points
        ax.plot(pl["x"].iloc[0], pl["y"].iloc[0], "go", markersize=8, label="Start")
        ax.plot(pl["x"].iloc[-1], pl["y"].iloc[-1], "ro", markersize=8, label="End")

    ax.legend()
    plt.savefig(ws / f"test_{name}.png")
    plt.show()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot, array_snapshot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
