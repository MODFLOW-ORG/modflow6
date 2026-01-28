"""
Minimal test for rotated quad-refined cells in PRT.

This test verifies that particle tracking through rotated quad-refined
cells correctly transforms coordinates. The bug being tested occurred
when composing/inverting transforms: a pure translation inversion was
incorrectly applying rotation to the stored origin.

The test uses gridgen to create a quad-refined DISV grid, then rotates
all vertices to create cells with non-zero sinrot. With cells that have
rotation AND quad subcells with non-zero origins, the bug causes
coordinate corruption.
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
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtrotrect"
cases = [simname]

# Grid parameters for base structured grid
# Use large coordinates to amplify the compose bug effect
# The bug causes larger errors when coordinates are far from origin
nlay = 1
nrow = 10
ncol = 10
Lx = 100000.0  # Large grid to trigger compose bug
Ly = 100000.0
delr = Lx / ncol
delc = Ly / nrow
top = 10.0
botm = [0.0]

# Rotation angle applied to vertices after gridgen creates the grid
# Non-zero rotation triggers the compose bug in quad cells
ROTATION_ANGLE = 30.0  # degrees

# Expected grid bounds (after 30 degree rotation)
# Grid starts at origin (0,0) before rotation. After rotation:
# - Corner (0,0) -> (0,0)
# - Corner (Lx,0) -> (Lx*cos30, Lx*sin30) = (86602, 50000)
# - Corner (0,Ly) -> (-Ly*sin30, Ly*cos30) = (-50000, 86602)
# - Corner (Lx,Ly) -> (Lx*cos30-Ly*sin30, Lx*sin30+Ly*cos30) = (36602, 136602)
# Add small margin for numerical tolerance
margin = 1000.0
XMIN = -Ly * np.sin(np.radians(ROTATION_ANGLE)) - margin  # -50000 - margin
XMAX = Lx * np.cos(np.radians(ROTATION_ANGLE)) + margin   # 86602 + margin
YMIN = -margin  # 0 - margin
YMAX = Lx * np.sin(np.radians(ROTATION_ANGLE)) + Ly * np.cos(np.radians(ROTATION_ANGLE)) + margin  # 136602 + margin


def rotate_point(x, y, angle_deg):
    """Rotate point around origin by angle in degrees."""
    angle_rad = np.radians(angle_deg)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    return x * cos_a - y * sin_a, x * sin_a + y * cos_a


def get_gridprops(test):
    """
    Create quad-refined grid using gridgen, then rotate vertices.

    Gridgen creates proper quad-refined cells with correct connectivity.
    Rotating the vertices afterward gives cells with non-zero sinrot,
    which is required to trigger the compose bug.
    """
    workspace = test.workspace
    targets = test.targets

    # Create base MODFLOW-2005 model for gridgen
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

    # Create gridgen workspace
    gridgen_ws = workspace / "gridgen"
    gridgen_ws.mkdir(parents=True, exist_ok=True)

    # Create Gridgen object
    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=targets["gridgen"],
    )

    # Add refinement polygon in center of grid to create quad-refined cells
    # Higher refinement levels create cells with 8 vertices (4 corners + 4 edge midpoints)
    # at the interface between different refinement zones
    polygon = [[(30000, 30000), (30000, 70000), (70000, 70000), (70000, 30000), (30000, 30000)]]
    refinement_levels = 2
    g.add_refinement_features(
        [polygon], "polygon", refinement_levels, range(nlay)
    )
    g.build(verbose=False)

    # Get gridprops from gridgen
    gridprops = g.get_gridprops_disv()

    # Now rotate all vertices by ROTATION_ANGLE degrees
    # This gives cells with non-zero sinrot, triggering the bug
    vertices = gridprops["vertices"]
    rotated_vertices = []
    for v in vertices:
        idx, x, y = v[0], v[1], v[2]
        rx, ry = rotate_point(x, y, ROTATION_ANGLE)
        rotated_vertices.append([idx, rx, ry])
    gridprops["vertices"] = rotated_vertices

    # Update cell centers to match rotated vertices
    cell2d = gridprops["cell2d"]
    rotated_cell2d = []
    for cell in cell2d:
        icell = cell[0]
        xc, yc = cell[1], cell[2]
        rxc, ryc = rotate_point(xc, yc, ROTATION_ANGLE)
        # Keep the rest of the cell data (nvert, vertex indices)
        rotated_cell2d.append([icell, rxc, ryc] + cell[3:])
    gridprops["cell2d"] = rotated_cell2d

    return gridprops


def build_sim(idx, test):
    """Build combined GWF-PRT simulation."""
    name = cases[idx]
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    ws = test.workspace

    # Get grid properties from gridgen (with rotated vertices)
    gridprops = get_gridprops(test)

    # Create simulation with both GWF and PRT models
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=test.targets["mf6"], sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])

    # --- GWF model ---
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
    # Find cells at corners after rotation
    ncpl = gridprops["ncpl"]
    chd_data = [
        [(0, 0), top],  # first cell, high head
        [(0, ncpl - 1), 0.0],  # last cell, low head
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

    # --- PRT model ---
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    flopy.mf6.ModflowPrtdisv(
        prt,
        length_units="FEET",
        **gridprops,
    )

    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)

    # Release particle from center of grid
    # Local coordinates (before rotation)
    x_local = 50000.0
    y_local = 50000.0
    z_local = 5.0

    # Transform to model coordinates (after rotation)
    x_model, y_model = rotate_point(x_local, y_local, ROTATION_ANGLE)

    # Find cell containing release point
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

    # Track file
    track_file = f"{prtname}.trk"
    track_csv = f"{prtname}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[track_file],
        trackcsv_filerecord=[track_csv],
    )

    # GWF-PRT exchange
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwfname,
        exgmnameb=prtname,
        filename=f"{gwfname}.gwfprt",
    )

    # Explicit model solution for PRT
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test):
    return build_sim(idx, test)


def check_output(idx, test):
    name = cases[idx]
    prtname = get_model_name(name, "prt")
    ws = test.workspace

    # Load track CSV
    track_csv = pd.read_csv(ws / f"{prtname}.trk.csv")

    print(f"Track data ({len(track_csv)} points):")
    print(track_csv[["x", "y", "z", "icell", "ireason"]].to_string())

    # Check all coordinates are within reasonable bounds
    x_ok = (track_csv["x"] >= XMIN) & (track_csv["x"] <= XMAX)
    y_ok = (track_csv["y"] >= YMIN) & (track_csv["y"] <= YMAX)
    z_ok = (track_csv["z"] >= -0.1) & (track_csv["z"] <= top + 0.1)

    bad_points = track_csv[~(x_ok & y_ok & z_ok)]
    if len(bad_points) > 0:
        print(f"\nFound {len(bad_points)} points outside grid bounds:")
        print(bad_points[["x", "y", "z"]].to_string())
        print(f"\nExpected bounds: x=[{XMIN:.0f}, {XMAX:.0f}], y=[{YMIN:.0f}, {YMAX:.0f}]")

    assert len(bad_points) == 0, (
        f"Found {len(bad_points)} track points outside grid bounds. "
        f"This indicates the coordinate transform composition bug."
    )

    # Check for the compose bug: with the bug, coordinates become erratic
    # The particle should move generally in a consistent direction (toward low head)
    # Check that x coordinates don't jump backwards by more than half a cell
    cell_size = delr  # ~10000 in this test
    x_vals = track_csv["x"].values
    for i in range(1, len(x_vals)):
        dx = x_vals[i] - x_vals[i - 1]
        # Large backward jumps indicate coordinate corruption
        if dx < -cell_size:
            print(f"\nDetected large backward x jump at point {i}:")
            print(f"  x[{i-1}] = {x_vals[i-1]:.2f}")
            print(f"  x[{i}] = {x_vals[i]:.2f}")
            print(f"  dx = {dx:.2f}")
            assert False, (
                f"Detected large backward x jump ({dx:.0f}) at point {i}. "
                f"This indicates the coordinate transform composition bug."
            )

    # Verify we have multiple track points (particle moved through cells)
    assert len(track_csv) >= 3, f"Expected at least 3 track points, got {len(track_csv)}"


def plot_output(idx, test):
    """Plot the grid, head distribution, and particle pathline."""
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
    ax.set_title(f"Rotated Quad-Refined Grid - {name}")

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
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
