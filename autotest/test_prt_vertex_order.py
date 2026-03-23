"""
Tests vertex ordering in quad-refined DISV grids.

This tests a fix to the logic that identifies the four corner points
of a rectangular cell that has one or more quad-refined neighbors, and
therefore one or more mid-face vertices. It specifically addresses the
case where a mid-face vertex is listed first in a cell2d block entry
for a cell.

Two cases are compared: vertex ordering in which a corner point is
always first, and a shifted ordering in which a mid-face point is first
for some cells. They should give identical particle tracks.
"""

from collections.abc import Iterable
from itertools import repeat
from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.binaryfile import CellBudgetFile, HeadFile
from flopy.utils.gridgen import Gridgen
from flopy.utils.gridintersect import GridIntersect
from framework import TestFramework, write_input
from prt_test_utils import check_budget_data, check_track_data, get_model_name
from shapely.geometry import LineString, MultiPoint

simname = "prtvord"
cases = [simname]

# Model parameters
nlay = 1
nrow = 7
ncol = 7
Lx = 3500.0
Ly = 3500.0
delr = Lx / ncol
delc = Ly / nrow
top = 400.0
botm = [0.0]
nper = 1
perlen = 1000.0
nstp = 1
tsmult = 1.0
porosity = 0.1
kh = 200.0
kv = 120.0
icelltype = [1]
rch = 0.005
rch_iface = 6
rch_iflowface = -1
riv_h = 320.0
riv_z = 318.0
riv_c = 1.0e5
riv_iface = 6
riv_iflowface = -1
wel_coords = [(1718.45, 1781.25)]
wel_q = [-150000.0]


def get_gridprops(test, shift=False):
    workspace = test.workspace
    targets = test.targets

    # Create base grid
    ms = flopy.modflow.Modflow()
    dis = flopy.modflow.ModflowDis(
        ms,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Create Gridgen workspace
    gridgen_ws = workspace / "gridgen"
    gridgen_ws.mkdir(parents=True, exist_ok=True)

    # Create Gridgen object
    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=targets["gridgen"],
    )

    # Add polygon for each refinement level
    outer_polygon = [[(500, 500), (500, 3000), (3000, 3000), (3000, 500), (500, 500)]]
    g.add_refinement_features([outer_polygon], "polygon", 1, range(nlay))

    middle_polygon = [
        [(1000, 1000), (1000, 2500), (2500, 2500), (2500, 1000), (1000, 1000)]
    ]
    g.add_refinement_features([middle_polygon], "polygon", 2, range(nlay))

    inner_polygon = [
        [(1500, 1500), (1500, 2000), (2000, 2000), (2000, 1500), (1500, 1500)]
    ]
    g.add_refinement_features([inner_polygon], "polygon", 3, range(nlay))

    # Build the grid
    g.build(verbose=False)
    disv_props = g.get_gridprops_disv()

    # Optionally shift vertex ordering
    if shift:
        cell2d = disv_props["cell2d"]
        for cell in cell2d:
            nv = cell[3] - 1  # number of vertices
            vlist = cell[4 : 4 + nv]
            # Rotate vertex list: move last vertex to front
            vlist.insert(0, vlist[-1])
            cell[4 : 4 + nv + 1] = vlist

    return disv_props


def build_gwf_sim(idx, test, mf6, shift=False):
    name = cases[idx]
    gwf_ws = test.workspace / "gwf"
    gwf_ws.mkdir(parents=True, exist_ok=True)
    gwf_name = get_model_name(name, "gwf")

    # Get grid properties
    disv_props = get_gridprops(test, shift=shift)

    # Create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name,
        exe_name=mf6,
        version="mf6",
        sim_ws=gwf_ws,
    )

    # Create TDIS
    tdis_rc = [(perlen, nstp, tsmult)]
    flopy.mf6.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # Create GWF model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        model_nam_file=f"{gwf_name}.nam",
    )
    gwf.name_file.save_flows = True

    # Create DISV
    flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units="FEET",
        **disv_props,
    )

    # GridIntersect for boundary conditions
    ix = GridIntersect(gwf.modelgrid, rtree=True)

    # Initial conditions
    flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=riv_h)

    # Node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=[("xt3d")],
        icelltype=icelltype,
        k=kh,
        k33=kv,
        save_saturation=True,
        save_specific_discharge=True,
    )

    # Recharge
    flopy.mf6.ModflowGwfrcha(
        gwf,
        recharge=rch,
        auxiliary=["iface", "iflowface"],
        aux=[rch_iface, rch_iflowface],
    )

    # Well
    welcells = ix.intersects(MultiPoint(wel_coords), dataframe=True)
    welcells = [row.cellid for row in welcells.itertuples()]
    welspd = [[(0, icpl), wel_q[i]] for i, icpl in enumerate(welcells)]
    flopy.mf6.ModflowGwfwel(gwf, print_input=True, stress_period_data=welspd)

    # River
    riverline = [(Lx - 1.0, Ly), (Lx - 1.0, 0.0)]
    rivcells = ix.intersects(LineString(riverline), dataframe=True)
    rivcells = [row.cellid for row in rivcells.itertuples()]
    rivspd = [
        [(0, icpl), riv_h, riv_c, riv_z, riv_iface, riv_iflowface] for icpl in rivcells
    ]
    flopy.mf6.ModflowGwfriv(
        gwf,
        stress_period_data=rivspd,
        auxiliary=[("iface", "iflowface")],
    )

    # Output control
    headfile = f"{gwf_name}.hds"
    budgetfile = f"{gwf_name}.cbb"
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="oc",
        budget_filerecord=budgetfile,
        head_filerecord=headfile,
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # IMS
    ims = flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1.0e-5,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=100,
        inner_dvclose=1.0e-6,
        rcloserecord=0.1,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=0.99,
    )
    sim.register_ims_package(ims, [gwf.name])

    return sim


def reverse_budgetfile(fpth, rev_fpth, tdis):
    f = CellBudgetFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)


def reverse_headfile(fpth, rev_fpth, tdis):
    f = HeadFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)


def build_prt_sim(idx, test, gwf_sim, mf6, shift=False):
    name = cases[idx]
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    prt_ws.mkdir(parents=True, exist_ok=True)
    prt_name = get_model_name(name, "prt")
    gwf_name = get_model_name(name, "gwf")

    # Get grid properties
    disv_props = get_gridprops(test, shift=shift)

    # Get GWF model for grid intersection
    gwf = gwf_sim.get_model(gwf_name)
    ix = GridIntersect(gwf.modelgrid, rtree=True)

    # Get well cell for particle release
    welcells = ix.intersects(MultiPoint(wel_coords), dataframe=True)
    welcells = [row.cellid for row in welcells.itertuples()]

    # Create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=prt_name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # Create TDIS
    tdis_rc = [(perlen, nstp, tsmult)]
    flopy.mf6.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # Create PRT model
    prt = flopy.mf6.ModflowPrt(
        sim,
        modelname=prt_name,
        model_nam_file=f"{prt_name}.nam",
    )

    # Create DISV
    flopy.mf6.ModflowGwfdisv(
        prt,
        length_units="FEET",
        **disv_props,
    )

    # MIP
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # PRP - release particles from well cell
    # Create particle data in local coordinates
    nodew = welcells[0]
    pcoord = np.array(
        [
            [0.000, 0.625, 0.500],
            [1.000, 0.375, 0.500],
            [0.375, 0.000, 0.500],
            [0.375, 1.000, 0.500],
        ]
    )
    partdata = flopy.modpath.ParticleData(
        [nodew for _ in range(pcoord.shape[0])],
        structured=False,
        localx=pcoord[:, 0],
        localy=pcoord[:, 1],
        localz=pcoord[:, 2],
        drape=0,
    )
    # Convert to PRP format (converts local to global coordinates)
    releasepts = list(partdata.to_prp(gwf.modelgrid))

    prpname = f"{prt_name}_1.prp"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=prpname,
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        exit_solve_tolerance=1e-5,
        extend_tracking=True,
    )

    # Output control
    budgetfile = f"{prt_name}.bud"
    trackfile = f"{prt_name}.trk"
    trackcsvfile = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budgetfile,
        track_filerecord=trackfile,
        trackcsv_filerecord=trackcsvfile,
        saverecord=[("BUDGET", "ALL")],
    )

    # FMI - use reversed budget and head files for backward tracking
    headfile_bkwd = f"{gwf_name}_bkwd.hds"
    budgetfile_bkwd = f"{gwf_name}_bkwd.cbb"
    pd = [
        ("GWFHEAD", Path(f"../{gwf_ws.name}/{headfile_bkwd}")),
        ("GWFBUDGET", Path(f"../{gwf_ws.name}/{budgetfile_bkwd}")),
    ]
    flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

    # EMS
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test, shift=False):
    gwf_sim = build_gwf_sim(idx, test, test.targets["mf6"], shift=shift)
    prt_sim = build_prt_sim(idx, test, gwf_sim, test.targets["mf6"], shift=shift)
    return gwf_sim, prt_sim


def check_output(idx, test, shift, snapshot):
    """Check PRT output and compare endpoints to snapshot."""
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    # Check output files exist
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()

    # Load MF6 pathlines
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # Check budget data
    check_budget_data(prt_ws / f"{name}_prt.lst", perlen, nper)

    # Check track data
    check_track_data(
        track_bin=prt_ws / prt_track_file,
        track_hdr=prt_ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=prt_ws / prt_track_csv_file,
    )

    # Extract endpoints and compare to snapshot
    mf6_eps = mf6_pls[mf6_pls.ireason == 3]
    assert snapshot == mf6_eps.round(2).to_records(index=False)


def plot_output(idx, test):
    """Plot PRT pathlines with head and velocity vectors."""
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # Load pathlines
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # Load head
    gwf_head_file = f"{gwf_name}.hds"
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()

    # Load specific discharge
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax.set_aspect("equal")

    # Plot MF6 pathlines
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    pmv.plot_grid(alpha=0.25)
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white", alpha=0.5)
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 PRT pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    plt.tight_layout()
    plt.savefig(prt_ws / f"test_{name}.png")
    plt.show()


@pytest.mark.snapshot
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("shift", [False, True])
def test_mf6model(idx, name, shift, function_tmpdir, targets, array_snapshot, plot):
    """Test PRT with normal and shifted vertex ordering."""

    def run(test):
        """Custom run sequence to handle file reversal for backward tracking."""
        sims = test.build(test)
        sims = sims if isinstance(sims, Iterable) else [sims]
        sims = [sim for sim in sims if sim]  # filter Nones
        test.sims = sims
        nsims = len(sims)
        test.buffs = list(repeat(None, nsims))

        write_input(*sims, overwrite=test.overwrite, verbose=test.verbose)

        gwf_sim, prt_sim = test.sims

        # Run GWF simulation
        success, buff = gwf_sim.run_simulation(silent=True, report=True)
        assert success, "GWF model run failed"

        # Reverse budget and head files for backward tracking
        gwf_ws = test.workspace / "gwf"
        gwf_name = get_model_name(name, "gwf")
        headfile = gwf_ws / f"{gwf_name}.hds"
        budgetfile = gwf_ws / f"{gwf_name}.cbb"
        headfile_bkwd = gwf_ws / f"{gwf_name}_bkwd.hds"
        budgetfile_bkwd = gwf_ws / f"{gwf_name}_bkwd.cbb"

        reverse_headfile(headfile, headfile_bkwd, gwf_sim.tdis)
        reverse_budgetfile(budgetfile, budgetfile_bkwd, gwf_sim.tdis)

        # Run PRT simulation
        success, buff = prt_sim.run_simulation(silent=False, report=True)
        assert success, "PRT model run failed"

        if test.check:
            if test.verbose:
                print("Checking results")
            test.check(test)

        if test.plot:
            if test.verbose:
                print("Plotting results")
            test.plot(test)

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t, shift=shift),
        check=lambda t: check_output(idx, t, shift=shift, snapshot=array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    run(test)
