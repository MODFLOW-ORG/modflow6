"""
Tests particle tracking in "dry" conditions.

The PRP package provides the `DRY` option
to specify how particles should behave in
dry conditions when the flow model enables
the Newton formulation.

This test case is adapted from the example
simulation provided by @javgs-bd in
https://github.com/MODFLOW-ORG/modflow6/issues/2014.
"""

import os
from warnings import warn

import flopy
import matplotlib.colors as clt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtdry"
cases = [
    # expect termination with status 7 immediately in 1st time step
    simname,
    # expect all particles to be draped to water table and tracked
    simname + "_drape",
    # the rest of the test cases activate newton and test the behavior
    # of the DRY option. with drop, expect all particles to be dropped
    # to the highest active cell below and then to be tracked as usual.
    simname + "_drop",
    # expect termination with status 7 immediately in 1st time step
    simname + "_stop",
    # expect particles to remain in release positions until the water
    # table rises to meet them
    simname + "_stay",
]

nper = 3
perlen = [1, 300, 1000]
nstp = [1, 30, 1]
tsmult = [1, 1.1, 1]

Lx = 100.0
Ly = 100.0

nlay = 3
ncol = 20
nrow = 20
arot = 0

delr = Lx / ncol
delc = Ly / nrow

strt = 1600

z_top = 1600
thickness = 10
xorigin, yorigin = (0, 0)

z_bot = np.zeros((nlay, nrow, ncol))
for i in range(nlay):
    z_bot[i, :, :] = z_top - (i + 1) * thickness

period_data = []
for i in range(nper):
    period_data.append((perlen[i], nstp[i], tsmult[i]))

user_time = 100.0

# particles are released in layer 0
offsets = [
    (-1, 0, 0),
    (-1, -1, 0),
    (-1, 1, 0),
    (-1, 0, -1),
    (-1, 0, 1),
]


def build_gwf_sim(name, gwf_ws, mf6, newton=False):
    sim = flopy.mf6.MFSimulation(
        sim_name=simname, exe_name=mf6, version="mf6", continue_=True, sim_ws=gwf_ws
    )

    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        model_nam_file=gwf_name + ".nam",
        newtonoptions="NEWTON" if newton else None,
        save_flows=True,
    )

    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="days", nper=nper, perioddata=period_data
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="METERS",
        xorigin=xorigin,
        yorigin=yorigin,
        angrot=arot,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=z_top,
        botm=z_bot,
        filename=gwf.name + ".dis",
        pname="dis",
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        csv_outer_output_filerecord="ims_outer.csv",
        csv_inner_output_filerecord="ims_inner.csv",
        outer_dvclose=1e-5,
        outer_maximum=100,
        under_relaxation="DBD",
        under_relaxation_gamma=0.01,
        under_relaxation_theta=0.7,
        under_relaxation_kappa=0.01,
        under_relaxation_momentum=0.0,
        inner_maximum=100,
        inner_dvclose=1e-6,
        rcloserecord=0.1,
        linear_acceleration="BICGSTAB",
        relaxation_factor=0.99,
        number_orthogonalizations=2,
        reordering_method="NONE",
        pname=gwf.name,
    )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf.name + ".cbb",
        head_filerecord=gwf.name + ".hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, strt=strt)

    bdry_name = "upstr.chd"

    filename = gwf.name + "." + bdry_name
    pname = bdry_name

    stress_period_data = {}
    lay = 1
    row, col = (0, 0)
    h_upstr = 1587

    for i in range(nper):
        stress_period_data[i] = [lay, row, col, h_upstr, bdry_name]

    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=stress_period_data,
        filename=filename,
        pname=pname,
        boundnames=True,
        save_flows=True,
    )

    bdry_name = "dwstr.chd"

    filename = gwf.name + "." + bdry_name
    pname = bdry_name

    stress_period_data = {}
    lay = 1
    row, col = (nrow - 1, ncol - 1)
    h_dwstr = 1582

    for i in range(nper):
        stress_period_data[i] = [lay, row, col, h_dwstr, bdry_name]

    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=stress_period_data,
        filename=filename,
        pname=pname,
        boundnames=True,
        save_flows=True,
    )

    bdry_name = "inj.wel"

    filename = gwf.name + "." + bdry_name
    pname = bdry_name

    stress_period_data = {}

    lay = 1
    row, col = (int(nrow / 4), int(ncol / 4))
    q = 100

    for i in range(1, nper):
        stress_period_data[i] = [lay, row, col, q, bdry_name]

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=stress_period_data,
        pname=pname,
        boundnames=True,
        save_flows=True,
    )

    icelltype = 1

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=icelltype,
        k=0.5,
        k33=0.1,
        save_specific_discharge=True,
        save_saturation=True,
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0001,
        sy=0.1,
        steady_state={0: True, 2: True},
        transient={1: True},
    )
    return sim


def build_prt_sim(name, gwf, prt_ws, mf6, drape=False, dry_tracking_method=False):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, version="mf6", continue_=True, sim_ws=prt_ws
    )

    gwf_ws = gwf.model_ws

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(
        sim,
        modelname=prt_name,
        model_nam_file=prt_name + ".nam",
    )

    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="days", nper=nper, perioddata=period_data
    )

    dis = flopy.mf6.ModflowPrtdis(
        prt,
        length_units="METERS",
        xorigin=xorigin,
        yorigin=yorigin,
        angrot=arot,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=z_top,
        botm=z_bot,
        filename=prt_name + ".dis",
        pname="dis",
    )

    porosity = 0.1

    flopy.mf6.ModflowPrtmip(prt, porosity=porosity, pname="mip")

    lay = 1
    row, col = (int(nrow / 4), int(ncol / 4))
    prp_cells = [(lay + k, row + i, col + j) for k, i, j in offsets]

    prp_coords = []
    prp_data = []

    for i, (lay, row, col) in enumerate(prp_cells):
        x, y, z = (
            gwf.modelgrid.xyzcellcenters[0][row, col],
            gwf.modelgrid.xyzcellcenters[1][row, col],
            gwf.modelgrid.xyzcellcenters[2][lay, row, col],
        )

        prp_lst = [i, (lay, row, col), x, y, z]

        prp_data.append(prp_lst)
        prp_coords.append((x, y, z))

    flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=len(prp_data),
        packagedata=prp_data,
        nreleasetimes=1,
        releasetimes=[(0.0,)],
        drape=drape,
        dry_tracking_method=dry_tracking_method,
        pname="prp",
        filename="tracking_1.prp",
        print_input=True,
    )

    budgetfile = prt.name + ".bud"
    trackfile = prt.name + ".trk"
    trackcsvfile = prt.name + ".csv"

    flopy.mf6.ModflowPrtoc(
        prt,
        budget_filerecord=[budgetfile],
        track_filerecord=[trackfile],
        trackcsv_filerecord=[trackcsvfile],
        saverecord=[("BUDGET", "ALL")],
        pname="oc",
        ntracktimes=1 if "stay" in name else 0,
        tracktimes=[(user_time,)] if "stay" in name else None,
    )

    rel_prt_folder = os.path.relpath(gwf_ws, start=prt_ws)

    packagedata = [
        ("GWFHEAD", f"{rel_prt_folder}/{gwf.name}.hds"),
        ("GWFBUDGET", f"{rel_prt_folder}/{gwf.name}.cbb"),
    ]

    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=packagedata,
    )

    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=prt.name + ".ems")

    return sim


def build_models(idx, test, newton, drape=False, dry_tracking_method=False):
    gwf_sim = build_gwf_sim(
        test.name, test.workspace / "gwf", test.targets["mf6"], newton=newton
    )
    prt_sim = build_prt_sim(
        test.name,
        gwf_sim.get_model(),
        test.workspace / "prt",
        test.targets["mf6"],
        drape=drape,
        dry_tracking_method=dry_tracking_method,
    )
    return gwf_sim, prt_sim


def check_output(idx, test, snapshot):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    hds_file = gwf.name + ".hds"
    trackcsv_file = prt_name + ".csv"
    trackcsv_path = prt_ws / trackcsv_file
    pls = pd.read_csv(trackcsv_path)
    strtpts = pls[pls.ireason == 0]

    # compare to expected results
    places = 1 if "drop" in name else 2
    actual = pls.drop(["name", "icell"], axis=1).round(places).reset_index(drop=True)
    nparts = len(offsets)  # number of particles

    if "drape" in name:
        assert len(actual[actual.ireason == 0]) == nparts  # release
    elif "drop" in name:
        # ignore particle 4, it terminates early when
        # mf6 is built with optimization=2 with ifort
        actual = actual.drop(actual[actual.irpt == 4].index)
        nparts -= 1
        assert len(actual[actual.ireason == 0]) == nparts  # release
    elif "stop" in name:
        assert len(actual[actual.ireason == 0]) == nparts  # release
    elif "stay" in name:
        assert len(actual[actual.ireason == 0]) == nparts  # release
        assert len(actual[actual.t == user_time]) == nparts  # user time
    else:
        # immediate termination, permanently unreleased
        assert len(actual) == nparts

    # in all cases, all particles should have a termination event
    assert len(actual[actual.ireason == 3]) == nparts

    # snapshot comparison
    assert snapshot == actual.to_records(index=False)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    hds_file = gwf.name + ".hds"
    trackcsv_file = prt_name + ".csv"
    trackcsv_path = prt_ws / trackcsv_file
    pls = pd.read_csv(trackcsv_path)
    strtpts = pls[pls.ireason == 0]

    def plot_pathlines_and_timeseries(
        ax, mg, ibd, pathlines, timeseries, plottitle, layer=1
    ):
        ax.set_aspect("equal")
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=layer)
        mm.plot_grid(color=(0.4, 0.4, 0.4, 0.5), lw=0.2)
        mm.plot_grid(color=(0.4, 0.4, 0.4, 0.5), lw=0.2)
        mm.plot_bc("WEL", plotAll=True)
        mm.plot_bc("CHD", plotAll=True)
        mm.plot_grid(lw=0.5)
        # mm.plot_array(gwf.output.head().get_data())
        v = mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")
        plt.scatter(pathlines["x"], pathlines["y"])
        mm.plot_pathline(pathlines, layer="all", colors=["blue"], lw=0.75)
        ax.set_title(plottitle, fontsize=12)
        ax.scatter(strtpts.x, strtpts.y)
        from shapely.geometry import Polygon

        cellids = [105, 85, 125, 104, 106]
        polys = [Polygon(gwf.modelgrid.get_cell_vertices(ic)) for ic in cellids]
        mm.plot_shapes(polys, alpha=0.2)
        plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    cmapbd = clt.ListedColormap(["b", "r"])
    chdcells = []
    for pckge in gwf.get_package("chd"):
        chd_nodes = pckge.stress_period_data.get_data(0).cellid
        for item in chd_nodes:
            chdcells.append(item)

    welcells = []
    pckg = gwf.get_package("wel")
    wel_nodes = pckg.stress_period_data.get_data(1).cellid
    for item in wel_nodes:
        welcells.append(item)

    # identify the boundary locations
    ibd = np.zeros((nlay, nrow, ncol), dtype=int)
    ilay, irow, icol = zip(*chdcells)
    ibd[ilay, irow, icol] = 1
    ilay, irow, icol = zip(*welcells)
    ibd[ilay, irow, icol] = 2
    ibd = np.ma.masked_equal(ibd, 0)

    plot_pathlines_and_timeseries(ax, gwf.modelgrid, ibd, pls, None, name)

    plot_3d = True
    if plot_3d:
        try:
            import pyvista as pv
            from flopy.export.vtk import Vtk
        except:
            warn("Couldn't make 3D plots, need pyvista/vtk")
            return

        vert_exag = 1
        vtk = Vtk(model=gwf, binary=False, vertical_exageration=vert_exag, smooth=False)
        vtk.add_model(gwf)
        vtk.add_pathline_points(pls.to_records(index=False))
        gwf_mesh, prt_mesh = vtk.to_pyvista()
        axes = pv.Axes(show_actor=False, actor_scale=2.0, line_width=5)
        p = pv.Plotter(window_size=[700, 700])
        p.enable_anti_aliasing()
        p.add_mesh(gwf_mesh, opacity=0.025, style="wireframe")
        p.add_mesh(
            prt_mesh,
            point_size=8,
            line_width=2.5,
            smooth_shading=True,
            color="blue",
        )
        p.show()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, array_snapshot, plot):
    dry_tracking_methods = ["drop", "stop", "stay"]
    if any(t in name for t in dry_tracking_methods):
        dry_tracking_method = name[-4:]
    else:
        dry_tracking_method = None
    newton = any(t in name for t in dry_tracking_methods)
    drape = "drape" in name
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t, newton, drape, dry_tracking_method),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
