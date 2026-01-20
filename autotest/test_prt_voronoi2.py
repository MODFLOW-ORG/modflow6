"""
Tests a PRT model on the Voronoi grid demonstrated
in Flopy's Voronoi example:

https://flopy.readthedocs.io/en/latest/Notebooks/dis_voronoi_example.html

Particles are released from the center of the plume
(i.e. the constant concentration cell) used in the
transport model.
"""

from pathlib import Path

import flopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.discretization import VertexGrid
from framework import TestFramework
from prt_test_utils import get_model_name
from test_prt_voronoi1 import build_gwf_sim, get_grid

simname = "prtvor2"
cases = [simname]
xmin = 0.0
xmax = 2000.0
ymin = 0.0
ymax = 1000.0
top = 1.0
botm = [0.0]
angle_min = 30
area_max = 1000.0
delr = area_max**0.5
nlay = 1
ncol = xmax / delr
nrow = ymax / delr
nodes = ncol * nrow
porosity = 0.1


def build_gwt_sim(name, gwf_ws, gwt_ws, targets):
    ws = Path(gwt_ws)
    gwf_name = get_model_name(name, "gwf")
    gwt_name = get_model_name(name, "gwt")

    # create grid
    grid = get_grid(ws / "grid", targets)
    vgrid = VertexGrid(**grid.get_gridprops_vertexgrid(), nlay=1)
    ibd = np.zeros(vgrid.ncpl, dtype=int)
    # gi = GridIntersect(vgrid)

    # identify release cell
    # point = Point((500, 500))
    # release_cells = gi.intersect(point)["cellids"]
    release_cells = [1294]
    release_cells = np.array(list(release_cells))

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", perioddata=[[100 * 365.0, 100, 1.0]]
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwt_name, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="simple",
        linear_acceleration="bicgstab",
        outer_dvclose=1.0e-6,
        inner_dvclose=1.0e-6,
    )
    disv_gridprops = grid.get_disv_gridprops()
    nlay = 1
    top = 1.0
    botm = [0.0]
    disv = flopy.mf6.ModflowGwtdisv(
        gwt, nlay=nlay, **disv_gridprops, top=top, botm=botm
    )
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=0.2)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=5.0, ath1=0.5)
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    cnclist = [
        [(0, release_cells[0]), 1.0],
    ]
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt, maxbound=len(cnclist), stress_period_data=cnclist, pname="CNC-1"
    )
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    flopy.mf6.ModflowGwtfmi(
        gwt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, {"release": release_cells}


def build_prt_sim(name, gwf_ws, prt_ws, targets, cell_ids):
    prt_ws = Path(prt_ws)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    # create grid
    grid = get_grid(prt_ws / "grid", targets)
    gridprops = grid.get_gridprops_vertexgrid()
    vgrid = VertexGrid(**gridprops, nlay=1)
    ibd = np.zeros(vgrid.ncpl, dtype=int)

    # identify cells on left edge
    left_cells = cell_ids["left"]
    left_cells = np.array(list(left_cells))
    ibd[left_cells] = 1

    # identify cells on right edge
    right_cells = cell_ids["right"]
    right_cells = np.array(list(right_cells))
    ibd[right_cells] = 2

    # identify release cell
    release_cells = cell_ids["release"]
    release_cells = np.array(list(release_cells))

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=targets["mf6"], sim_ws=prt_ws
    )
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    disv = flopy.mf6.ModflowGwfdisv(
        prt, nlay=nlay, **grid.get_disv_gridprops(), top=top, botm=botm
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    sddata = flopy.modpath.CellDataType(columncelldivisions=1, rowcelldivisions=1)
    data = flopy.modpath.NodeParticleData(subdivisiondata=sddata, nodes=[release_cells])
    prpdata = list(data.to_prp(prt.modelgrid))
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        boundnames=True,
        stop_at_weak_sink=True,  # currently required for this problem
        exit_solve_tolerance=1e-10,
        extend_tracking=True,
    )
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    grb_file = gwf_ws / f"{gwf_name}.disv.grb"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFGRID", grb_file),
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(idx, test):
    gwf_sim, gwf_cell_ids = build_gwf_sim(test.name, test.workspace, test.targets)
    gwt_sim, gwt_cell_ids = build_gwt_sim(
        test.name, test.workspace, test.workspace / "gwt", test.targets
    )
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets,
        {**gwf_cell_ids, **gwt_cell_ids},
    )
    return gwf_sim, gwt_sim, prt_sim


def check_output(idx, test):
    name = test.name
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(name, "prt")
    gwfsim, gwtsim, prtsim = test.sims

    # get gwf output
    gwf = gwfsim.get_model()
    head = gwf.output.head().get_data()
    bdobj = gwf.output.budget()
    spdis = bdobj.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get gwt output
    gwt = gwtsim.get_model()
    conc = gwt.output.concentration().get_data()

    # get prt output
    prt_track_csv_file = f"{prt_name}.prp.trk.csv"
    pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # pathlines should be west -> east and no elevation change
    assert len(pls.irpt.unique()) == 3
    assert pls.x.max() > 1950
    assert np.allclose(pls[pls.irpt == 1].z, 0.166666)
    assert np.allclose(pls[pls.irpt == 2].z, 0.5)
    assert np.allclose(pls[pls.irpt == 3].z, 0.833333)


def plot_output(idx, test):
    name = test.name
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(name, "prt")
    gwfsim, gwtsim, prtsim = test.sims

    # get gwf output
    gwf = gwfsim.get_model()
    head = gwf.output.head().get_data()
    bdobj = gwf.output.budget()
    spdis = bdobj.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get gwt output
    gwt = gwtsim.get_model()
    conc = gwt.output.concentration().get_data()

    # get prt output
    prt_track_csv_file = f"{prt_name}.prp.trk.csv"
    pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # plot in 2d with mpl
    fig = plt.figure(figsize=(16, 10))
    ax = plt.subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pmv.plot_grid(alpha=0.25)
    pmv.plot_ibound(alpha=0.5)
    # headmesh = pmv.plot_array(head, alpha=0.25)
    # headctr = pmv.contour_array(head, levels=np.linspace(0, 1, 9), colors="black")
    # plt.clabel(headctr)
    # plt.colorbar(headmesh, shrink=0.25, ax=ax, label="Head", location="right")
    concmesh = pmv.plot_array(conc, cmap="jet")
    concctr = pmv.contour_array(conc, levels=(0.0001, 0.001, 0.01, 0.1), colors="y")
    plt.clabel(concctr)
    plt.colorbar(concmesh, shrink=0.25, ax=ax, label="Concentration", location="right")

    handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker=">",
            linestyle="",
            label="Specific discharge",
            color="grey",
            markerfacecolor="gray",
        ),
    ]
    ax.legend(handles=handles, loc="lower right")
    pmv.plot_vector(qx, qy, normalize=True, alpha=0.25)
    mf6_plines = pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        title = "DISV voronoi grid particle tracks"
        pl.plot(
            title=title,
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color="black",
        )
    plt.show()
    plt.savefig(prt_ws / f"{name}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
