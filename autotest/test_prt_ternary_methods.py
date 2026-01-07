"""
Tests a PRT model on the vertex grid demonstrated
at the end of Flopy's triangular mesh example:

https://flopy.readthedocs.io/en/latest/Notebooks/dis_triangle_example.html

There are two scenarios, both of which release
particles from the right border of the grid. In
the 1st case flow is left to right, in the 2nd
flow is top right to bottom left.

Runtime is benchmarked with pytest-benchmark.
The EXIT_SOLVE_METHOD option is used to select
root- finding methods for runtime comparison.
"""

from pathlib import Path

import flopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import get_model_name
from test_prt_triangle import (
    botm,
    build_gwf_sim,
    get_tri,
    nlay,
    porosity,
    top,
)

simname = "prtter"
cases = [
    f"{simname}br",
    f"{simname}ch",
]
methods = [
    1,  # brent
    2,  # chandrupatla
]


def build_prt_sim(idx, name, gwf_ws, prt_ws, targets, exit_solve_tolerance=1e-5):
    prt_ws = Path(prt_ws)
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=targets["mf6"], sim_ws=prt_ws
    )
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    tri = get_tri(prt_ws / "grid", targets)
    cell2d = tri.get_cell2d()
    vertices = tri.get_vertices()
    xcyc = tri.get_xcyc()
    ncpl = tri.ncpl
    nvert = tri.nvert
    dis = flopy.mf6.ModflowGwfdisv(
        prt,
        nlay=nlay,
        ncpl=ncpl,
        nvert=nvert,
        top=top,
        botm=botm,
        vertices=vertices,
        cell2d=cell2d,
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)
    prpdata = [
        # particle index, (layer, cell index), x, y, z
        (0, (0, 88), 95, 92, 0.5),
        (1, (0, 86), 96, 86, 0.5),
    ]
    prp_track_file = f"{prtname}.prp.trk"
    prp_track_csv_file = f"{prtname}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        boundnames=True,
        stop_at_weak_sink=True,  # currently required for this problem
        dev_exit_solve_method=methods[idx],
        exit_solve_tolerance=exit_solve_tolerance,
        extend_tracking=True,
    )
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )
    gwf_budget_file = gwf_ws / f"{gwfname}.cbc"
    gwf_head_file = gwf_ws / f"{gwfname}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(idx, test, exit_solve_tolerance=1e-7):
    gwf_sim = build_gwf_sim(test.name, test.workspace, test.targets, ["left", "botm"])
    prt_sim = build_prt_sim(
        idx,
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets,
        exit_solve_tolerance,
    )
    return gwf_sim, prt_sim


def plot_output(idx, test):
    name = test.name
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)

    # get gwf output
    gwf = gwf_sim.get_model()
    head = gwf.output.head().get_data()
    bdobj = gwf.output.budget()
    spdis = bdobj.get_data(text="DATA-SPDIS")[0]
    qx, qy, _ = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get prt output
    prt_name = get_model_name(name, "prt")
    prt_track_csv_file = f"{prt_name}.prp.trk.csv"
    pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)
    endpts = pls.sort_values("t").groupby(["imdl", "iprp", "irpt", "trelease"]).tail(1)

    # plot in 2d with mpl
    fig = plt.figure(figsize=(16, 10))
    ax = plt.subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pmv.plot_grid(alpha=0.25)
    pmv.plot_ibound(alpha=0.5)
    headmesh = pmv.plot_array(head, alpha=0.25)
    cv = pmv.contour_array(head, levels=np.linspace(0, 1, 9), colors="black")
    plt.clabel(cv)
    plt.colorbar(headmesh, shrink=0.25, ax=ax, label="Head", location="right")
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
    if "wel" in name:
        handles.append(
            mpl.lines.Line2D(
                [0], [0], marker="o", linestyle="", label="Well", markerfacecolor="red"
            )
        )
    ax.legend(handles=handles, loc="lower right")
    pmv.plot_vector(qx, qy, normalize=True, alpha=0.25)
    if "wel" in name:
        pmv.plot_bc(ftype="WEL")
    mf6_plines = pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        title = "DISV voronoi grid particle tracks"
        if "welp" in name:
            title += ": pumping wells"
        elif "weli" in name:
            title += ": injection wells"
        pl.plot(
            title=title,
            kind="line",
            linestyle="--",
            marker="o",
            markersize=2,
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color="black",
        )
    xc, yc = (
        gwf.modelgrid.get_xcellcenters_for_layer(0),
        gwf.modelgrid.get_ycellcenters_for_layer(0),
    )
    for i in range(gwf.modelgrid.ncpl):
        x, y = xc[i], yc[i]
        if i == 1639:
            color = "green"
            ms = 10
        else:
            color = "grey"
            ms = 2
        ax.plot(x, y, "o", color=color, alpha=0.25, ms=ms)
        ax.annotate(str(i + 1), (x, y), color="grey", alpha=0.5)

    plt.show()
    plt.savefig(prt_ws / f"{name}.png")


def check_output(idx, test, snapshot):
    name = test.name
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)

    # get gwf output
    gwf = gwf_sim.get_model()
    head = gwf.output.head().get_data()
    bdobj = gwf.output.budget()
    spdis = bdobj.get_data(text="DATA-SPDIS")[0]
    qx, qy, _ = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get prt output
    prt_name = get_model_name(name, "prt")
    prt_track_csv_file = f"{prt_name}.prp.trk.csv"
    pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)
    endpts = pls.sort_values("t").groupby(["imdl", "iprp", "irpt", "trelease"]).tail(1)

    # check pathline shape and endpoints
    assert pls.shape == (57, 16)
    assert endpts.shape == (2, 16)
    assert set(endpts.icell) == {111, 112}

    # check pathlines against snapshot
    assert snapshot == pls.drop("name", axis=1).round(3).to_records(index=False)


@pytest.mark.snapshot
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, benchmark, array_snapshot, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    benchmark(test.run)
