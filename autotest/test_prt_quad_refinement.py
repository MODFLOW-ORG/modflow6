"""
Tests quad-refined grids with 1 and 2 refinement levels.
The 2-level case has an abrupt transition between levels
with no smoothing.

The flow system is based on the FloPy README example.
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
from prt_test_utils import (
    FlopyReadmeCase,
    check_budget_data,
    check_track_data,
)

simname = "prtqref01"
cases = [simname]


def get_model_name(idx, mdl):
    return f"{cases[idx]}_{mdl}"


def get_gridprops(test, **kwargs):
    workspace = test.workspace
    targets = test.targets

    Lx = 1000.0
    Ly = 1000.0
    nlay = 1
    nrow = 10
    ncol = 10
    delr = Lx / ncol
    delc = Ly / nrow
    top = 10
    botm = [0]

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

    # create Gridgen workspace
    gridgen_ws = workspace / "gridgen"
    gridgen_ws.mkdir(parents=True, exist_ok=True)

    refinement_levels = kwargs.pop("refinement_levels")

    # create Gridgen object
    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=targets["gridgen"],
        **kwargs,
    )

    # add polygon for each refinement level
    polygon = [[(300, 300), (300, 700), (700, 700), (700, 300), (300, 300)]]
    g.add_refinement_features([polygon], "polygon", refinement_levels, range(nlay))
    g.build(verbose=False)
    return g.get_gridprops_disv()


def build_mf6_sim(idx, test, **kwargs):
    # model names
    gwf_name = get_model_name(idx, "gwf")
    prt_name = get_model_name(idx, "prt")

    tracking_method = kwargs.pop("tracking_method")

    # create refined grid
    gridprops = get_gridprops(test, **kwargs)

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, exe_name="mf6", version="mf6", sim_ws=test.workspace
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]], nper=1
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, model_nam_file=gwf_name)
    gwf.name_file.save_flows = True
    disv = flopy.mf6.ModflowGwfdisv(
        gwf,
        length_units="FEET",
        **gridprops,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=5)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_saturation=True,
        save_specific_discharge=True,
    )
    spd = {
        0: [[(0, 0), 1.0], [(0, gwf.modelgrid.intersect(950, 50)), 0.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=spd,
    )
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_budget_file,
        head_filerecord=gwf_head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    ims = flopy.mf6.ModflowIms(sim)
    sim.register_solution_package(ims, [gwf.name])

    # create prt model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdisv(prt, length_units="FEET", **gridprops)
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)
    rpts = [(50, 950), (45, 945), (55, 955)]
    rpts = [
        [i, 0, prt.modelgrid.intersect(x, y), x, y, 5.0]
        for i, (x, y) in enumerate(rpts)
    ]
    flopy.mf6.ModflowPrtprp(
        prt,
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(rpts),
        packagedata=rpts,
        perioddata={0: ["FIRST"]},
        dev_forceternary=tracking_method == "ternary",
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
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name,
        exgmnameb=prt_name,
        filename=f"{gwf_name}.gwfprt",
    )
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test, **kwargs):
    return build_mf6_sim(idx, test, **kwargs)


def check_output(idx, test, snapshot):
    name = test.name
    gwf_ws = Path(test.workspace)
    mp7_ws = gwf_ws / "mp7"
    gwf_name = get_model_name(idx, "gwf")
    prt_name = get_model_name(idx, "prt")
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)
    prt = sim.get_model(prt_name)
    mg = gwf.modelgrid
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(gwf_ws / prt_track_csv_file).replace(
        r"^\s*$", np.nan, regex=True
    )

    # extract head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(gwf_name)
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # check mf6 output files exist
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (gwf_ws / prt_track_file).is_file()
    assert (gwf_ws / prt_track_csv_file).is_file()

    # check budget data were written to mf6 prt list file
    check_budget_data(
        gwf_ws / f"{name}_prt.lst",
        FlopyReadmeCase.perlen,
        FlopyReadmeCase.nper,
    )

    # check mf6 prt particle track data were written to binary/CSV files
    check_track_data(
        track_bin=gwf_ws / prt_track_file,
        track_hdr=gwf_ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=gwf_ws / prt_track_csv_file,
    )

    # extract endpoints and compare to snapshot
    mf6_eps = mf6_pls[mf6_pls.ireason == 3]
    assert snapshot == mf6_eps.round(2).to_records(index=False)


def plot_output(idx, test):
    name = test.name
    gwf_ws = Path(test.workspace)
    gwf_name = get_model_name(idx, "gwf")
    prt_name = get_model_name(idx, "prt")
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)
    mg = gwf.modelgrid
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(gwf_ws / prt_track_csv_file).replace(
        r"^\s*$", np.nan, regex=True
    )

    # extract head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(gwf_name)
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # plot nodes
    xc, yc = (mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0))
    for i in range(mg.ncpl):
        x, y = xc[i], yc[i]
        ax.annotate(str(i + 1), (x, y), color="grey", alpha=0.5)

    # view/save plot
    plt.show()
    plt.savefig(gwf_ws / f"test_{name}.png")


@pytest.mark.snapshot
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize("levels", [1, 2])
@pytest.mark.parametrize("method", ["pollock", "ternary"])
def test_mf6model(
    idx, name, function_tmpdir, targets, levels, method, array_snapshot, plot
):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(
            idx,
            t,
            tracking_method=method,
            refinement_levels=levels,
            smoothing_level_vertical=levels,
            smoothing_level_horizontal=levels,
        ),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
