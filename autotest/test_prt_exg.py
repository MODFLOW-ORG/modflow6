"""
Test GWF and PRT models in the same simulation
with an exchange.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.
Particles are released from the top left cell.

Results are compared against a MODPATH 7 model.

This test includes four cases: one which gives
boundnames to particles, one which does not, a
third in which the flow model has a uniformly
active idomain array while the tracking model
does not, and a final case in which flow and
tracking model have different IDOMAIN arrays,
both non-uniform, where the active region is
the same size but consists of different cells.
Both latter cases should be caught as errors.
"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import (
    FlopyReadmeCase,
    check_budget_data,
    check_track_data,
)

simname = "prtexg01"
cases = [simname, f"{simname}bnms", f"{simname}idmu", f"{simname}idmn"]


def get_model_name(idx, mdl):
    return f"{cases[idx]}_{mdl}"


def build_mf6_sim(idx, test):
    # create simulation
    name = cases[idx]
    sim = FlopyReadmeCase.get_gwf_sim(name, test.workspace, test.targets["mf6"])
    gwf = sim.get_model()

    # create prt model
    prt_name = get_model_name(idx, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    # create prt discretization
    idomain = np.ones(
        (FlopyReadmeCase.nlay, FlopyReadmeCase.nrow, FlopyReadmeCase.ncol)
    )
    if "idm" in name:
        # add an inactive cell to
        # tracking model idomain
        idomain[-1, -1, -1] = 0
    if "idmn" in name:
        # add a (different) inactive
        # cell to flow model idomain
        gwf_idomain = idomain.copy()
        gwf_idomain[-1, -1, -1] = 1
        gwf_idomain[0, 0, 0] = 0
        gwf.dis.idomain = gwf_idomain
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
        idomain=idomain,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # create prp package
    rpts = (
        [r + [str(r[0] + 1)] for r in FlopyReadmeCase.releasepts_prt]
        if "bnms" in name
        else FlopyReadmeCase.releasepts_prt
    )
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(rpts),
        packagedata=rpts,
        perioddata={0: ["FIRST"]},
        boundnames="bnms" in name,
        extend_tracking=True,
    )

    # create output control package
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        # dev_dump_event_trace=True
    )

    # create a flow model interface. should be ignored.
    # Mike Fienen reported that an early version of PRT
    # crashed if FMI is provided when using an exchange.
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            # garbage paths
            ("GWFHEAD", "heads.hds"),
            ("GWFBUDGET", "budget.cbc"),
        ],
    )

    # create exchange
    gwf_name = get_model_name(idx, "gwf")
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name,
        exgmnameb=prt_name,
        filename=f"{gwf_name}.gwfprt",
    )

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_mp7_sim(idx, ws, mp7, gwf):
    partdata = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in FlopyReadmeCase.releasepts_mp7],
        localx=[p[1] for p in FlopyReadmeCase.releasepts_mp7],
        localy=[p[2] for p in FlopyReadmeCase.releasepts_mp7],
        localz=[p[3] for p in FlopyReadmeCase.releasepts_mp7],
        timeoffset=0,
        drape=0,
    )
    mp7_name = get_model_name(idx, "mp7")
    pg = flopy.modpath.ParticleGroup(
        particlegroupname="G1",
        particledata=partdata,
        filename=f"{mp7_name}.sloc",
    )
    mp = flopy.modpath.Modpath7(
        modelname=mp7_name,
        flowmodel=gwf,
        exe_name=mp7,
        model_ws=ws,
        headfilename=f"{gwf.name}.hds",
        budgetfilename=f"{gwf.name}.bud",
    )
    mpbas = flopy.modpath.Modpath7Bas(mp, porosity=FlopyReadmeCase.porosity)
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        particlegroups=[pg],
    )
    return mp


def build_models(idx, test):
    mf6sim = build_mf6_sim(idx, test)
    gwf_name = get_model_name(idx, "gwf")
    gwf = mf6sim.get_model(gwf_name)
    mp7sim = build_mp7_sim(idx, test.workspace / "mp7", test.targets["mp7"], gwf)
    return mf6sim, None if "idm" in test.name else mp7sim


def check_output(idx, test):
    from flopy.plot.plotutil import to_mp7_pathlines

    name = test.name
    gwf_ws = Path(test.workspace)
    mp7_ws = gwf_ws / "mp7"

    # model names
    gwf_name = get_model_name(idx, "gwf")
    prt_name = get_model_name(idx, "prt")
    mp7_name = get_model_name(idx, "mp7")

    # extract model objects
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)
    prt = sim.get_model(prt_name)

    # extract model grid
    mg = gwf.modelgrid

    if "idm" in name:
        return

    # check mf6 output files exist
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (gwf_ws / prt_track_file).is_file()
    assert (gwf_ws / prt_track_csv_file).is_file()

    # check mp7 output files exist
    mp7_pathline_file = f"{mp7_name}.mppth"
    assert (mp7_ws / mp7_pathline_file).is_file()

    # load head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(gwf_name)
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mf6 pathline results
    mf6_pls = pd.read_csv(gwf_ws / prt_track_csv_file).replace(
        r"^\s*$", np.nan, regex=True
    )

    # make sure pathline dataframe has "name" column
    assert "name" in mf6_pls

    # check boundname values
    if "bnms" in name:
        # boundnames should be release point numbers (so pandas parses them as ints)
        assert np.array_equal(mf6_pls["name"].to_numpy(), mf6_pls["irpt"].to_numpy())
    else:
        # no boundnames given so check for defaults
        assert pd.isna(mf6_pls["name"]).all()

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

    # convert mf6 pathlines to mp7 format
    mf6_pls = to_mp7_pathlines(mf6_pls)

    # sort both dataframes by particleid and time
    mf6_pls.sort_values(by=["particleid", "time"], inplace=True)
    mp7_pls.sort_values(by=["particleid", "time"], inplace=True)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pls["sequencenumber"]
    del mf6_pls["particleidloc"]
    del mf6_pls["xloc"]
    del mf6_pls["yloc"]
    del mf6_pls["zloc"]
    del mp7_pls["sequencenumber"]
    del mp7_pls["particleidloc"]
    del mp7_pls["xloc"]
    del mp7_pls["yloc"]
    del mp7_pls["zloc"]

    # compare mf6 / mp7 pathline data
    assert mf6_pls.shape == mp7_pls.shape
    assert np.allclose(mf6_pls, mp7_pls, atol=1e-3)


def plot_output(idx, test):
    name = test.name
    gwf_ws = Path(test.workspace)
    mp7_ws = gwf_ws / "mp7"

    # model names
    gwf_name = get_model_name(idx, "gwf")
    prt_name = get_model_name(idx, "prt")
    mp7_name = get_model_name(idx, "mp7")

    # extract model objects
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)

    # extract model grid
    mg = gwf.modelgrid

    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(gwf_name)
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mf6 pathline results
    mf6_pls = pd.read_csv(gwf_ws / prt_track_csv_file).replace(
        r"^\s*$", np.nan, regex=True
    )

    # setup plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
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
            ax=ax[0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # plot mp7 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pls.groupby(["particleid"])
    for ipl, (pid, pl) in enumerate(mp7_plines):
        pl.plot(
            title="MP7 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[1],
            legend=False,
            color=cm.plasma(ipl / len(mp7_plines)),
        )

    # view/save plot
    plt.show()
    plt.savefig(gwf_ws / f"{name}.png")


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
        xfail="idm" in name,
    )
    test.run()
