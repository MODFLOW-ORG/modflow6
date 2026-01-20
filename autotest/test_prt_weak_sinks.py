"""
GWF and PRT models run in separate simulations
via flow model interface.

The grid is a 10x10 square with 2 layers, based
on the flow system provided in the FloPy readme.
There is a well near the middle of the grid in
the top layer, which pumps at a very low rate.
Two test cases are defined, one with particle
release package (PRP) option STOP_AT_WEAK_SINK
disabled and one with the option enabled.

Particles are released from the top left cell.
With the STOP_AT_WEAK_SINK option enabled, the
well is expected to capture one particle. With
STOP_AT_WEAK_SINK disabled, the well no longer
captures the particle.

Results are compared against a MODPATH 7 model,
using WeakSinkOption 1 (pass-through) when the
STOP_AT_WEAK_SINK option is disabled, and when
it is enabled using WeakSinkOption 2 (stop-at).
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
    get_ireason_code,
    get_model_name,
)

simname = "prtfmi04"
cases = [simname, f"{simname}saws"]


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # output files
    gwfname = f"{name}_gwf"
    prtname = f"{name}_prt"
    gwf_budget_file = gwf_ws / f"{gwfname}.bud"
    gwf_head_file = gwf_ws / f"{gwfname}.hds"
    prt_track_file = prt_ws / f"{prtname}.trk"
    prt_track_csv_file = prt_ws / f"{prtname}.trk.csv"

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    pd = (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult)
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=FlopyReadmeCase.nper,
        perioddata=[pd],
    )

    # create prt model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # create prp package
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(FlopyReadmeCase.releasepts_prt),
        packagedata=FlopyReadmeCase.releasepts_prt,
        perioddata={0: ["FIRST"]},
        stop_at_weak_sink="saws" in name,
        extend_tracking=True,
    )

    # create output control package
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_mp7_sim(name, ws, mp7, gwf):
    mp7_name = f"{name}_mp7"
    partdata = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in FlopyReadmeCase.releasepts_mp7],
        localx=[p[1] for p in FlopyReadmeCase.releasepts_mp7],
        localy=[p[2] for p in FlopyReadmeCase.releasepts_mp7],
        localz=[p[3] for p in FlopyReadmeCase.releasepts_mp7],
        timeoffset=0,
        drape=0,
    )
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
    mpbas = flopy.modpath.Modpath7Bas(
        mp,
        porosity=FlopyReadmeCase.porosity,
    )
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        particlegroups=[pg],
        weaksinkoption="stop_at" if "saws" in name else "pass_through",
    )

    return mp


def get_different_rows(source_df, new_df):
    """Returns just the rows from the new dataframe that differ
    from the source dataframe
    """
    merged_df = source_df.merge(new_df, indicator=True, how="outer")
    changed_rows_df = merged_df[merged_df["_merge"] == "right_only"]
    return changed_rows_df.drop("_merge", axis=1)


def build_models(idx, test):
    # build gwf model
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )
    # add wel package
    gwf = gwf_sim.get_model()
    wells = [
        # k, i, j, q
        (0, 4, 4, -0.1),
    ]
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        maxbound=len(wells),
        save_flows=True,
        stress_period_data={0: wells, 1: wells},
    )

    # build prt model
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )

    # build mp7 model
    mp7_sim = build_mp7_sim(test.name, test.workspace / "mp7", test.targets["mp7"], gwf)
    return gwf_sim, prt_sim, mp7_sim


def check_output(idx, test):
    from flopy.plot.plotutil import to_mp7_pathlines

    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # check output files exist
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (mp7_ws / mp7_pathline_file).is_file()

    # if STOP_AT_WEAK_SINK disabled,
    # check for an extra datum when particle exited weak sink
    wksk_irsn = get_ireason_code("WEAKSINK")
    assert len(mf6_pls[mf6_pls["ireason"] == wksk_irsn]) == (
        1 if "saws" not in name else 0
    )
    # then drop the row so comparison will succeed below
    mf6_pls.drop(mf6_pls[mf6_pls["ireason"] == wksk_irsn].index, inplace=True)

    # make sure all mf6 pathline data have correct model and PRP index (1)
    def all_equal(col, val):
        a = col.to_numpy()
        return a[0] == val and (a[0] == a).all()

    assert all_equal(mf6_pls["imdl"], 1)
    assert all_equal(mf6_pls["iprp"], 1)

    # check budget data were written to mf6 prt list file
    check_budget_data(
        prt_ws / f"{name}_prt.lst",
        FlopyReadmeCase.perlen,
        FlopyReadmeCase.nper,
    )

    # check mf6 prt particle track data were written to binary/CSV files
    check_track_data(
        track_bin=prt_ws / prt_track_file,
        track_hdr=prt_ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=prt_ws / prt_track_csv_file,
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

    # drop node number column because prt and mp7 disagree on a few
    del mf6_pls["node"]
    del mp7_pls["node"]

    # compare mf6 / mp7 pathline data
    assert mf6_pls.shape == mp7_pls.shape
    assert np.allclose(mf6_pls, mp7_pls, atol=1e-3)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    pmv.plot_bc("WEL")
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
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    pmv.plot_bc("WEL")
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
    )
    test.run()
