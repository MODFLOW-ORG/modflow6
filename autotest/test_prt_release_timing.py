"""
Test cases exercising release timing option.
including RELEASETIMES block release options
as well as period-block release options.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

Particles are released from the top left cell.

Results are compared against a MODPATH 7 model.
"""

from pathlib import Path
from typing import Optional

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.plot.plotutil import to_mp7_pathlines
from flopy.utils import PathlineFile
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from prt_test_utils import (
    FlopyReadmeCase,
    all_equal,
    check_budget_data,
    check_track_data,
    get_model_name,
    get_partdata,
)

simname = "prtrelt"
cases = [
    # options block options
    f"{simname}sgl",  # RELEASETIMES block: 0.5
    f"{simname}dbl",  # RELEASETIMES block: 0.5 and 0.6
    f"{simname}open",  # RELEASETIMES block: 0.5 and 0.6, OPEN/CLOSE
    # period block options
    f"{simname}all",  # ALL
    f"{simname}frac",  # ALL FRACTION 0.5, expect removal warning
    f"{simname}frst",  # FIRST
    f"{simname}stps",  # STEPS 1
    f"{simname}freq",  # FREQUENCY 1 and RELEASE_TIME_FREQUENCY 0.2
    # both releasetimes block and period block options
    f"{simname}both",  # RELEASETIMES block: 0.1; also FIRST
    f"{simname}dupe",  # RELEASETIMES block: 0.0: also FIRST, expect consolidation
    # test an absurdly high RELEASE_TIME_TOLERANCE
    f"{simname}tol",
    # test fill-forward with empty period block
    f"{simname}fill",  # FIRST in period 0, fill-forward to period 1, empty period 2
    # test different period block configs per period
    f"{simname}multi",  # FIRST in period 0, ALL in period 1, FIRST in period 2
]


def get_perioddata(name, periods=1) -> Optional[dict]:
    if "sgl" in name or "dbl" in name or "open" in name or "tol" in name:
        return None

    if "fill" in name:
        return {
            0: [("FIRST",)],  # Period 0: release on first time step
            # Period 1: omitted to test fill-forward
            2: [],  # Period 2: empty block to stop fill-forward
        }

    if "multi" in name:
        return {
            0: [("FIRST",)],  # Period 0: release on first time step only
            1: [("ALL",)],  # Period 1: release on all time steps
            2: [("FIRST",)],  # Period 2: release on first time step only
        }

    opt = []
    if "frst" in name or "both" in name or "dupe" in name:
        opt.append(("FIRST",))
    elif "all" in name:
        opt.append(("ALL",))
    elif "frac" in name:
        opt.append(("ALL",))
        opt.append(("FRACTION", 0.5))
    elif "stps" in name:
        opt.append(("STEPS", 1))
    elif "freq" in name:
        opt.append(("FREQUENCY", 1))
    else:
        opt.append(None)

    if opt[0] is None:
        raise ValueError(f"Invalid period option: {name}")

    return {i: opt for i in range(periods)}


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    # Use 3 periods for fill-forward and multi-period tests, 1 period for others
    nper = 3 if ("fill" in name or "multi" in name) else FlopyReadmeCase.nper
    if "multi" in name:
        # For multi test: period 1 has 5 steps so ALL is different from FIRST
        perioddata = [
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            ),  # Period 0: 1 time step
            (
                FlopyReadmeCase.perlen,
                5,
                FlopyReadmeCase.tsmult,
            ),  # Period 1: 5 time steps
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            ),  # Period 2: 1 time step
        ]
    elif "fill" in name:
        perioddata = [
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            )
        ] * nper
    else:
        perioddata = [
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            )
        ]
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=perioddata,
    )

    # create prt model
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

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

    # convert mp7 particledata to prt release points
    partdata = get_partdata(prt.modelgrid, FlopyReadmeCase.releasepts_mp7)
    releasepts = list(partdata.to_prp(prt.modelgrid))

    # check release points match expectation
    assert np.allclose(FlopyReadmeCase.releasepts_prt, releasepts)

    # create prp package
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    pdat = get_perioddata(prt_name)
    releasetimes = (
        [(0.0,)]
        if ("sgl" in name or "dupe" in name)
        else (
            [(0.1,)]
            if "both" in name
            else (
                [(0.0,), (0.1,)]
                if ("dbl" in name or "open" in name or "tol" in name)
                else None
            )
        )
    )
    releasetimes_path = prt_ws / "releasetimes.txt"
    if "open" in name:
        with open(releasetimes_path, "w") as f:
            for t in releasetimes:
                f.write(str(t[0]) + "\n")
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata=pdat,
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        nreleasetimes=(
            1
            if "sgl" in prt_name
            else 2
            if ("dbl" in name or "open" in name or "tol" in name)
            else None
        ),
        releasetimes=f"open/close {releasetimes_path.name}"
        if "open" in name
        else releasetimes
        if (
            "sgl" in name
            or "dbl" in name
            or "both" in name
            or "dupe" in name
            or "tol" in name
        )
        else None,
        release_time_frequency=0.2 if "freq" in name else None,
        print_input=True,
        extend_tracking=True,
        release_time_tolerance=0.2 if "tol" in name else None,
    )

    # create output control package
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
    gwf_name = get_model_name(name, "gwf")
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
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
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_mp7_sim(name, ws, mp7, gwf):
    partdata = get_partdata(gwf.modelgrid, FlopyReadmeCase.releasepts_mp7)
    mp7_name = get_model_name(name, "mp7")

    # Configure release times for multi-period tests to match MF6
    if "fill" in name:
        # Fill-forward test: release at t=0.0 and t=1.0
        # Format: [count, [list of times]]
        releasedata = [2, [0.0, 1.0]]
        pg = flopy.modpath.ParticleGroup(
            particlegroupname="G1",
            particledata=partdata,
            filename=f"{mp7_name}.sloc",
            releasedata=releasedata,
        )
    elif "multi" in name:
        # Multi-period test: release at 7 times
        # Format: [count, [list of times]]
        releasedata = [7, [0.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]]
        pg = flopy.modpath.ParticleGroup(
            particlegroupname="G1",
            particledata=partdata,
            filename=f"{mp7_name}.sloc",
            releasedata=releasedata,
        )
    else:
        # Default: release at start (t=0.0)
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


def build_models(test):
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )

    # For fill-forward and multi-period tests, update GWF simulation to use 3 periods
    if "fill" in test.name or "multi" in test.name:
        tdis = gwf_sim.get_package("tdis")
        tdis.nper = 3
        if "multi" in test.name:
            # For multi test: period 1 has 5 steps so ALL is different from FIRST
            tdis.perioddata = [
                (
                    FlopyReadmeCase.perlen,
                    FlopyReadmeCase.nstp,
                    FlopyReadmeCase.tsmult,
                ),  # Period 0: 1 time step
                (
                    FlopyReadmeCase.perlen,
                    5,
                    FlopyReadmeCase.tsmult,
                ),  # Period 1: 5 time steps
                (
                    FlopyReadmeCase.perlen,
                    FlopyReadmeCase.nstp,
                    FlopyReadmeCase.tsmult,
                ),  # Period 2: 1 time step
            ]
        else:
            tdis.perioddata = [
                (
                    FlopyReadmeCase.perlen,
                    FlopyReadmeCase.nstp,
                    FlopyReadmeCase.tsmult,
                )
            ] * 3

    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets["mf6"],
    )
    mp7_sim = build_mp7_sim(
        test.name,
        test.workspace / "mp7",
        test.targets["mp7"],
        gwf_sim.get_model(),
    )
    return gwf_sim, prt_sim, mp7_sim


def check_output(test, snapshot):
    name = test.name
    ws = test.workspace
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
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load head, budget and intercell flows from gwf model
    head = gwf.output.head().get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, _ = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # apply reference time to mp7 results (mp7 reports relative times)
    mp7_pls["time"] = mp7_pls["time"]

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # check output files exist
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (prt_ws / prp_track_file).is_file()
    assert (prt_ws / prp_track_csv_file).is_file()
    assert (mp7_ws / mp7_pathline_file).is_file()

    # check list file for logged release configuration
    list_file = prt_ws / f"{prt_name}.lst"
    assert list_file.is_file()
    lines = open(list_file).readlines()
    lines = [l.strip() for l in lines]
    if "frac" in name:
        # FRACTION no longer supported
        return
    else:
        li = lines.index("PARTICLE RELEASE FOR PRP 1")
        assert "RELEASE SCHEDULE:" in lines[li + 1]

    # make sure pathline df has "name" (boundname) column and empty values
    assert "name" in mf6_pls
    assert (mf6_pls["name"] == "").all()

    # make sure all mf6 pathline data have correct model and PRP index (1)
    assert all_equal(mf6_pls["imdl"], 1)
    assert all_equal(mf6_pls["iprp"], 1)

    # check budget data were written to mf6 prt list file
    # Multi-period tests use 3 periods, others use 1
    nper = 3 if ("fill" in name or "multi" in name) else FlopyReadmeCase.nper
    check_budget_data(
        prt_ws / f"{name}_prt.lst",
        FlopyReadmeCase.perlen,
        nper,
    )

    # check mf6 prt particle track data were written to binary/CSV files
    # and that different formats are equal
    for track_csv in [
        prt_ws / prt_track_csv_file,
        prt_ws / prp_track_csv_file,
    ]:
        check_track_data(
            track_bin=prt_ws / prt_track_file,
            track_hdr=prt_ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
            track_csv=track_csv,
        )

    # compare pathlines with snapshot
    assert snapshot == mf6_pls.drop("name", axis=1).round(3).to_records(index=False)

    # check release timing for fill-forward case
    if "fill" in name:
        # Period 0: FIRST at time = 0.0
        # Period 1: fill-forward, time = 1.0
        # Period 2: empty block should stop releases
        release_times = sorted(mf6_pls["trelease"].unique())
        expected_release_times = [0.0, FlopyReadmeCase.perlen]  # [0.0, 1.0]
        assert len(release_times) == len(expected_release_times)
        assert np.allclose(release_times, expected_release_times)

    # check release timing for case with several period blocks
    if "multi" in name:
        # Period 0: FIRST at time = 0.0
        # Period 1: ALL with 5 time steps at times 1.0, 1.2, 1.4, 1.6, 1.8
        # Period 2: FIRST at time = 2.0
        release_times = sorted(mf6_pls["trelease"].unique())
        expected_release_times = [
            0.0,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]
        assert len(release_times) == len(expected_release_times)
        assert np.allclose(release_times, expected_release_times)

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

    # compare mf6 / mp7 pathline results, todo check mass
    if "dbl" in name or "open" in name or "both" in name:
        # 2 release times, expect double mp7's data size
        assert len(mf6_pls) == 2 * len(mp7_pls)
    elif "freq" in name:
        # the "freq" case uses both time step frequency
        # in the period block and RELEASE_TIME_FREQUENCY
        # in the options block, setting the former to 1
        # and the latter to 0.2, so we expect 5 times as
        # many particles (first time t=0 is deduplicated)
        assert len(mf6_pls) == 5 * len(mp7_pls)
    elif "fill" in name or "multi" in name:
        # MP7 has duplicate records for particles released
        # at t = 1.0. these particles stop without moving,
        # but MP7 gives them 3 records, 2 in kper = 1 and
        # another in kper = 2. PRT terminates them at the
        # end of kper = 1 as per our philosophy that they
        # have been tracked under the prior period's flow
        # system, the next shouldn't "get credit" for it.
        assert len(mf6_pls) == len(mp7_pls) - 9
    else:
        # the rest of the cases should match mp7 results,
        # with duplicate times debounced in "dupe"/"tol".
        assert np.allclose(mf6_pls, mp7_pls, atol=1e-3)


def plot_output(test):
    name = test.name
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load head, budget and intercell flows from gwf model
    head = gwf.output.head().get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, _ = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # apply reference time to mp7 results (mp7 reports relative times)
    mp7_pls["time"] = mp7_pls["time"]

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # get time range for consistent colormap across both plots
    min_time = min(mf6_pls["t"].min(), mp7_pls["time"].min())
    max_time = max(mf6_pls["t"].max(), mp7_pls["time"].max())

    # plot mf6 pathlines in map view, colorcoded by time
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(head[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        # Plot lines connecting points
        ax[0].plot(pl["x"], pl["y"], "k-", linewidth=0.5, alpha=0.3)
        # Plot points colorcoded by time
        sc = ax[0].scatter(
            pl["x"],
            pl["y"],
            c=pl["t"],
            cmap="viridis",
            vmin=min_time,
            vmax=max_time,
            s=20,
            edgecolors="black",
            linewidth=0.5,
        )
    ax[0].set_title("MF6 pathlines (colored by time)")

    # plot mp7 pathlines in map view, colorcoded by time
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(head[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pls.groupby(["particleid"])
    for ipl, (pid, pl) in enumerate(mp7_plines):
        # Plot lines connecting points
        ax[1].plot(pl["x"], pl["y"], "k-", linewidth=0.5, alpha=0.3)
        # Plot points colorcoded by time
        sc = ax[1].scatter(
            pl["x"],
            pl["y"],
            c=pl["time"],
            cmap="viridis",
            vmin=min_time,
            vmax=max_time,
            s=20,
            edgecolors="black",
            linewidth=0.5,
        )
    ax[1].set_title("MP7 pathlines (colored by time)")

    # add colorbar to show time scale
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=0.05, aspect=40)
    cbar.set_label("Time")

    # view/save plot
    plt.show()
    plt.savefig(prt_ws / f"{name}.png")


@requires_pkg("syrupy")
@pytest.mark.snapshot
@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets, array_snapshot, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t),
        check=lambda t: check_output(t, array_snapshot),
        plot=lambda t: plot_output(t) if plot else None,
        targets=targets,
        compare=None,
        # expect case using FRACTION to fail
        xfail=[False, True, False] if "frac" in name else False,
    )
    test.run()
