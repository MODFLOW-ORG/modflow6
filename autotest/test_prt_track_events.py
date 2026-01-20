"""
This test exercises TRACKEVENT options to check
that tracking event selection works as expected.

GWF and PRT models run in separate simulations.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme,
except for 2 inactive cells in the bottom left
and top right corners.

The flow system is similar to test_prt_fmi01.py.
Particles are split across two release packages,
and the grid has an inactive region this time,
to check cell numbers recorded in pathline data
are converted from reduced to user node numbers.
This is verified with FloPy by intersecting path
points with the grid then computing node numbers.

Particles are released from the top left cell.

Pathlines are compared with a MODPATH 7 model.
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
    get_model_name,
)

simname = "prtevnt"
cases = [
    f"{simname}all",
    f"{simname}rel",
    f"{simname}tsit",
    f"{simname}tstp",
    f"{simname}term",
    f"{simname}wksk",
    f"{simname}mult",
    f"{simname}trts",
    f"{simname}open",
]
releasepts_prt = {
    "a": [
        # index, k, i, j, x, y, z
        [i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5]
        for i in range(4)
    ],
    "b": [
        # index, k, i, j, x, y, z
        [i, 0, 0, 0, float(f"0.{i + 5}"), float(f"9.{i + 5}"), 0.5]
        for i in range(5)
    ],
}
releasepts_mp7 = {
    "a": [
        # node number, localx, localy, localz
        (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
        for i in range(4)
    ],
    "b": [
        # node number, localx, localy, localz
        (0, float(f"0.{i + 5}"), float(f"0.{i + 5}"), 0.5)
        for i in range(5)
    ],
}
tracktimes = list(np.linspace(0, 50, 1000))


# function to create idomain from grid dimensions
def create_idomain(nlay, nrow, ncol):
    idmn = np.ones((nlay, nrow, ncol), dtype=int)
    idmn[0, 0, 9] = 0
    idmn[0, 9, 0] = 0
    return idmn


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=FlopyReadmeCase.nper,
        perioddata=[
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            )
        ],
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
        idomain=create_idomain(
            FlopyReadmeCase.nlay, FlopyReadmeCase.nrow, FlopyReadmeCase.ncol
        ),
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # create a prp package for groups a and b
    prps = [
        flopy.mf6.ModflowPrtprp(
            prt,
            pname=f"prp_{grp}",
            filename=f"{prt_name}_{grp}.prp",
            nreleasepts=len(releasepts_prt[grp]),
            packagedata=releasepts_prt[grp],
            perioddata={0: ["FIRST"]},
            extend_tracking=True,
        )
        for grp in ["a", "b"]
    ]

    def get_oc() -> list[str]:
        prt_track_file = f"{prt_name}.trk"
        prt_track_csv_file = f"{prt_name}.trk.csv"
        if "all" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
            )
        elif "rel" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_release=True,
            )
        elif "tsit" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_exit=True,
            )
        elif "tstp" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_timestep=True,
            )
        elif "wksk" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_weaksink=True,
            )
        elif "term" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_terminate=True,
            )
        elif "mult" in name:
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_release=True,
                track_terminate=True,
            )
        elif "trts" in name or "open" in name:
            tracktimes_path = prt_ws / "tracktimes.txt"
            if "open" in name:
                with open(tracktimes_path, "w") as f:
                    for t in tracktimes:
                        f.write(str(t) + "\n")
            return flopy.mf6.ModflowPrtoc(
                prt,
                pname="oc",
                track_filerecord=[prt_track_file],
                trackcsv_filerecord=[prt_track_csv_file],
                track_usertime=True,
                ntracktimes=len(tracktimes),
                tracktimes=f"open/close {tracktimes_path.name}"
                if "open" in name
                else [(t,) for t in tracktimes],
            )

    oc = get_oc()

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
    mp7_name = get_model_name(name, "mp7")
    pgs = [
        flopy.modpath.ParticleGroup(
            particlegroupname=f"group_{grp}",
            particledata=flopy.modpath.ParticleData(
                partlocs=[p[0] for p in releasepts_mp7[grp]],
                localx=[p[1] for p in releasepts_mp7[grp]],
                localy=[p[2] for p in releasepts_mp7[grp]],
                localz=[p[3] for p in releasepts_mp7[grp]],
                timeoffset=0,
                drape=0,
            ),
            filename=f"{mp7_name}_{grp}.sloc",
        )
        for grp in ["a", "b"]
    ]
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
        particlegroups=pgs,
    )

    return mp


def build_models(idx, test):
    # build gwf model
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )
    # add idomain
    gwf = gwf_sim.get_model()
    dis = gwf.get_package("DIS")
    dis.idomain = create_idomain(
        FlopyReadmeCase.nlay, FlopyReadmeCase.nrow, FlopyReadmeCase.ncol
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

    # get head, budget, and spdis results from GWF model
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
    # convert zero-based to one-based
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # check output files exist
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (mp7_ws / mp7_pathline_file).is_file()

    # check pathlines total size
    expected_len = 0
    if "all" in name:
        expected_len = len(mp7_pls)
    if "rel" in name:
        expected_len += len(releasepts_prt["a"]) + len(releasepts_prt["b"])
    if "term" in name:
        expected_len += len(releasepts_prt["a"]) + len(releasepts_prt["b"])
    if "tsit" in name:
        expected_len += len(mp7_pls) - 2 * (
            len(releasepts_prt["a"]) + len(releasepts_prt["b"])
        )
    if "tstp" in name:
        pass
    if "wksk" in name:
        pass
    if "trts" in name or "open" in name:
        expected_len += 5324
    if "mult" in name:
        expected_len += 2 * (len(releasepts_prt["a"]) + len(releasepts_prt["b"]))
    assert len(mf6_pls) == expected_len

    # make sure mf6 pathline data have correct
    #   - model index (1)
    #   - PRP index (1 or 2, depending on release point index)
    def all_equal(col, val):
        a = col.to_numpy()
        return a[0] == val and (a[0] == a).all()

    if len(mf6_pls) > 0:
        assert all_equal(mf6_pls["imdl"], 1)
        assert set(mf6_pls[mf6_pls["iprp"] == 1]["irpt"].unique()) == set(range(1, 5))
        assert set(mf6_pls[mf6_pls["iprp"] == 2]["irpt"].unique()) == set(range(1, 6))

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

    # check that cell numbers are correct
    for i, row in list(mf6_pls.iterrows()):
        # todo debug final cell number disagreement
        if row.ireason == 3:  # termination
            continue

        x, y, z, t, ilay, icell = (row.x, row.y, row.z, row.t, row.ilay, row.icell)
        k, i, j = mg.intersect(x, y, z)
        nn = mg.get_node([k, i, j]) + 1
        neighbors = mg.neighbors(nn)
        assert np.isclose(nn, icell, atol=1) or any((nn - 1) == n for n in neighbors), (
            f"nn comparison failed: expected {nn}, got {icell}"
        )
        assert ilay == (k + 1) == 1

    if "all" in name:
        # convert mf6 pathlines to mp7 format
        mf6_pls = to_mp7_pathlines(mf6_pls)

        # drop columns for which there is no direct correspondence between mf6 and mp7
        del mf6_pls["sequencenumber"]
        del mf6_pls["particleidloc"]
        del mf6_pls["xloc"]
        del mf6_pls["yloc"]
        del mf6_pls["zloc"]
        del mf6_pls["node"]
        del mp7_pls["sequencenumber"]
        del mp7_pls["particleidloc"]
        del mp7_pls["xloc"]
        del mp7_pls["yloc"]
        del mp7_pls["zloc"]
        del mp7_pls["node"]

        # sort both dataframes
        cols = ["x", "y", "z", "time"]
        mf6_pls = mf6_pls.sort_values(by=cols)
        mp7_pls = mp7_pls.sort_values(by=cols)

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

    # get head, budget, and spdis results from GWF model
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
    # convert zero-based to one-based
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # set up plots
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
            marker="o",
            markersize=2,
            linestyle="None",
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
