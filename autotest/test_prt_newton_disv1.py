"""
Tests particle tracking on a vertex (DISV) grid
that reduces to a regular grid. This exercises
PRT's ability to detect when a vertex grid can
be solved via Pollock's method applied to quad-
refined cells, instead of the new ternary method
which applies more generally to polygonal cells.

The simulation uses the Newton formulation.

Several cases are provided:
    - default: No user-specified tracking times, MP7 in pathline mode.
    - bprp: Mismatching cell IDs in PRP input, expect PRT to catch and reject these.
    - trts: User-specified tracking times, some falling exactly on boundaries between
            time steps. PRT and MP7 should both assign the datum to the prior time step.

This test case also exercises the `print_input`
option which enables logging for the package's
particle release settings to the listing file.
"""

from os import environ
from platform import system

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.plot.plotutil import to_mp7_pathlines
from flopy.utils import EndpointFile, PathlineFile
from flopy.utils.binaryfile import HeadFile
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import (
    all_equal,
    check_budget_data,
    check_track_data,
    get_partdata,
    has_default_boundnames,
    plot_nodes_and_vertices,
)

simname = "prtdisv1"
cases = [f"{simname}", f"{simname}bprp", f"{simname}trts"]

# model info
nlay = 1
nrow = 10
ncol = 10
ncpl = nrow * ncol
delr = 1.0
delc = 1.0
nper = 1
perlen = 10
nstp = 5
tsmult = 1.0
tdis_rc = [(perlen, nstp, tsmult)]
top = 25.0
botm = [20.0]
strt = 20
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97
porosity = 0.1
tracktimes = list(np.linspace(0, 19, 20))


# vertex grid properties
disvkwargs = get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, botm)

# release points in mp7 format
releasepts_mp7 = [
    # node number, localx, localy, localz
    (i * 10, 0.5, 0.5, 0.5)
    for i in range(10)
]


def build_gwf_sim(idx, ws, mf6):
    gwf_name = f"{cases[idx]}_gwf"
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, version="mf6", exe_name=mf6, sim_ws=ws
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwf_name, newtonoptions="NEWTON", save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    disv = flopy.mf6.ModflowGwfdisv(gwf, **disvkwargs)

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )

    # constant head boundary
    spd = {
        0: [[(0, 0), 1.0, 1.0], [(0, 99), 0.0, 0.0]],
        # 1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=spd,
        auxiliary=["concentration"],
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_name}.cbc",
        head_filerecord=f"{gwf_name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwf_name}.oc",
    )

    # Print human-readable heads
    obs_lst = []
    for k in np.arange(0, 1, 1):
        for i in np.arange(40, 50, 1):
            obs_lst.append(["obs_" + str(i + 1), "head", (k, i)])

    obs_dict = {f"{gwf_name}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obs_dict)

    return sim


def build_prt_sim(idx, gwf_ws, prt_ws, mf6):
    # create simulation
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create prt model
    prt_name = f"{cases[idx]}_prt"
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    # create prt discretization
    disv = flopy.mf6.ModflowGwfdisv(prt, **disvkwargs)

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # convert mp7 particledata to prt release points
    partdata = get_partdata(prt.modelgrid, releasepts_mp7)
    releasepts = list(partdata.to_prp(prt.modelgrid))
    if "bprp" in name:
        # wrong cell index, point is in cell (0, 0)
        releasepts[0] = (0, (0, 1), 0.5, 9.5, 22.5)

    # create prp package
    for i in range(2):
        prp_track_file = f"{prt_name}{i}.prp.trk"
        prp_track_csv_file = f"{prt_name}{i}.prp.trk.csv"
        flopy.mf6.ModflowPrtprp(
            prt,
            pname=f"prp{i}",
            filename=f"{prt_name}{i}.prp",
            nreleasepts=len(releasepts),
            packagedata=releasepts,
            perioddata={0: ["FIRST"]},
            track_filerecord=[prp_track_file],
            trackcsv_filerecord=[prp_track_csv_file],
            stop_at_weak_sink=False,
            boundnames=True,
            print_input=True,
            dev_forceternary=i == 1,
            exit_solve_tolerance=1e-10,
            extend_tracking=True,
        )

    # create output control package
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    if "trts" in name:
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            track_filerecord=[prt_track_file],
            trackcsv_filerecord=[prt_track_csv_file],
            track_release=True,
            track_terminate=True,
            track_exit=True,
            track_usertime=True,
            ntracktimes=len(tracktimes) if "trts" in name else None,
            tracktimes=[(t,) for t in tracktimes] if "trts" in name else None,
        )
    else:
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            track_filerecord=[prt_track_file],
            trackcsv_filerecord=[prt_track_csv_file],
        )

    # create the flow model interface
    gwf_name = f"{cases[idx]}_gwf"
    gwf_budget_file = gwf_ws / f"{gwf_name}.cbc"
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


def build_mp7_sim(idx, ws, mp7, gwf):
    name = cases[idx]
    partdata = get_partdata(gwf.modelgrid, releasepts_mp7)
    mp7_name = f"{cases[idx]}_mp7"
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
        budgetfilename=f"{gwf.name}.cbc",
    )
    mpbas = flopy.modpath.Modpath7Bas(mp, porosity=porosity)
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype=(
            "combined" if ("trts" in name or "trtf" in name) else "pathline"
        ),
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        timepointdata=[20, tracktimes],
        particlegroups=[pg],
    )

    return mp


def build_models(idx, test):
    gwf_sim = build_gwf_sim(idx, test.workspace, test.targets["mf6"])
    prt_sim = build_prt_sim(
        idx, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )
    mp7_sim = build_mp7_sim(
        idx, test.workspace / "mp7", test.targets["mp7"], gwf_sim.get_model()
    )
    return gwf_sim, prt_sim, mp7_sim


def compare_output(name, mf6_pls, mp7_pls, mp7_eps):
    mf6_eps = mf6_pls[mf6_pls.ireason == 3]  # get prt endpoints
    mp7_eps = to_mp7_pathlines(mp7_eps)  # convert mp7 pathlines to mp7 format
    mf6_pls = to_mp7_pathlines(mf6_pls)  # convert mf6 pathlines to mp7 format
    mf6_eps = to_mp7_pathlines(mf6_eps)  # convert mf6 endpoints to mp7 format

    # sort both dataframes by particleid and time
    mf6_pls.sort_values(by=["particleid", "time"], inplace=True)
    mp7_pls.sort_values(by=["particleid", "time"], inplace=True)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pls["sequencenumber"]
    del mf6_pls["particleidloc"]
    del mf6_pls["xloc"]
    del mf6_pls["yloc"]
    del mf6_pls["zloc"]
    del mf6_pls["node"]  # node numbers reversed in y direction in mp7
    del mp7_pls["sequencenumber"]
    del mp7_pls["particleidloc"]
    del mp7_pls["xloc"]
    del mp7_pls["yloc"]
    del mp7_pls["zloc"]
    del mp7_pls["node"]

    if "bprp" not in name:
        # both methods should be (nearly) identical to mp7
        mf6_pls_plck = mf6_pls[mf6_pls.particlegroup == 1]
        mf6_pls_tern = mf6_pls[mf6_pls.particlegroup == 2]

        assert mf6_pls_plck.shape == mp7_pls.shape
        assert np.allclose(mf6_pls_plck, mp7_pls, atol=1e-3)

        del mf6_pls_tern["particleid"]
        del mf6_pls_tern["particlegroup"]
        del mp7_pls["particleid"]
        del mp7_pls["particlegroup"]

        if system() == "Windows" and environ.get("FC", None) in ["ifort", "ifx"]:
            # TODO resolve compiler differences produced by
            # ifort/ifx on Windows... ifort/ifx do not even
            # agree completely with default flags.
            assert abs(len(mf6_pls_tern) - len(mp7_pls)) < 5
        else:
            assert mf6_pls_tern.shape == mp7_pls.shape
            assert np.allclose(mf6_pls_tern, mp7_pls, atol=1e-3)


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = f"{name}_gwf"
    prt_name = f"{name}_prt"
    mp7_name = f"{name}_mp7"
    gwf_sim = test.sims[0]
    prt_sim = test.sims[1]
    gwf = gwf_sim.get_model(gwf_name)
    prt = prt_sim.get_model(prt_name)
    mg = gwf.modelgrid

    gwf_budget_file = f"{gwf_name}.cbc"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"
    mp7_endpoint_file = f"{mp7_name}.mpend"

    # if invalid release points, check for error message
    if "bprp" in name:
        buff = test.buffs[1]
        assert any("Error: release point" in l for l in buff)
        return

    # check mf6 output files exist
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    for i in range(2):
        prp_track_file = f"{prt_name}{i}.prp.trk"
        prp_track_csv_file = f"{prt_name}{i}.prp.trk.csv"
        assert (prt_ws / prp_track_file).is_file()
        assert (prt_ws / prp_track_csv_file).is_file()

    # check mp7 output files exist
    assert (mp7_ws / mp7_pathline_file).is_file()
    assert (mp7_ws / mp7_endpoint_file).is_file()

    # check list file for logged release configuration
    list_file = prt_ws / f"{prt_name}.lst"
    assert list_file.is_file()
    lines = open(list_file).readlines()
    lines = [l.strip() for l in lines]
    for iprp in range(1, 3):
        i = lines.index(f"PARTICLE RELEASE FOR PRP {iprp}")
        assert "RELEASE SCHEDULE:" in lines[i + 1]

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mp7 endpoint results
    epf = EndpointFile(mp7_ws / mp7_endpoint_file)
    mp7_eps = pd.DataFrame(epf.get_destination_endpoint_data(range(mg.nnodes)))
    # convert zero-based to one-based indexing in mp7 results
    mp7_eps["particlegroup"] = mp7_eps["particlegroup"] + 1
    mp7_eps["node"] = mp7_eps["node"] + 1
    mp7_eps["k"] = mp7_eps["k"] + 1

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # make sure pathline df has "name" (boundname) column and default values
    assert "name" in mf6_pls
    assert has_default_boundnames(mf6_pls)

    # make sure all mf6 pathline data have correct model index
    assert all_equal(mf6_pls["imdl"], 1)

    # check budget data were written to mf6 prt list file
    check_budget_data(prt_ws / f"{name}_prt.lst", perlen, nper, nstp)

    # check mf6 prt particle track data were written to binary/CSV files
    # and that different formats are equal
    check_track_data(
        track_bin=prt_ws / prt_track_file,
        track_hdr=str(prt_ws / prt_track_file).replace(".trk", ".trk.hdr"),
        track_csv=prt_ws / prt_track_csv_file,
    )

    # extract head, budget, and specific discharge results from GWF model
    head = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # compare prt and mp7 results
    compare_output(name, mf6_pls, mp7_pls, mp7_eps)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = f"{name}_gwf"
    prt_name = f"{name}_prt"
    mp7_name = f"{name}_mp7"
    gwf_sim = test.sims[0]
    prt_sim = test.sims[1]
    gwf = gwf_sim.get_model(gwf_name)
    prt = prt_sim.get_model(prt_name)
    mg = gwf.modelgrid

    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"
    mp7_endpoint_file = f"{mp7_name}.mpend"

    # extract head, budget, and specific discharge results from GWF model
    head = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mp7 endpoint results
    epf = EndpointFile(mp7_ws / mp7_endpoint_file)
    mp7_eps = pd.DataFrame(epf.get_destination_endpoint_data(range(mg.nnodes)))
    # convert zero-based to one-based indexing in mp7 results
    mp7_eps["particlegroup"] = mp7_eps["particlegroup"] + 1
    mp7_eps["node"] = mp7_eps["node"] + 1
    mp7_eps["k"] = mp7_eps["k"] + 1

    # setup plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(head[0], alpha=0.2)
    pmv.plot_vector(qx, qy, normalize=True, color="white")

    # plot labeled nodes and vertices
    plot_nodes_and_vertices(gwf, gwf.modelgrid, None, gwf.modelgrid.ncpl, ax[0])
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 pathlines",
            linestyle="--" if "trst" in name else None,
            marker="o",
            markersize=4,
            x="x",
            y="y",
            ax=ax[0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # plot mp7 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(head[0], alpha=0.2)
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
            marker="o",
            markersize=4,
        )

    # view/save plot
    plt.show()
    plt.savefig(prt_ws / f"{name}.png")


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if (plot and "bprp" not in name) else None,
        targets=targets,
        compare=None,
        xfail=[False, "bprp" in name, False],
    )
    test.run()
