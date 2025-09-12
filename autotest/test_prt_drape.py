"""
Tests the "drape" option for the PRP package,
which moves particles released into dry cells
to the top-most active cell below, if any.

The grid is a 10x10 square with 2 layers, based
on the flow system in autotest/test_gwf_rch01.py.

Particles are released from the top left cell.
"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import (
    all_equal,
    check_track_data,
    get_model_name,
)

simname = "prtdrape"
cases = [simname, f"{simname}_drp"]
nlay, nrow, ncol = 2, 1, 5
chdheads = [25.0]
nper = len(chdheads)
perlen = nper * [0.01]
nstp = nper * [1]
tsmult = nper * [1.0]
delr = delc = 1.0
strt = [[[25.0, 25.0, 75.0, 25.0, 25.0], [25.0, 25.0, 75.0, 25.0, 25.0]]]
strt = np.array(strt, dtype=float)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97
tdis_rc = [(perlen[i], nstp[i], tsmult[i]) for i in range(nper)]
porosity = 0.1
releasepts = [
    # particle index, k, i, j, x, y, z
    [i, 0, 0, 1, float(f"1.{i + 1}"), float(f"0.{i + 1}"), 75]
    for i in range(5)
] + [[i, 0, 0, 3, float(f"3.{i + 1}"), float(f"0.{i + 1}"), 75] for i in range(5, 9)]


def build_gwf_sim(name, ws, mf6):
    ws = Path(ws)
    gwf_name = get_model_name(name, "gwf")

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name, version="mf6", exe_name=mf6, sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # set ims csv files
    csv0 = f"{gwf_name}.outer.ims.csv"
    csv1 = f"{gwf_name}.inner.ims.csv"

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=csv0,
        csv_inner_output_filerecord=csv1,
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

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=100.0,
        botm=[50.0, 0.0],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
        icelltype=1,
        k=1.0,
    )

    sto = flopy.mf6.ModflowGwfsto(gwf, ss=1.0e-5, sy=0.1)

    # chd files
    chdspd = {}
    for kper, chdval in enumerate(chdheads):
        chdspd[kper] = [
            [(nlay - 1, 0, 0), chdval],
            [(nlay - 1, 0, ncol - 1), chdval],
        ]
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.1)

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

    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    prt_ws = Path(prt_ws)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=prt_name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create prt model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    # create prt discretization
    dis = flopy.mf6.ModflowGwfdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=100.0,
        botm=[50.0, 0.0],
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # create prp package
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        drape="drp" in name,
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
    )

    # create the flow model interface
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


def build_models(idx, test):
    gwf_sim = build_gwf_sim(test.name, test.workspace, test.targets["mf6"])
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )
    return gwf_sim, prt_sim


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid
    drape = "drp" in name

    # check mf6 output files exist
    gwf_budget_file = f"{gwf_name}.cbc"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (prt_ws / prp_track_file).is_file()
    assert (prt_ws / prp_track_csv_file).is_file()

    # make sure all mf6 pathline data have correct model and PRP index (1)
    assert all_equal(mf6_pls["imdl"], 1)
    assert all_equal(mf6_pls["iprp"], 1)

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

    if drape:
        assert mf6_pls.shape[0] == 36
    else:
        # expect no movement without drape
        assert mf6_pls.shape[0] == 9
        # istatus=8 permanently unreleased
        assert mf6_pls.istatus.eq(8).all()


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid
    drape = "drp" in name

    # check mf6 output files exist
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

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
            title=f"MF6 pathlines{' (drape)' if drape else ''}",
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
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
