"""
1D unconfined model with a constant head boundary on the left, and a pumping well,
either in the middle or the right. Modified from `test_prt_iflowface_stopzone.py`.

A single particle is placed near the left side of the grid. IFLOWFACE is assigned
to the well's right face. In the first test case, the well is internal to the grid.
In the second case, the well is at the edge of the grid.

This test compares PRT and MP7 models. MODPATH 7 ignores internal IFLOWFACE
assignments, terminating the particle on the entry face when the well is in
the interior of the grid. PRT respects internal IFLOWFACE assignments, and
tracks the particle through the internal well cell, terminating it when it
reaches the boundary face. MP7 and PRT agree when the well is at the edge
of the grid and assigned boundary face is external.
"""

import os

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.modpath import (
    ParticleData,
    ParticleGroup,
)
from flopy.plot.plotutil import to_mp7_pathlines
from flopy.utils import EndpointFile, PathlineFile
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import get_model_name

pd.set_option("display.max_columns", None)

simname = "prtiff"
cases = [f"{simname}int", f"{simname}ext"]

top = 30.0
botm = 0
nlay = 1
nrow = 1
ncol = 10
nnodes = nlay * nrow * ncol
delr = 100.0
delc = 100.0

particledata = ParticleData(
    partlocs=[(0, 0, i) for i in range(1, 2)],
    structured=True,
    particleids=list(range(1, 2)),
    localx=0.5,
    localy=0.5,
    localz=0.5,
)
porosity = 0.1


def build_gwf_sim(name, ws, mf6):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, version="mf6", sim_ws=ws)
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)]
    )
    ims = flopy.mf6.modflow.mfims.ModflowIms(sim, pname="ims", complexity="SIMPLE")
    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=10.0)
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        icelltype=1,
        k=1,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )
    wel_period_array = [((0, 0, 7 if "int" in name else 9), -50.0, 3)]
    wel_period = {0: wel_period_array}
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        save_flows=f"{gwf_name}.cbc",
        stress_period_data=wel_period,
        auxiliary=["iflowface"],
    )
    rch1 = flopy.mf6.ModflowGwfrcha(
        gwf,
        recharge=5e-4,
        auxiliary=["iflowface"],
        aux={0: [-1]},
    )
    chd_rec = [((0, 0, 0), 30.0)]
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(chd_rec),
        stress_period_data=chd_rec,
        save_flows=True,
    )
    headfile = f"{gwf_name}.hds"
    head_filerecord = [headfile]
    budgetfile = f"{gwf_name}.cbb"
    budget_filerecord = [budgetfile]
    saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
    printrecord = [("HEAD", "LAST")]
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        pname="oc",
        saverecord=saverecord,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        printrecord=printrecord,
    )
    return sim


def build_prt_sim(
    name,
    gwf,
    prt_ws,
    mf6,
    iffmeth="STOP",
    extend=False,
):
    prt_name = get_model_name(name, "prt")
    sim = flopy.mf6.MFSimulation(sim_name=prt_name, exe_name=mf6, sim_ws=prt_ws)
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=1,
        perioddata=[(1, 1, 1)],
    )
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=prt_name, model_nam_file=f"{prt_name}.nam"
    )
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)
    prpdata = list(particledata.to_prp(gwf.modelgrid, localz=True))
    prp = flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        local_z=True,
        exit_solve_tolerance=1e-5,
        extend_tracking=extend,
    )
    budgetfile_prt = f"{prt_name}.cbc"
    trackfile_prt = f"{prt_name}.trk"
    trackcsvfile_prt = f"{prt_name}.trk.csv"
    budget_record = [budgetfile_prt]
    track_record = [trackfile_prt]
    trackcsv_record = [trackcsvfile_prt]
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        saverecord=[("BUDGET", "ALL")],
        track_release=True,
        track_terminate=True,
        track_exit=True,
    )
    gwf_ws = gwf.model_ws
    rel_prt_folder = os.path.relpath(gwf_ws, start=prt_ws)
    fmi_pd = [
        ("GWFGRID", f"{rel_prt_folder}/{gwf.name}.dis.grb"),
        ("GWFHEAD", f"{rel_prt_folder}/{gwf.name}.hds"),
        ("GWFBUDGET", f"{rel_prt_folder}/{gwf.name}.cbb"),
    ]
    fmi = flopy.mf6.ModflowPrtfmi(prt, packagedata=fmi_pd)
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_mp7_sim(
    name,
    ws,
    mp7,
    gwf,
    extend=False,
):
    mp7_name = get_model_name(name, "mp7")
    pg = ParticleGroup(
        particledata=particledata,
    )
    mp = flopy.modpath.Modpath7(
        modelname=mp7_name,
        flowmodel=gwf,
        exe_name=mp7,
        model_ws=ws,
    )
    defaultiface = {"WEL": 2, "RCH": 6}
    mpbas = flopy.modpath.Modpath7Bas(
        mp,
        porosity=porosity,
        defaultiface=defaultiface,
    )
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        weaksinkoption="pass_through",
        referencetime=(0, 0, 0.0),
        stoptimeoption="extend" if extend else "total",
        particlegroups=[pg],
    )
    return mp


def build_models(
    idx,
    test,
    extend=False,
):
    gwf_sim = build_gwf_sim(
        name=test.name,
        ws=test.workspace / "gwf",
        mf6=test.targets["mf6"],
    )
    gwf = gwf_sim.get_model(get_model_name(test.name, "gwf"))
    prt_sim = build_prt_sim(
        name=test.name,
        gwf=gwf,
        prt_ws=test.workspace / "prt",
        mf6=test.targets["mf6"],
        extend=extend,
    )
    mp7_sim = build_mp7_sim(
        name=test.name,
        ws=test.workspace / "mp7",
        mp7=test.targets["mp7"],
        gwf=gwf,
        extend=extend,
    )
    return gwf_sim, prt_sim, mp7_sim


def compare_output(name, mf6_pls, mp7_pls, mp7_eps, tolerance=1e-3):
    mf6_eps = mf6_pls[(mf6_pls.ireason == 3)]
    mp7_eps = to_mp7_pathlines(mp7_eps)
    mf6_pls = to_mp7_pathlines(mf6_pls)
    mf6_eps = to_mp7_pathlines(mf6_eps)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    columns_to_remove = [
        "sequencenumber",
        "particleidloc",
        "xloc",
        "yloc",
        "zloc",
        "node",
        "stressperiod",
        "timestep",
        "xloc0",
        "yloc0",
        "zloc0",
        "node0",
        "time0",
        "zone0",
        "status",
        "initialcellface",
        "cellface",
        "zone",
        "x0",
        "y0",
        "z0",
        "k0",
    ]

    def drop_cols(df, columns):
        for col in columns:
            if col in df.columns:
                del df[col]

    drop_cols(mf6_pls, columns_to_remove)
    drop_cols(mp7_pls, columns_to_remove)
    drop_cols(mf6_eps, columns_to_remove)
    drop_cols(mp7_eps, columns_to_remove)

    assert np.isclose(mf6_eps.x, 800 if "int" in name else 1000)
    assert np.isclose(mp7_eps.x, 700 if "int" in name else 1000)
    assert np.isclose(mf6_eps.y, 50)
    assert np.isclose(mp7_eps.y, 50)
    assert np.isclose(mf6_eps.z, 12.4488 if "int" in name else 10.6019, atol=tolerance)
    assert np.isclose(mp7_eps.z, 12.4488 if "int" in name else 10.6019, atol=tolerance)


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace / "gwf"
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
    mp7_endpoint_file = f"{mp7_name}.mpend"

    headfile = HeadFile(gwf_ws / gwf_head_file)
    hds = headfile.get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    epf = EndpointFile(mp7_ws / mp7_endpoint_file)
    mp7_eps = pd.DataFrame(epf.get_destination_endpoint_data(range(mg.nnodes)))
    mp7_eps["particlegroup"] = mp7_eps["particlegroup"] + 1
    mp7_eps["node"] = mp7_eps["node"] + 1
    mp7_eps["k"] = mp7_eps["k"] + 1

    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    compare_output(
        name, mf6_pls, mp7_pls, mp7_eps, tolerance=1e-1 if "ext" in test.name else 1e-3
    )


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace / "gwf"
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
    mp7_endpoint_file = f"{mp7_name}.mpend"

    headfile = HeadFile(gwf_ws / gwf_head_file)
    hds = headfile.get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    epf = EndpointFile(mp7_ws / mp7_endpoint_file)
    mp7_eps = pd.DataFrame(epf.get_destination_endpoint_data(range(mg.nnodes)))
    mp7_eps["particlegroup"] = mp7_eps["particlegroup"] + 1
    mp7_eps["node"] = mp7_eps["node"] + 1
    mp7_eps["k"] = mp7_eps["k"] + 1

    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)
    mf6_eps = to_mp7_pathlines(mf6_pls[mf6_pls.ireason == 3])

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    fig.tight_layout(pad=3.0)
    for a in ax.ravel():
        a.set_aspect("equal")

    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax[0][0])
    pmv.plot_grid()
    pmv.plot_bc("CHD", alpha=0.4)
    pmv.plot_bc("WEL", alpha=0.4)
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6, map view",
            kind="line",
            x="x",
            y="y",
            ax=ax[0][0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
            lw=2,
        )

    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax[0][1])
    pmv.plot_grid()
    pmv.plot_bc("CHD", alpha=0.4)
    pmv.plot_bc("WEL", alpha=0.4)
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pls.groupby(["particleid"])
    for ipl, (pid, pl) in enumerate(mp7_plines):
        pl.plot(
            title="MP7, map view",
            kind="line",
            x="x",
            y="y",
            ax=ax[0][1],
            legend=False,
            color=cm.plasma(ipl / len(mp7_plines)),
            lw=2,
        )

    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax[1][0], line={"row": 0})
    pxs.plot_grid()
    pxs.plot_bc("CHD", alpha=0.4)
    pxs.plot_bc("WEL", alpha=0.4)
    pxs.plot_array(hds[0], alpha=0.1)
    pxs.plot_vector(qx, qy, qz, normalize=True, color="white")
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6, cross section",
            kind="line",
            x="x",
            y="z",
            ax=ax[1][0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
            lw=2,
        )

    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax[1][1], line={"row": 0})
    pxs.plot_grid()
    pxs.plot_bc("CHD", alpha=0.4)
    pxs.plot_bc("WEL", alpha=0.4)
    pxs.plot_array(hds[0], alpha=0.1)
    pxs.plot_vector(qx, qy, qz, normalize=True, color="white")
    for ipl, (pid, pl) in enumerate(mp7_plines):
        pl.plot(
            title="MP7, cross section",
            kind="line",
            x="x",
            y="z",
            ax=ax[1][1],
            legend=False,
            color=cm.plasma(ipl / len(mp7_plines)),
            lw=2,
        )

    plt.show()
    plt.savefig(gwf_ws / f"{name}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(
            idx,
            t,
            extend=True,
        ),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
