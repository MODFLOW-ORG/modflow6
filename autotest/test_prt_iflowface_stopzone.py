import os

import flopy
import geopandas as gpd
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

simname = "prt2358"
cases = {
    f"{simname}a": {
        "stop_at_weak_sink": False,
        "istopzone": None,
        "iflowface": None,
        "iface": None,
    },
    f"{simname}b": {
        "stop_at_weak_sink": False,
        "istopzone": None,
        "iflowface": -1,
        "iface": 6,
    },
    f"{simname}c": {
        "stop_at_weak_sink": False,
        "istopzone": -1,
        "iflowface": -1,
        "iface": None,
    },
    f"{simname}d": {
        "stop_at_weak_sink": False,
        "istopzone": 1,
        "iflowface": -1,
        "iface": None,
    },
    f"{simname}e": {
        "stop_at_weak_sink": True,
        "istopzone": None,
        "iflowface": None,
        "iface": None,
    },
}

# grid info
top = 20.0
botm = 0
nlay = 1
nrow = 1
ncol = 10
nnodes = nlay * nrow * ncol
delr = 100.0
delc = 100.0

# chd info
chd_rec = [((0, 0, 0), 10.0)]

# particle tracking info
particledata = ParticleData(
    partlocs=[(0, 0, i) for i in range(nnodes)],
    structured=True,
    particleids=list(range(nnodes)),
    localx=0.5,
    localy=0.5,
    localz=0.5,
)
porosity = 0.1


def build_gwf_sim(name, ws, mf6, iflowface=None):
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
    # Create the node property flow package
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        icelltype=1,
        k=1,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )

    if iflowface is not None:
        riv_period_array = [
            ((0, 0, 5), 9.99, 1e6, 8.99, iflowface),
            ((0, 0, 6), 9.99, 1e6, 8.98, iflowface),
            ((0, 0, 7), 9.99, 1e6, 8.97, iflowface),
            ((0, 0, 8), 9.99, 1e6, 8.96, iflowface),
            ((0, 0, 9), 9.99, 1e6, 8.95, iflowface),
        ]
        auxiliary = ["iflowface"]
    else:
        riv_period_array = [
            ((0, 0, 5), 9.99, 1e6, 8.99),
            ((0, 0, 6), 9.99, 1e6, 8.98),
            ((0, 0, 7), 9.99, 1e6, 8.97),
            ((0, 0, 8), 9.99, 1e6, 8.96),
            ((0, 0, 9), 9.99, 1e6, 8.95),
        ]
        auxiliary = None
    riv_period = {0: riv_period_array}
    riv = flopy.mf6.ModflowGwfriv(
        gwf,
        save_flows=f"{gwf_name}.cbc",
        stress_period_data=riv_period,
        auxiliary=auxiliary,
    )
    rch1 = flopy.mf6.ModflowGwfrcha(gwf, recharge=5e-4)

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
    name, gwf, prt_ws, mf6, iflowface=None, istopzone=None, stop_at_weak_sink=False
):
    prt_name = get_model_name(name, "prt")
    sim = flopy.mf6.MFSimulation(sim_name=prt_name, exe_name=mf6, sim_ws=prt_ws)
    # Instantiate the MODFLOW 6 temporal discretization package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=1,
        perioddata=[(1, 1, 1)],
    )
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=prt_name, model_nam_file=f"{prt_name}.nam"
    )

    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
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

    if istopzone is not None:
        izone_array = np.zeros((nlay, nrow, ncol), dtype=int)
        izone_array[:, 0, 0] = istopzone
        izone_array[:, 0, 5:] = istopzone
    else:
        izone_array = None

    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity, izone=izone_array)

    prpdata = list(particledata.to_prp(gwf.modelgrid, localz=True))  # [-1:]
    # prpdata[0] = [0, *prpdata[0][1:]]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        local_z=True,
        exit_solve_tolerance=1e-5,
        stop_at_weak_sink=stop_at_weak_sink,
        extend_tracking=True,
        istopzone=istopzone,
    )
    # Instantiate the MODFLOW 6 prt output control package
    budgetfile_prt = f"{prt_name}.cbc"
    trackfile_prt = f"{prt_name}.trk"
    trackcsvfile_prt = f"{prt_name}.trk.csv"
    budget_record = [budgetfile_prt]
    track_record = [trackfile_prt]
    trackcsv_record = [trackcsvfile_prt]
    # track positions every year for 100 years
    track_nyears = 100
    tracktimes = np.linspace(0, track_nyears * 365.25, track_nyears + 1)
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        ntracktimes=0,  # len(tracktimes),
        tracktimes=None,  # [(t,) for t in tracktimes],
        saverecord=[("BUDGET", "ALL")],
    )

    gwf_ws = gwf.model_ws
    rel_prt_folder = os.path.relpath(gwf_ws, start=prt_ws)

    # Instantiate the MODFLOW 6 prt flow model interface
    fmi_pd = [
        ("GWFHEAD", f"{rel_prt_folder}/{gwf.name}.hds"),
        ("GWFBUDGET", f"{rel_prt_folder}/{gwf.name}.cbb"),
    ]
    flopy.mf6.ModflowPrtfmi(prt, packagedata=fmi_pd)

    # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_mp7_sim(name, ws, mp7, gwf, iface=None, stop_at_weak_sink=False):
    # make an equivalent MP7 simulation
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
    defaultiface = dict() if iface is None else {"CHD": iface, "RIV": iface}
    mpbas = flopy.modpath.Modpath7Bas(
        mp,
        porosity=porosity,
        defaultiface=defaultiface,
    )
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="combined",
        trackingdirection="forward",
        budgetoutputoption="summary",
        weaksinkoption="stop_at" if stop_at_weak_sink else "pass_through",
        referencetime=(0, 0, 0.0),
        stoptimeoption="specified",
        stoptime=1e15,
        timepointdata=[20, np.array([365.25])],
        particlegroups=[pg],
    )
    return mp


def build_models(
    idx, test, iflowface=None, iface=None, istopzone=None, stop_at_weak_sink=False
):
    gwf_sim = build_gwf_sim(
        name=test.name,
        ws=test.workspace / "gwf",
        mf6=test.targets["mf6"],
        iflowface=iflowface,
    )
    gwf = gwf_sim.get_model(get_model_name(test.name, "gwf"))
    prt_sim = build_prt_sim(
        name=test.name,
        gwf=gwf,
        prt_ws=test.workspace / "prt",
        mf6=test.targets["mf6"],
        iflowface=iflowface,
        istopzone=istopzone,
        stop_at_weak_sink=stop_at_weak_sink,
    )
    mp7_sim = build_mp7_sim(
        name=test.name,
        ws=test.workspace / "mp7",
        mp7=test.targets["mp7"],
        gwf=gwf,
        iface=iface,
        stop_at_weak_sink=stop_at_weak_sink,
    )
    return gwf_sim, prt_sim, mp7_sim


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

    # check mf6 output files exist
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"
    mp7_endpoint_file = f"{mp7_name}.mpend"

    # extract head, budget, and specific discharge results from GWF model
    headfile = HeadFile(gwf_ws / gwf_head_file)
    hds = headfile.get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

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
    mf6_eps = to_mp7_pathlines(mf6_pls[mf6_pls.ireason == 3])

    # make a geopackage to reproduce figure
    gdf = gpd.GeoDataFrame(
        mf6_pls, geometry=gpd.points_from_xy(mf6_pls["x"], mf6_pls["y"])
    )
    gdf.to_file(prt_ws / f"{prt_name}.trk.csv.gpkg", index=False)
    modelgrid_gdf = gwf.modelgrid.geo_dataframe
    riv_period_array = gwf.riv.stress_period_data.get_data()[0]
    modelgrid_gdf.loc[[rec[0][2] for rec in chd_rec], "bc"] = "CHD"
    modelgrid_gdf.loc[[rec[0][2] for rec in riv_period_array], "bc"] = "RIV"
    modelgrid_gdf["head"] = headfile.get_data(totim=1)[0, 0, :]
    modelgrid_gdf.to_file(gwf_ws / "grid.gpkg")


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

    # check mf6 output files exist
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"
    mp7_endpoint_file = f"{mp7_name}.mpend"

    # extract head, budget, and specific discharge results from GWF model
    headfile = HeadFile(gwf_ws / gwf_head_file)
    hds = headfile.get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

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
    mf6_eps = to_mp7_pathlines(mf6_pls[mf6_pls.ireason == 3])

    # setup plot
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    fig.tight_layout(pad=3.0)
    for a in ax.ravel():
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0][0])
    pmv.plot_grid()
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

    # plot mp7 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0][1])
    pmv.plot_grid()
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

    # plot mf6 pathlines in cross section
    pxs = flopy.plot.PlotCrossSection(modelgrid=mg, ax=ax[1][0], line={"row": 0})
    pxs.plot_grid()
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

    # plot mp7 pathlines in cross section
    pxs = flopy.plot.PlotCrossSection(modelgrid=mg, ax=ax[1][1], line={"row": 0})
    pxs.plot_grid()
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

    # view/save plot
    plt.show()
    plt.savefig(gwf_ws / f"{name}.png")


@pytest.mark.parametrize("idx, name", enumerate(list(cases.keys())[:1]))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    case = cases[name]
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(
            idx,
            t,
            iflowface=case["iflowface"],
            iface=case["iface"],
            istopzone=case["istopzone"],
            stop_at_weak_sink=case["stop_at_weak_sink"],
        ),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
