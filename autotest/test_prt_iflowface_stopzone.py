import os

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from flopy.utils import binaryfile as bf
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prt2358"
cases = [
    simname,  # stop_at_weak_sink=False, No istopzone, No iface or iflowface
    f"{simname}if6iff-1", # stop_at_weak_sink=False, No istopzone, iface=6, iflowface=-1
    f"{simname}isz-1iff-1", # stop_at_weak_sink=False, istopzone=-1, iflowface=-1
    f"{simname}isz1iff-1", # stop_at_weak_sink=False, istopzone=1, iflowface=-1
    f"{simname}saws", # stop_at_weak_sink=True, No istopzone, No iflowface
]
top = 20.0
botm = 0
chd_rec = [((0, 0, 0), 10.0)]


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
        nlay=1,
        nrow=1,
        ncol=10,
        delr=100,
        delc=100,
        top=top,
        botm=botm,
    )
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=10)
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

    nlay, nrow, ncol = gwf.dis.botm.array.shape
    # get external files for model bottoms
    # (so that PRT can reference those directly)
    botm_record = gwf.dis.botm.get_record()
    idomain_record = gwf.dis.idomain.get_record()
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        length_units=gwf.dis.length_units.array,
        delr=gwf.dis.delr.array,
        delc=gwf.dis.delc.array,
        top=top,
        botm=botm,
        # idomain=1,
        xorigin=gwf.dis.xorigin.array,
        yorigin=gwf.dis.yorigin.array,
    )

    if istopzone is not None:
        izone_array = np.zeros((nlay, nrow, ncol), dtype=int)
        izone_array[:, 0, 0] = istopzone
        izone_array[:, 0, 5:] = istopzone
    else:
        izone_array = None

    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=0.1, izone=izone_array)
    # Instantiate the MODFLOW 6 prt particle release point (prp) package
    # release particles at each row, column location
    prt_i, prt_j = np.indices(gwf.modelgrid.shape[1:])
    prt_start_x, prt_start_y = gwf.modelgrid.get_local_coords(
        gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters
    )

    prt_start_data = pd.DataFrame(
        {
            "irptno": np.arange(len(prt_start_x.ravel())),
            "k": 0,
            "i": prt_i.ravel(),
            "j": prt_j.ravel(),
            "xrpt": prt_start_x.ravel(),
            "yrpt": prt_start_y.ravel(),
            # start particles near top of saturation in cell with local_z=True
            "zrpt": 0.0,
            "boundname": [f"prt_{i}_{j}" for i, j in zip(prt_i.ravel(), prt_j.ravel())],
        }
    )

    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp",
        nreleasepts=len(prt_start_data),
        packagedata=prt_start_data.to_records(index=False).tolist(),
        local_z=True,
        boundnames=True,
        perioddata={0: ["FIRST"]},
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

    mp = flopy.modpath.Modpath7.create_mp7(
        modelname=mp7_name,
        trackdir="forward",
        exe_name=mp7,
        flowmodel=gwf,
        model_ws=ws,
        rowcelldivisions=1,
        columncelldivisions=1,
        layercelldivisions=1,
        nodes=np.arange(gwf.modelgrid.top.size),
    )
    # use a dedicated TDIS file
    # either for different time periods
    # or because mp7 is incompatible with some MODFLOW 6 TDIS options
    # tdis_file = f"{mp7_name}.tdis"  # f'{mp7_workspace.name}/{mp7_name}.tdis'
    # mp.tdis_file = tdis_file

    if stop_at_weak_sink:
        mp.mpsim.weaksinkoption = 2
    else:
        mp.mpsim.weaksinkoption = 1

    # set the reference time for when the particles start
    # referencetimeOption:
    # 1 = a value of time
    # 2 = stress period, time step, and relative time position within the time step
    mp.mpsim.referencetimeOption = 2
    mp.mpsim.referencetime = (0, 0, 0.0)

    # StopTimeOption for how long the particles can go
    # 1 = start or end of MODFLOW sim
    # 2 = until termination
    # 3 = specific value
    mp.mpsim.stoptimeoption = 3
    mp.mpsim.stoptime = 1e15

    # TimePointOption 1
    # TimePointCount, TimePointInterval
    mp.mpsim.timepointoption = 1
    mp.mpsim.timepointdata = [20, np.array([365.25])]

    mp.mpbas.porosity = 0.1
    # set default iface for RIV and SFR
    if iface is not None:
        mp.mpbas.defaultiface["CHD"] = iface
        mp.mpbas.defaultiface["RIV"] = iface
    else:
        mp.mpbas.defaultiface = dict()
    mp.mpbas.defaultifacecount = len(mp.mpbas.defaultiface)
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

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
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

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    import pdb; pdb.set_trace()
    
    # make a geopackage to reproduce figure
    gdf = gpd.GeoDataFrame(
        prt_results, geometry=gpd.points_from_xy(prt_results["x"], prt_results["y"])
    )
    gdf.to_file(prt_ws / f"{prt_name}.trk.csv.gpkg", index=False)
    modelgrid_gdf = gwf.modelgrid.geo_dataframe
    riv_period_array = gwf.riv.stress_period_data.get_data()[0]
    modelgrid_gdf.loc[[rec[0][2] for rec in chd_rec], "bc"] = "CHD"
    modelgrid_gdf.loc[[rec[0][2] for rec in riv_period_array], "bc"] = "RIV"
    modelgrid_gdf["head"] = headfile.get_data(totim=1)[0, 0, :]
    modelgrid_gdf.to_file(prt_ws / "grid.gpkg")

    prt_results.loc[prt_results["name"] == "PRT_0_3"]
    import pdb; pdb.set_trace()
    mp7_results.loc[mp7_results["name"] == "PRT_0_3"]


@pytest.mark.parametrize("idx, name", enumerate(cases[:1]))
def test_mf6model(idx, name, function_tmpdir, targets):
    iflowface = None
    iface = None
    istopzone = None
    stop_at_weak_sink = False
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(
            idx,
            t,
            iflowface=iflowface,
            iface=iface,
            istopzone=istopzone,
            stop_at_weak_sink=stop_at_weak_sink,
        ),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()
