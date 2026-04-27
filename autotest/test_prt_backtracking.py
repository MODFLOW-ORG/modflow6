"""
Test backward particle tracking with PRT-FMI's BACKWARD option.
Use the Flopy README example, but start particles in the bottom
right corner. They should terminate in the top left. Test cases
include one- and multiple-period simulations, both steady state,
as well as a simulation with multiple time steps per PRT period.
"""

import flopy
import flopy.utils.binaryfile as bf
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import (
    FlopyReadmeCase,
    get_model_name,
)

simname = "prtbktr"
cases = [
    simname + "_1per",  # 1 period, 1 step
    simname + "_2per",  # 2 periods, 1 step each
    simname + "_2stp",  # 2 periods, GWF 1 step, PRT 2 steps
]


def build_gwf_sim(name, ws, mf6, nper=1):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[
            (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult)
        ]
        * nper,
    )

    gwfname = f"{name}_gwf"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
        top=FlopyReadmeCase.top,
        botm=FlopyReadmeCase.botm,
    )

    flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

    flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
    )

    flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data={0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]]},
        auxiliary=["concentration"],
        save_flows=True,
    )

    gwf_budget_file = f"{gwfname}.bud"
    gwf_head_file = f"{gwfname}.hds"
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_budget_file,
        head_filerecord=gwf_head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    flopy.mf6.ModflowIms(sim)

    return sim


def build_models(idx, test):
    name = test.name

    mf6 = test.targets["mf6"]
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"

    gwf_name = get_model_name(name, "gwf")

    if idx == 0:
        prt_nper, prt_nstp = 1, 1
        gwf_sim = FlopyReadmeCase.get_gwf_sim(name, gwf_ws, mf6)
    elif idx == 1:
        prt_nper, prt_nstp = 2, 1
        gwf_sim = build_gwf_sim(name, gwf_ws, mf6, nper=2)
    else:
        prt_nper, prt_nstp = 2, 2
        gwf_sim = build_gwf_sim(name, gwf_ws, mf6, nper=2)

    head_file = gwf_ws / f"{gwf_name}.hds"
    budget_file = gwf_ws / f"{gwf_name}.bud"

    prt_sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        prt_sim,
        pname="tdis",
        time_units="DAYS",
        nper=prt_nper,
        perioddata=[
            (
                FlopyReadmeCase.perlen,
                prt_nstp,
                FlopyReadmeCase.tsmult,
            )
        ]
        * prt_nper,
    )

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(prt_sim, modelname=prt_name, save_flows=True)

    flopy.mf6.modflow.mfprtdis.ModflowPrtdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )

    flopy.mf6.ModflowPrtmip(
        prt,
        pname="mip",
        porosity=FlopyReadmeCase.porosity,
    )

    pd = [
        ("GWFHEAD", head_file),
        ("GWFBUDGET", budget_file),
    ]
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=pd,
        backward=True,
    )

    releasepts = [
        # particle index, k, i, j, x, y, z
        [i, 0, 9, 9, float(f"9.{i + 1}"), float(f"0.{i + 1}"), 0.5]
        for i in range(9)
    ]

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
        stop_at_weak_sink=False,
        exit_solve_tolerance=1e-5,
        extend_tracking=True,
    )

    prt_budgetfile = f"{prt_name}.bud"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budgetfile],
        track_filerecord=[f"{prt_name}.trk"],
        trackcsv_filerecord=[f"{prt_name}.trk.csv"],
        printrecord=[("BUDGET", "ALL")],
        saverecord=[("BUDGET", "ALL")],
        dev_dump_event_trace=True,
    )

    ems = flopy.mf6.ModflowEms(
        prt_sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    prt_sim.register_solution_package(ems, [prt.name])

    return gwf_sim, prt_sim


def check_output(idx, test, snapshot):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    prt_sim = test.sims[1]
    gwf = gwf_sim.get_model(gwf_name)
    prt = prt_sim.get_model(prt_name)
    mg = prt.modelgrid

    prt_track_csv_file = prt_ws / f"{prt_name}.prp.trk.csv"
    assert prt_track_csv_file.exists()

    pls = pd.read_csv(prt_track_csv_file)
    assert len(pls) > 0

    assert snapshot == pls.drop("name", axis=1).round(1).to_records(index=False)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    prt_sim = test.sims[1]
    gwf = gwf_sim.get_model(gwf_name)
    prt = prt_sim.get_model(prt_name)
    mg = prt.modelgrid

    # check mf6 output files exist
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    # extract head, budget, and specific discharge results from GWF model
    hds = bf.HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    prt_budget_file = prt_ws / f"{prt_name}.bud"
    prt_bud = flopy.utils.CellBudgetFile(prt_budget_file, precision="double")

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

    plt.show()


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot, snapshot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
