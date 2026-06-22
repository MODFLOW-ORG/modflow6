"""
Test a PRT model exchange-coupled to an ATS-enabled GWF model. Verify
that PRT "stages" particle state during each attempted solve and only
"commits" it after a successful solve. Also verify that PRT re-solves
on Picard iterations, since each iteration may yield different flows.

There are four cases, covering both buffer strategies (memory and scratch
file) and both solver configurations:
  mxiter=1: ATS retry only (no Picard loop)
  mxiter=2: Picard loop enabled

The GWF model is nonlinear (unconfined, icelltype=1), so with mxiter=2
the GWF results are different on the 2nd iteration, and PRT also gives
different particle trajectories.
"""

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from flopy.utils.binaryfile import CellBudgetFile, HeadFile
from framework import TestFramework
from prt_test_utils import get_model_name
from test_gwf_ats01 import build_gwf_sim

simname = "prtatsexg"

cases = {
    "mxiter=1, memory buffer": (simname, 1, False),
    "mxiter=2, memory buffer": (simname + "pic", 2, False),
    "mxiter=1, scratch buffer": (simname, 1, True),
    "mxiter=2, scratch buffer": (simname + "pic", 2, True),
}


def build_mf6_sim(name, ws, mf6, mxiter=1, scratch_buffer=False):
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    sim = build_gwf_sim(name, ws, mf6)

    # create prt model
    gwf = sim.get_model()
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)
    grid = gwf.modelgrid

    # create prt discretization
    flopy.mf6.modflow.mfprtdis.ModflowPrtdis(
        prt,
        pname="dis",
        nlay=grid.nlay,
        nrow=grid.nrow,
        ncol=grid.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=0.1)

    # create prp package
    rpts = [
        # particle index, k, i, j, x, y, z
        [i, 0, 0, 0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5]
        for i in range(9)
    ]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(rpts),
        packagedata=rpts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )

    # create output control package
    prt_budget_file = f"{prt_name}.cbb"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        scratch_buffer=scratch_buffer,
        printrecord=[("BUDGET", "ALL")],
        saverecord=[("BUDGET", "ALL")],
    )

    # create exchange
    gwf_name = get_model_name(name, "gwf")
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

    if mxiter > 1:
        sim.name_file.mxiter.set_data(mxiter, key=0)

    return sim


def build_models(idx, test, mxiter, scratch):
    return build_mf6_sim(
        test.name,
        test.workspace,
        test.targets["mf6"],
        mxiter=mxiter,
        scratch_buffer=scratch,
    )


def check_output(idx, test, snapshot):
    name = test.name
    ws = test.workspace
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    gwf_head_file = ws / f"{gwf_name}.hds"
    gwf_budget_file = ws / f"{gwf_name}.cbc"

    gwf_hds = HeadFile(gwf_head_file)
    gwf_cbb = CellBudgetFile(gwf_budget_file)
    gwf_times = set(gwf_hds.get_times())
    assert gwf_times == set(gwf_cbb.get_times())

    prt_track_csv_file = ws / f"{prt_name}.trk.csv"
    prt_budget_file = ws / f"{prt_name}.cbb"

    prt_pls = pd.read_csv(prt_track_csv_file)
    prt_cbb = CellBudgetFile(prt_budget_file)
    assert gwf_times == set(prt_cbb.get_times())
    assert (prt_pls["ireason"] == 0).sum() == 9

    # Compare particle tracks to snapshot. mxiter=1 and mxiter=2 produce
    # different trajectories because the GWF model is nonlinear (unconfined,
    # icelltype=1): with mxiter=2 PRT re-tracks on each Picard iteration with
    # updated GWF flows, so the final paths reflect more converged flow results.
    actual = prt_pls.drop("name", axis=1).round(3)
    assert snapshot == actual.to_records(index=False)


def plot_output(idx, test):
    name = test.name
    ws = test.workspace
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # check mf6 output files exist
    gwf_head_file = ws / f"{gwf_name}.hds"
    prt_track_csv_file = ws / f"{prt_name}.trk.csv"
    prt_budget_file = ws / f"{prt_name}.cbb"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    prt_pls = pd.read_csv(ws / prt_track_csv_file, na_filter=False)
    prt_bud = flopy.utils.CellBudgetFile(prt_budget_file, precision="double")

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = prt_pls.groupby(["iprp", "irpt", "trelease"])
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

    # view/save plot
    plt.show()
    plt.savefig(ws / f"{name}.png")


@pytest.mark.snapshot
@pytest.mark.parametrize("idx, case", enumerate(cases.values()), ids=cases.keys())
def test_mf6model(idx, case, function_tmpdir, targets, array_snapshot, plot):
    name, mxiter, scratch = case
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t, mxiter, scratch),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
