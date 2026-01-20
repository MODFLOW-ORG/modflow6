"""
Test coordinate checking functionality for PRT release points.

Tests the COORDINATE_CHECK_METHOD option which can be set to:
- 'eager' (default): check at release time
- 'none': skip release coordinate validation

Two test cases:
1. checks on ('eager'): Should catch mismatching points/cellids with informative error
2. no checks ('none'): Should allow mismatches and let tracking proceed. This produces
incorrect pathlines where particles jump to the specified cells after reporting initial
positions at the release coordinates.
"""

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from flopy.utils import HeadFile
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtchk"
cases = [
    f"{simname}_eager",
    f"{simname}_none",
]


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
        top=FlopyReadmeCase.top,
        botm=FlopyReadmeCase.botm,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    releasepts = [
        # particle index, k, i, j, x, y, z (wrong coordinates for cell 0,0,0)
        [0, 0, 0, 0, 5.5, 5.5, 0.5],
    ]

    # create prp package
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"

    if "eager" in name:
        coordinate_check_method = "eager"
    else:
        coordinate_check_method = "none"

    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: [("FIRST",)]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        coordinate_check_method=coordinate_check_method,
        print_input=True,
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


def build_models(test):
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets["mf6"],
    )
    return gwf_sim, prt_sim


def check_output(test):
    name = test.name

    if "eager" in name:
        buff = test.buffs[1]
        assert any("Error: release point" in l for l in buff)


def plot_output(test):
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


@pytest.mark.developmode
@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build_models,
        check=check_output,
        plot=plot_output if plot else None,
        targets=targets,
        compare=None,
        xfail=[False, "eager" in name],
    )
    test.run()
