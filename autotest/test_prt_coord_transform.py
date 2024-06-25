"""
Tests ability to transform release point coordinates from
global, local, or local offset coordinate representations
to model coordinates, as well as to transform output data
from model to global cooordinates.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

Test cases are defined for each coordinate system option.
These are:

 - local_xy
 - local_xy_offset
 - dev_global_xy

"""

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridintersect import GridIntersect
from framework import TestFramework
from prt_test_utils import (
    DEFAULT_EXIT_SOLVE_TOL,
    FlopyReadmeCase,
    get_model_name,
)
from shapely.geometry import Point

simname = "prtcrd"
cases = [f"{simname}mdl", f"{simname}lcl", f"{simname}ofs"]  # todo global
nodes = list(range(100))


def get_start_points():
    return np.add(
        np.transpose(np.array(np.meshgrid(range(10), range(10))).reshape(2, -1)),
        0.5,
    )


def get_prp(prt):
    gi = GridIntersect(prt.modelgrid, method="vertex", rtree=True)
    if "mdl" in prt.name:
        startpts = get_start_points()
        releasepts = [
            (i, (0, *gi.intersect(Point(rpt))[0].cellids), *rpt, 0.5)
            for i, rpt in enumerate(startpts)
        ]
    elif "lcl" in prt.name:
        releasepts = [
            (nn, tuple([*prt.modelgrid.get_lrc([nn])[0]]), 0.5, 0.5, 0.5)
            for nn in nodes
        ]
    elif "ofs" in prt.name:
        releasepts = [
            (nn, tuple([*prt.modelgrid.get_lrc([nn])[0]]), 0.0, 0.0, 0.0)
            for nn in nodes
        ]
    elif "gbl" in prt.name:
        # todo global
        pass

    prp_track_file = f"{prt.name}_prt.prp.trk"
    prp_track_csv_file = f"{prt.name}_prt.prp.trk.csv"
    return flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt.name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        local_xy="lcl" in prt.name,
        local_xy_offset="ofs" in prt.name,
        exit_solve_tolerance=DEFAULT_EXIT_SOLVE_TOL,
        extend_tracking=True,
    )


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
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)

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

    # create prp package
    prp = get_prp(prt)

    # create output control package
    prt_budget_file = f"{prt_name}.bud"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        saverecord=[("BUDGET", "ALL")],
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


def build_models(idx, test):
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_sim = FlopyReadmeCase.get_gwf_sim(test.name, gwf_ws, test.targets["mf6"])
    prt_sim = build_prt_sim(test.name, gwf_ws, prt_ws, test.targets["mf6"])
    return gwf_sim, prt_sim


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwf_sim, prt_sim = test.sims[:2]
    gwf = gwf_sim.get_model()
    prt = prt_sim.get_model()
    hds = gwf.output.head().get_data()
    mg = gwf.modelgrid
    pathlines = pd.read_csv(prt_ws / f"{prt.name}.trk.csv")
    startpts = pathlines[pathlines.ireason == 0]
    endpts = pathlines[pathlines.ireason == 3]

    # check start points
    assert np.allclose(
        np.sort(startpts[["x", "y"]].to_numpy(), axis=0),
        np.sort(np.array(get_start_points()), axis=0),
    )

    # setup plot
    plot_data = False
    if plot_data:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
        for a in ax:
            a.set_aspect("equal")

        # plot pathline in map view
        pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
        pmv.plot_grid()
        pmv.plot_array(hds[0], alpha=0.1, cmap="jet")
        pmv.plot_pathline(pathlines)

        # view/save plot
        plt.show()
        plt.savefig(gwf_ws / f"test_{simname}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()
