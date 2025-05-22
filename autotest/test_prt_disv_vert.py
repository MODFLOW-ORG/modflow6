"""
This reproduces a ternary method regression
where the particle's local z coordinate was
improperly set to the bottom of the cell in
an attempt to clamp the z coordinate to the
unit interval.

This is a DISV grid which reduces to a DIS
grid; we run the PRT model twice, once with
pollock's method and once with the ternary
method, and check that the results are equal.
"""

from os import environ

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.binaryfile import HeadFile
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from modflow_devtools.misc import is_in_ci
from prt_test_utils import get_model_name, get_partdata

simname = "prtdisvvert"
cases = [simname]


# model info
nlay = 2
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
botm = [20.0, 15.0]
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


def build_gwf_sim(name, ws, mf6):
    gwf_name = f"{name}_gwf"
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


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create prt model
    prt_name = f"{name}_prt"
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    # create prt discretization
    disv = flopy.mf6.ModflowGwfdisv(prt, **disvkwargs)

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # convert mp7 particledata to prt release points
    partdata = get_partdata(prt.modelgrid, releasepts_mp7)
    releasepts = list(partdata.to_prp(prt.modelgrid))

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
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
    gwf_name = f"{name}_gwf"
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
    gwf_sim = build_gwf_sim(test.name, test.workspace / "gwf", test.targets["mf6"])
    prt_sim = build_prt_sim(
        test.name, test.workspace / "gwf", test.workspace / "prt", test.targets["mf6"]
    )
    return gwf_sim, prt_sim


def check_output(idx, test, snapshot):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    if is_in_ci() and "gfortran" not in environ.get("FC", ""):
        return

    assert snapshot == mf6_pls.drop("name", axis=1).round(2).to_records(index=False)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

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
    plt.savefig(prt_ws / f"{name}.png")

    import pyvista as pv

    pv.set_plot_theme("document")
    pv.global_theme.allow_empty_mesh = True

    from flopy.export.vtk import Vtk

    axes = pv.Axes(show_actor=False, actor_scale=2.0, line_width=5)
    vtk = Vtk(model=gwf, binary=False, smooth=False)
    vtk.add_model(gwf)
    gwf_mesh = vtk.to_pyvista()

    p = pv.Plotter(
        window_size=[700, 700],
    )
    p.enable_anti_aliasing()
    p.add_mesh(gwf_mesh, opacity=0.1, style="wireframe")
    p.add_mesh(
        mf6_pls[mf6_pls.iprp == 1][["x", "y", "z"]].to_numpy(),
        color="red",
        label="pollock's method",
        point_size=15,
    )
    p.add_mesh(
        mf6_pls[mf6_pls.iprp == 2][["x", "y", "z"]].to_numpy(),
        color="green",
        label="ternary method",
        point_size=15,
    )
    p.show()


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, array_snapshot, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
