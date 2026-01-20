"""

Model based on Beg et al. 2022.  It consists of a single
segment discretized into 21 reaches.  There is inflow into
reach 1 based on time series input and outflow to a constant
stage cell in reach 21.  The simulation is for one day using
time steps of 600 seconds (144 time steps).  The results are
compared to a mf2005 run using the SWR package.  A plot can
optionally be created, which will also show results from a
HEC-RAS simulation.

"""

import flopy
import numpy as np
import pandas as pd
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = [
    "chf-beg2022",
]
data_path = project_root_path / "autotest/data/beg2022/"


def build_models(idx, test):
    sim_ws = test.workspace
    name = "chfmodel"
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=sim_ws,
        memory_print_option="all",
    )

    hr2sec = 60.0 * 60.0
    dt = 600  # seconds
    perlen = 24 * hr2sec
    nstp = perlen / dt
    perioddata = [(0.0, 1, 1.0), (perlen, nstp, 1.0)]
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="SECONDS", nper=len(perioddata), perioddata=perioddata
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        # no_ptcrecord=True,
        print_option="all",
        # under_relaxation="DBD",
        # under_relaxation_theta=0.9,
        # under_relaxation_kappa=0.0001,
        # under_relaxation_gamma=0.0,
        # backtracking_number=200,
        # backtracking_tolerance=1.1,
        # backtracking_reduction_factor=0.2,
        # backtracking_residual_limit=1.0,
        linear_acceleration="BICGSTAB",
        outer_dvclose=1.0e-4,
        inner_dvclose=1.0e-4,
    )
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    total_length = 21000
    dx = 1000.0
    nreach = int(total_length / dx)
    vertices = []
    vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
    cell1d = []
    for j in range(nreach):
        cell1d.append([j, 0.5, 2, j, j + 1])
    nodes = len(cell1d)
    nvert = len(vertices)

    slope = 1.0 / 10000.0
    x = np.linspace(dx / 2, total_length - dx / 2, nreach)
    z = (total_length - x) * slope

    disv1d = flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=nodes,
        nvert=nvert,
        width=40.0,
        bottom=z,
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )

    dfw = flopy.mf6.ModflowChfdfw(
        chf,
        central_in_space=True,
        print_flows=True,
        save_flows=True,
        manningsn=1.0 / 80.0,
        idcxs=0,
    )

    sto = flopy.mf6.ModflowChfsto(
        chf,
        save_flows=True,
        steady_state={0: True, 1: False},
        transient={0: False, 1: True},
    )

    water_depth = 4.0
    strt = z + water_depth
    ic = flopy.mf6.ModflowChfic(chf, strt=strt)

    xfraction = np.array([0.0, 0.0, 10.0, 15.0, 25.0, 30.0, 40.0, 40.0]) / 40.0
    height = [40.0, 10.0, 10.0, 0.0, 0.0, 10.0, 10.0, 40.0]
    npts = len(height)
    mannfraction = npts * [1.0]
    cxsdata = list(zip(xfraction, height, mannfraction))
    cxs = flopy.mf6.ModflowChfcxs(
        chf,
        nsections=1,
        npoints=npts,
        packagedata=[(0, npts)],
        crosssectiondata=cxsdata,
    )

    # output control
    oc = flopy.mf6.ModflowChfoc(
        chf,
        budget_filerecord=f"{name}.bud",
        stage_filerecord=f"{name}.stage",
        saverecord=[
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("STAGE", "LAST"),
            ("BUDGET", "ALL"),
        ],
    )

    # time, reach1 (cms)
    reach_inflow = [
        (0, 20.0),
        (2 * hr2sec, 20),
        (3 * hr2sec, 25),
        (4 * hr2sec, 20),
        (24 * hr2sec, 20),
    ]
    flwlist = [
        [(0,), "reach1"],
    ]
    flw = flopy.mf6.ModflowChfflw(
        chf,
        maxbound=len(flwlist),
        print_input=True,
        print_flows=True,
        stress_period_data=flwlist,
    )
    filename = name + ".flw.ts"
    time_series_namerecord = [("reach1")]
    interpolation_methodrecord = [("linearend")]
    flw.ts.initialize(
        filename=filename,
        timeseries=reach_inflow,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
    )

    chd = flopy.mf6.ModflowChfchd(
        chf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(nreach - 1, z[-1] + water_depth)],
    )

    return sim, None


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    hecras = data_path / "hecras0125.csv.cmp"
    df_hecras = pd.read_csv(hecras, index_col=False)
    print(df_hecras)

    swrdata = data_path / "mfswr0125.csv.cmp"
    df_mfswr = pd.read_csv(swrdata, index_col=False)
    df_mfswr = df_mfswr.loc[df_mfswr["RCHGRP"] == 21]
    df_mfswr = df_mfswr.loc[df_mfswr["TOTTIME"] > 86400]
    print(df_mfswr)

    fpth = test.workspace / "chfmodel.bud"
    budobj = flopy.utils.binaryfile.CellBudgetFile(fpth, precision="double")
    flowja = budobj.get_data(text="FLOW-JA-FACE")
    qstorage = budobj.get_data(text="STORAGE")
    qflw = budobj.get_data(text="FLW")
    qchd = budobj.get_data(text="CHD")

    qoutflow = []
    times = np.array(budobj.times)
    for ra in qchd:
        q = -ra[0]["q"]
        qoutflow.append(q)

    qinflow = []
    for ra in qflw:
        q = ra[0]["q"]
        qinflow.append(q)

    # plot upstream and downstream flow
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times / 60.0 / 60.0, qinflow, "r-", label="Inflow")
    ax.plot(times / 60.0 / 60.0, df_hecras["Flow Flow (CMS)"], "b-", label="HEC-RAS")
    x = df_mfswr["TOTTIME"] - 86400.0
    x = x / 60.0 / 60.0
    ax.plot(x, -df_mfswr["QCRFLOW"], "go:", mfc="none", label="MODFLOW-SWR")
    ax.plot(times / 60.0 / 60.0, qoutflow, "bo:", mfc="none", label="MODFLOW 6")
    ax.set_xlim(0, 24.0)
    ax.set_ylim(19, 26)
    plt.xlabel("time, in hours")
    plt.ylabel("flow, in meters cubed per second")
    plt.legend()
    fname = test.workspace / "chfmodel.flow.png"
    plt.savefig(fname)

    # read and plot stages
    fpth = test.workspace / "chfmodel.stage"
    qobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage = qobj.get_alldata()

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times / 60.0 / 60.0, stage[:, 0, 0, 0], "r-", label="Upstream")
    ax.plot(
        times / 60.0 / 60.0, stage[:, 0, 0, -1], "bo", mfc="none", label="Downstream"
    )
    # ax.set_xlim(0, 24.)
    # ax.set_ylim(19, 26)
    plt.xlabel("time, in hours")
    plt.ylabel("stage, in meters")
    plt.legend()
    fname = test.workspace / "chfmodel.stage.png"
    plt.savefig(fname)

    return


def check_output(idx, test):
    print("evaluating model...")

    # get MFSimulation from test
    sim = test.sims[0]

    # assign name
    name = "chfmodel"

    # read the binary grid file
    fpth = test.workspace / f"{name}.disv1d.grb"
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    ia = grb.ia
    ja = grb.ja
    assert ia.shape[0] == grb.ncells + 1, "ia in grb file is not correct size"

    # check to make sure stage file can be read
    fpth = test.workspace / f"{name}.stage"
    qobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage = qobj.get_alldata()

    # check to make sure budget file can be read
    fpth = test.workspace / f"{name}.bud"
    budobj = flopy.utils.binaryfile.CellBudgetFile(fpth, precision="double")
    flowja = budobj.get_data(text="FLOW-JA-FACE")
    qstorage = budobj.get_data(text="STORAGE")
    qflw = budobj.get_data(text="FLW")
    qchd = budobj.get_data(text="CHD")
    qresidual = np.zeros(grb.nodes)

    # compare the mf6 swr outflow to the mf2005 swr outflow
    swrdata = data_path / "mfswr0125.csv.cmp"
    df_mfswr = pd.read_csv(swrdata, index_col=False)
    df_mfswr = df_mfswr.loc[df_mfswr["RCHGRP"] == 21]
    df_mfswr = df_mfswr.loc[df_mfswr["TOTTIME"] >= 86400]
    qoutflow_mf2005 = -df_mfswr["QCRFLOW"].to_numpy()

    # create a list of outflows from mf6 simulation
    qoutflow_mf6 = []
    for ra in qchd:
        q = -ra[0]["q"]
        qoutflow_mf6.append(q)

    # check to make sure the difference in outflow between mf6 and
    # mf2005 is less than atol
    diff = np.abs(qoutflow_mf6 - qoutflow_mf2005)
    print(diff)
    print(f"max difference is: {diff.max()}")
    assert np.allclose(qoutflow_mf6, qoutflow_mf2005, atol=0.6)

    return


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
