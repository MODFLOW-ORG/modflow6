"""

SWR Test Problem 2 simulates two-dimensional overland flow using
a grid of rows and columns.  The SWR code was compared with results
from SWIFT2D, a USGS 2D overland flow simulator.  The CHF model is used
here to simulate a one-dimensional version of this same problem.  The problem
is set up so that once steady conditions are achieved, the depth in each reach
should be 1.0 m.

"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "chf-swrt2-01",
    "chf-swrt2-02",
]


def build_models(idx, test):
    dx = 500.0
    width = dx
    nreach = 11
    nper = 1
    perlen = [5040 * 2 * 60.0]  # 7 days (in seconds)
    nstp = [50]  # In SWR report nstp = [5040] and tsmult is 1.
    tsmult = [1.2]

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "chf"

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=f"{name}_sim", version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="SECONDS", nper=nper, perioddata=tdis_rc
    )

    # surface water model
    chfname = f"{name}_model"
    chf = flopy.mf6.ModflowChf(sim, modelname=chfname, save_flows=True)

    nouter, ninner = 100, 50
    hclose, rclose, relax = 1e-8, 1e-8, 1.0
    imschf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.9,
        under_relaxation_kappa=0.0001,
        under_relaxation_gamma=0.0,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        backtracking_number=5,
        backtracking_tolerance=1.0,
        backtracking_reduction_factor=0.3,
        backtracking_residual_limit=100.0,
        filename=f"{chfname}.ims",
    )
    sim.register_ims_package(imschf, [chf.name])

    vertices = []
    vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
    cell1d = []
    for j in range(nreach):
        cell1d.append([j, 0.5, 2, j, j + 1])
    nodes = len(cell1d)
    nvert = len(vertices)

    reach_bottom = np.linspace(1.05, 0.05, nreach)

    disv1d = flopy.mf6.ModflowChfdisv1D(
        chf,
        export_array_ascii=True,
        nodes=nodes,
        nvert=nvert,
        width=width,
        bottom=reach_bottom,
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )

    dfw = flopy.mf6.ModflowChfdfw(
        chf,
        export_array_ascii=True,
        print_flows=True,
        save_flows=True,
        manningsn=0.30,
        idcxs=None,
    )

    sto = flopy.mf6.ModflowChfsto(chf, save_flows=True)

    ic = flopy.mf6.ModflowChfic(
        chf,
        export_array_ascii=True,
        strt=2.05,
    )

    # output control
    oc = flopy.mf6.ModflowChfoc(
        chf,
        budget_filerecord=f"{chfname}.bud",
        stage_filerecord=f"{chfname}.stage",
        saverecord=[
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    # flw
    inflow_reach = 0
    qinflow = 23.570
    if idx == 0:
        flw = flopy.mf6.ModflowChfflw(
            chf,
            maxbound=1,
            print_input=True,
            print_flows=True,
            stress_period_data=[(inflow_reach, qinflow)],
        )

    if idx == 1:
        pcp = flopy.mf6.ModflowOlfpcp(
            chf,
            maxbound=1,
            print_input=True,
            print_flows=True,
            stress_period_data=[(inflow_reach, qinflow / dx / width)],
        )

    chd = flopy.mf6.ModflowChfchd(
        chf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(nreach - 1, 1.05)],
    )

    obs_data = {
        f"{chfname}.obs.csv": [
            ("OBS1", "STAGE", (1,)),
            ("OBS2", "STAGE", (5,)),
            ("OBS3", "STAGE", (8,)),
        ],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        chf,
        filename=f"{chfname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    return sim, None


def plot_output(idx, test):
    print("making plots...")
    import matplotlib.pyplot as plt

    fpth = test.workspace / "chf_model.obs.csv"
    obsvals = np.genfromtxt(fpth, names=True, delimiter=",")

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    for irch in [1, 2, 3]:
        ax.plot(
            obsvals["time"] / 3600.0,
            obsvals[f"OBS{irch}"],
            marker="o",
            mfc="none",
            mec="k",
            lw=0.0,
            label=f"MF6 reach {irch}",
        )
        # ax.plot(obsvals["time"], answer[f"STAGE00000000{irch:02d}"], "k-", label=f"SWR Reach {irch}")  # noqa
    ax.set_xlim(0, 30.0)
    ax.set_ylim(1.2, 2.4)
    plt.xlabel("time, in hours")
    plt.ylabel("stage, in meters")
    plt.legend()
    fname = test.workspace / "chf_model.obs.1.png"
    plt.savefig(fname)

    return


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    chfname = "chf_model"
    ws = test.workspace
    mfsim = flopy.mf6.MFSimulation.load(sim_ws=ws)

    # read binary stage file
    fpth = test.workspace / f"{chfname}.stage"
    sobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage_all = sobj.get_alldata()

    # at end of simulation, water depth should be 1.0 for all reaches
    chf = mfsim.get_model(chfname)
    depth = stage_all[-1] - chf.disv1d.bottom.array
    (
        np.allclose(depth, 1.0),
        f"Simulated depth at end should be 1, but found {depth}",
    )

    # ensure export array is working properly
    flist = [
        "disv1d.width",
        "disv1d.bottom",
        "disv1d.idomain",
        "dfw.manningsn",
        "ic.strt",
    ]
    files = [pl.Path(ws / f"{chfname}-{f}.txt") for f in flist]
    chf = test.sims[0].chf[0]
    for i, fpth in enumerate(files):
        assert fpth.is_file(), f"Expected file does not exist: {fpth.name}"
        a = np.loadtxt(fpth)
        array_name = flist[i][flist[i].index(".") + 1 :]
        package_name = flist[i][0 : flist[i].index(".")]
        package = getattr(chf, package_name)
        b = getattr(package, array_name).array
        assert np.allclose(a, b)
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
