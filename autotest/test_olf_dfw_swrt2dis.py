"""

SWR Test Problem 2 simulates two-dimensional overland flow using
a grid of rows and columns.  The SWR code was compared with results
from SWIFT2D, a USGS 2D overland flow simulator.  This version of
the problem uses the DIS Package as a 2d grid.  The problem is set up
so that once steady conditions are achieved, the depth in each reach
should be 1.0 m.

"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "olf-swrt2-dis01",
    "olf-swrt2-dis02",
]


def build_models(idx, test):
    dx = 500.0
    nper = 1
    perlen = [5040 * 2 * 60.0]  # 7 days (in seconds)
    nstp = [50]  # In SWR report nstp = [5040] and tsmult is 1.
    tsmult = [1.2]

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "olf"

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
    olfname = f"{name}_model"
    olf = flopy.mf6.ModflowOlf(sim, modelname=olfname, save_flows=True)

    nouter, ninner = 100, 50
    hclose, rclose, relax = 1e-6, 1e-8, 1.0
    imsolf = flopy.mf6.ModflowIms(
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
        # backtracking_number=5,
        # backtracking_tolerance=1.0,
        # backtracking_reduction_factor=0.3,
        # backtracking_residual_limit=100.0,
        filename=f"{olfname}.ims",
    )
    sim.register_ims_package(imsolf, [olf.name])

    nrow = 11
    ncol = 11
    bottom = np.empty((nrow, ncol), dtype=float)
    for i in range(nrow):
        bottom[i, :] = np.linspace(1.05, 0.05, nrow)

    dis = flopy.mf6.ModflowOlfdis2D(
        olf,
        nrow=nrow,
        ncol=ncol,
        delr=dx,
        delc=dx,
        bottom=bottom,
    )

    dfw = flopy.mf6.ModflowOlfdfw(
        olf,
        print_flows=True,
        save_flows=True,
        manningsn=0.30,
        idcxs=None,
    )

    sto = flopy.mf6.ModflowOlfsto(olf, save_flows=True)

    ic = flopy.mf6.ModflowOlfic(olf, strt=2.05)

    # output control
    oc = flopy.mf6.ModflowOlfoc(
        olf,
        budget_filerecord=f"{olfname}.bud",
        stage_filerecord=f"{olfname}.stage",
        saverecord=[
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    qinflow = 23.570

    if idx == 0:
        # apply inflow with flw package
        spd = [(i, 0, qinflow) for i in range(nrow)]
        flw = flopy.mf6.ModflowOlfflw(
            olf,
            maxbound=len(spd),
            print_input=True,
            print_flows=True,
            stress_period_data=spd,
        )

    if idx == 1:
        # apply inflow with pcp package
        spd = spd = [(i, 0, qinflow / dx**2) for i in range(nrow)]
        pcp = flopy.mf6.ModflowOlfpcp(
            olf,
            maxbound=len(spd),
            print_input=True,
            print_flows=True,
            stress_period_data=spd,
        )

    spd = [(i, ncol - 1, 1.05) for i in range(nrow)]
    chd = flopy.mf6.ModflowOlfchd(
        olf,
        maxbound=len(spd),
        print_input=True,
        print_flows=True,
        stress_period_data=spd,
    )

    obs_data = {
        f"{olfname}.obs.csv": [
            ("OBS1", "STAGE", (5, 1)),
            ("OBS2", "STAGE", (5, 5)),
            ("OBS3", "STAGE", (5, 8)),
        ],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        olf,
        filename=f"{olfname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    return sim, None


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    fpth = test.workspace / "olf_model.obs.csv"
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
            label=f"MF6 obs {irch}",
        )
        # ax.plot(obsvals["time"], answer[f"STAGE00000000{irch:02d}"], "k-", label=f"SWR Reach {irch}")  # noqa
    ax.set_xlim(0, 30.0)
    ax.set_ylim(1.2, 2.4)
    plt.xlabel("time, in hours")
    plt.ylabel("stage, in meters")
    plt.legend()
    fname = test.workspace / "olf_model.obs.1.png"
    plt.savefig(fname)

    return


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    olfname = "olf_model"
    ws = test.workspace
    mfsim = flopy.mf6.MFSimulation.load(sim_ws=ws)

    # read binary stage file
    fpth = test.workspace / f"{olfname}.stage"
    sobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage_all = sobj.get_alldata()

    # at end of simulation, water depth should be 1.0 for all reaches
    olf = mfsim.get_model(olfname)
    depth = stage_all[-1] - olf.dis.bottom.array
    assert np.allclose(depth, 1.0), (
        f"Simulated depth at end should be 1, but found {depth}"
    )


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
