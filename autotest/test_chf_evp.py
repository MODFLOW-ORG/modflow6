"""

Test the evaporation package for the channel model.
This example is tough to converge, because evaporation
removes all of the water in the cell.  The model did
not converge well until newton unrelaxation was added
and delta-bar-delta underrelaxation was also activated.

"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "chf-evp01",
]

# grid size
dx = 1.0
nreach = 2
total_length = dx * nreach
vertices = []
vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
cell1d = []
for j in range(nreach):
    cell1d.append([j, 0.5, 2, j, j + 1])
nodes = len(cell1d)
nvert = len(vertices)
xorigin = 0.0
yorigin = 0.0
angrot = 0.0

idomain = np.ones((nodes,), dtype=int)


def build_models(idx, test):
    perlen = [1]  # 1 second
    nstp = [10]
    tsmult = [1.0]
    nper = len(perlen)

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

    nouter, ninner = 100, 50
    hclose, rclose, relax = 1e-6, 1e-8, 1.0
    ims = flopy.mf6.ModflowIms(
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
    )

    add_chf_model_disv1d(sim)

    return sim, None


def add_chf_model_disv1d(sim):
    name = "channel"
    chf = flopy.mf6.ModflowChf(
        sim, modelname=name, save_flows=True, newtonoptions="UNDER_RELAXATION"
    )

    disv1d = flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=nodes,
        nvert=nvert,
        width=1.0,
        bottom=0.0,
        idomain=idomain,
        vertices=vertices,
        cell1d=cell1d,
        xorigin=xorigin,
        yorigin=yorigin,
        angrot=angrot,
    )

    dfw = flopy.mf6.ModflowChfdfw(
        chf,
        print_flows=True,
        save_flows=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        manningsn=0.035,
        idcxs=0,
    )

    xfraction = np.array([0.0, 1.0, 2.0])
    height = [1.0, 0.0, 1.0]
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

    ic = flopy.mf6.ModflowChfic(chf, strt=1.0)

    sto = flopy.mf6.ModflowChfsto(
        chf,
        steady_state={0: False},
        transient={0: True},
    )

    evp_rate = 1.0
    evp_spd = [
        (0, evp_rate),
        (1, evp_rate),
    ]

    pcp = flopy.mf6.ModflowChfevp(
        chf, maxbound=1, print_input=True, print_flows=True, stress_period_data=evp_spd
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
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    return


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    mfsim = test.sims[0]
    chf = mfsim.chf[0]

    sobj = chf.output.stage()
    times = np.array(sobj.times)
    ntimes = times.shape[0]
    nreach = 2
    stage = sobj.get_alldata().reshape((ntimes, nreach))

    bobj = chf.output.budget()
    evap_list = bobj.get_data(text="EVP")
    evap = []
    for rec in evap_list:
        evap.append(rec[0]["q"])
    evap = np.array(evap)

    fig, ax = plt.subplots(1, 1)
    ax.plot(times, stage[:, 0], label="stage")
    ax.plot(times, -evap, label="Qevap")
    plt.legend()

    fname = test.workspace / "results.png"
    plt.savefig(fname)


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    mfsim = test.sims[0]
    chf = mfsim.chf[0]

    # read binary stage file
    sobj = chf.output.stage()
    times = np.array(sobj.times)
    ntimes = times.shape[0]
    nreach = 2
    stage = sobj.get_alldata().reshape((ntimes, nreach))
    answer = np.linspace(0.9, 0, ntimes)
    assert np.allclose(stage[:, 0], answer, atol=1.0e-5), "stage is not correct"


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
