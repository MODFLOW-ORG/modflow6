"""

Test problem designed specifically for the DFW Package.

The following shows the stream bottom elevations for the
9-cell diffusive wave model.

    _       _      __ bottom elevation = 2.
  _| |_   _| |_    __ bottom elevation = 1.
_|     |_|     |_  __ bottom elevation = 0.
1 2 3 4 5 6 7 8 9

"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "chf-bowl01",
    "chf-bowl02",
]

inflow_location = ["middle", "left"]

strt = [0.0, 1, 2, 1, 0, 1, 2, 1, 0.0]  # start with dry bowl
# strt = [0.0, 1, 2, 2, 2, 2, 2, 1, 0.] # start bowl filled
reach_bottom = [0.0, 1, 2, 1, 0, 1, 2, 1, 0.0]


def build_models(idx, test):
    dx = 100.0
    nreach = 9
    nper = 1
    perlen = [86400.0 * 10]  # 10 days
    nstp = [100]
    tsmult = [1.0]

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
    chf = flopy.mf6.ModflowChf(
        sim,
        modelname=chfname,
        save_flows=True,
        newtonoptions="newton underrelaxation",
    )

    nouter, ninner = 200, 50
    hclose, rclose, relax = 1e-8, 1e-8, 1.0
    imschf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.95,
        under_relaxation_kappa=0.0001,
        under_relaxation_gamma=0.0,
        under_relaxation_momentum=0.0,
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
        csv_outer_output_filerecord=f"{chfname}.ims.outer.csv",
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

    disv1d = flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=nodes,
        nvert=nvert,
        width=1.0,
        bottom=reach_bottom,
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )

    dfw = flopy.mf6.ModflowChfdfw(
        chf,
        print_flows=True,
        save_flows=True,
        manningsn=0.035,
        idcxs=0,
    )

    # note: for specifying zero-based reach number, put reach number in tuple
    fname = f"{chfname}.zdg.obs.csv"
    zdg_obs = {
        fname: [
            ("OUTFLOW1", "ZDG", (0,)),
            ("OUTFLOW9", "ZDG", (nodes - 1,)),
        ],
        "digits": 10,
    }

    idcxs = 0  # use cross section 0
    width = 1.0
    slope = 1.0 / dx
    rough = 0.035
    spd = [((nreach - 1,), idcxs, width, slope, rough)]
    if inflow_location[idx] == "middle":
        spd.append((0, idcxs, width, slope, rough))
    zdg = flopy.mf6.ModflowChfzdg(
        chf,
        observations=zdg_obs,
        print_input=True,
        maxbound=len(spd),
        stress_period_data=spd,
    )

    sto = flopy.mf6.ModflowChfsto(chf, save_flows=True)

    ic = flopy.mf6.ModflowChfic(chf, strt=strt)

    xfraction = [0.0, 0.0, 1.0, 1.0]
    height = [100.0, 0.0, 0.0, 100.0]
    mannfraction = [1.0, 1.0, 1.0, 1.0]
    cxsdata = list(zip(xfraction, height, mannfraction))
    cxs = flopy.mf6.ModflowChfcxs(
        chf,
        nsections=1,
        npoints=4,
        packagedata=[(0, 4)],
        crosssectiondata=cxsdata,
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
    if inflow_location[idx] == "middle":
        inflow_cell = 4
    elif inflow_location[idx] == "left":
        inflow_cell = 0  # first cell
    else:
        raise Exception(f"invalid inflow location {inflow_location[idx]}")
    qinflow = dx * 1.0 / 86400.0
    flw = flopy.mf6.ModflowChfflw(
        chf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(inflow_cell, qinflow)],
    )

    return sim, None


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    chfname = "chf_model"

    # ensure outflow on left and right is the same
    fpth = test.workspace / f"{chfname}.zdg.obs.csv"
    obsvals = np.genfromtxt(fpth, names=True, delimiter=",")
    diff = obsvals["OUTFLOW1"] - obsvals["OUTFLOW9"]
    atol = 1.0e-6
    # This isn't working right now because the obs file is missing an E
    # when the numbers are very small.
    # if inflow_location[sim.idxsim] == "middle":
    #    assert np.allclose(diff, 0., atol=atol), f"{diff}"

    # read binary stage file
    fpth = test.workspace / f"{chfname}.stage"
    sobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage_all = sobj.get_alldata()
    for kstp, stage in enumerate(stage_all):
        print(kstp, stage.flatten())

    # burned in answer
    if inflow_location[idx] == "middle":
        stage_answer = [
            0.00610977,
            1.00610977,
            2.00610977,
            2.00610995,
            2.00610998,
            2.00610995,
            2.00610977,
            1.00610977,
            0.00610977,
        ]
    elif inflow_location[idx] == "left":
        stage_answer = [
            2.03013084,
            2.03013070,
            2.03013004,
            2.00928540,
            2.00928470,
            2.00928457,
            2.00928387,
            1.00928387,
            0.00928387,
        ]

    # print out the answer in a form that can be
    # dropped into this script, if necessary, as the answer
    for v in stage_all[-1].flatten():
        print(f"{v:18.8f},")

    msg = "Simulated stage does not match with the answer stored from a previous run."
    assert np.allclose(stage_all[-1].flatten(), stage_answer, atol=1.0e-5), msg


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
