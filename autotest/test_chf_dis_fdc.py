"""

Test the fdc input parameter for disv1d to shift the node location along the
cell.  Use a 2 cell model and assign different roughness values for cell 1
and cell 2.  Fix the head in both cells and make a hand calculation of the flow.
Compare the hand calculation of flow with the simulated flow.

"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "chf-dis-fdc",
]

# grid size
dx = 1000.0
nreach = 2
total_length = dx * nreach
vertices = []
vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
cell1d = []
fdc = [0.0, 1.0]
for j in range(nreach):
    cell1d.append([j, fdc[j], 2, j, j + 1])
nodes = len(cell1d)
nvert = len(vertices)
xorigin = 100.0
yorigin = 200.0
angrot = 25.0

idomain = np.ones((nodes,), dtype=int)
h0 = 0.5
h1 = 1.0
rough = [0.035, 0.35]


def build_models(idx, test):
    perlen = [1]  # 1 second
    nstp = [1]
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
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    disv1d = flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=nodes,
        nvert=nvert,
        width=50.0,
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
        manningsn=rough,
        idcxs=None,
    )

    ic = flopy.mf6.ModflowChfic(chf, strt=1.0)

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

    chd = flopy.mf6.ModflowChfchd(
        chf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(0, h0), (nodes - 1, h1)],
    )

    return


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    mfsim = test.sims[0]
    chf = mfsim.chf[0]

    stage = chf.output.stage().get_data().reshape((nodes,))

    pmv = flopy.plot.PlotMapView(model=chf)
    pmv.plot_array(stage, masked_values=[3e30])  # not working yet
    pmv.plot_grid()

    fname = test.workspace / "results.png"
    plt.savefig(fname)


def check_grb_disv1d(fpth):
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    assert grb.grid_type == "DISV1D", "grb grid type not DISV1D"
    assert grb.ncells == nodes, "grb ncells is incorrect"
    assert grb._datadict["NCELLS"] == nodes, "grb nodes is incorrect"
    assert grb.verts.shape == (nodes + 1, 2), "vertices shape is incorrect"
    assert grb.nja == 4, "nja in grb file is not 4"
    assert grb.xorigin == xorigin, "xorigin in grb file is not correct"
    assert grb.yorigin == yorigin, "yorigin in grb file is not correct"
    assert grb.angrot == angrot, "angrot in grb file is not correct"
    assert np.allclose(grb.bot.reshape((nodes,)), np.zeros((nodes,))), (
        "grb botm not correct"
    )
    cellx = np.array([0.0, 2 * dx])  # node centers pushed all the way to left and right
    celly = np.zeros(nreach)
    assert np.allclose(grb._datadict["CELLX"], cellx.flatten()), "cellx is not right"
    assert np.allclose(grb._datadict["CELLY"], celly.flatten()), "celly is not right"
    assert grb._datadict["IAVERT"].shape[0] == nodes + 1, "iavert size not right"
    assert grb._datadict["IAVERT"][-1] - 1 == grb._datadict["JAVERT"].shape[0], (
        "javert size not right"
    )
    assert grb.ia.shape[0] == grb.ncells + 1, "ia in grb file is not correct size"
    assert grb.ja.shape[0] == grb.nja, "ja in grb file is not correct size"
    assert np.allclose(grb.idomain.reshape((nodes,)), idomain.reshape((nodes,))), (
        "grb idomain not correct"
    )


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    modelname = "channel"
    ws = test.workspace
    mfsim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    chf = mfsim.get_model(modelname)

    # check binary grid file
    fpth = test.workspace / f"{modelname}.disv1d.grb"
    check_grb_disv1d(fpth)

    # read binary stage file
    stage = chf.output.stage().get_data().reshape((nodes,))
    assert stage[idomain == 1].max() == 1.0, "maximum stage should be 1.0"
    assert stage[idomain == 1].min() == 0.5, "minimum stage should be 0.5"

    # extract the simulated flow from the budget file
    fpth = test.workspace / f"{modelname}.bud"
    sobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    flowjaface = sobj.get_data(text="FLOW-JA-FACE")[-1].flatten()
    flow_sim = flowjaface[1]

    # make an independent calculate of the flow between two constant head
    # cells
    def get_cond_n(depth, width, rough, dhds):
        unitconv = 86400.0
        a = depth * width
        rh = depth
        conveyance = a * rh ** (2.0 / 3.0) / rough
        dhds_sqr = dhds**0.5
        c = unitconv * conveyance / dx / dhds_sqr
        return c

    cln = 1000.0
    clm = 1000.0
    dhds = abs(h1 - h0) / (cln + clm)
    depth_n = h0 - 0.0
    depth_m = h1 - 0.0
    if h0 > h1:
        depth_m = depth_n
    else:
        depth_n = depth_m

    width = 50.0
    rough_n, rough_m = rough
    cn = get_cond_n(depth_n, width, rough_n, dhds)
    cm = get_cond_n(depth_m, width, rough_m, dhds)
    cond = cn * cm / (cn + cm)
    flow = cond * (h1 - h0)
    print(f"{cn=} {cm=} {cond=}")
    print(f"Known flow is {flow} cubic meters per seconds")
    print(f"Simulated flow is {flow_sim} cubic meters per seconds")
    assert np.allclose(flow, flow_sim), "known flow and simulated flow not the same"


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
