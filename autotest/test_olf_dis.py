"""

Test the dis2d and disv2d discretization packages so that they support
IDOMAIN and have valid binary grid files.

"""

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework

cases = [
    "olf-dis2d",
    "olf-disv2d",
]

# grid size
dx = 1000.0
nrow = 10
ncol = 10
ncpl = nrow * ncol
xorigin = 75.0
yorigin = 500.0
angrot = 32.0
botm = np.zeros((nrow, ncol), dtype=float)
idomain = np.ones((nrow, ncol), dtype=int)
idomain[4:6, 4:6] = 0


def build_models(idx, test):
    perlen = [1]  # 1 second
    nstp = [1]
    tsmult = [1.0]
    nper = len(perlen)

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

    if idx == 0:
        add_olf_model_dis2d(sim)
    elif idx == 1:
        add_olf_model_disv2d(sim)

    return sim, None


def add_olf_model_dis2d(sim):
    name = "overland"
    olf = flopy.mf6.ModflowOlf(sim, modelname=name, save_flows=True)

    dis = flopy.mf6.ModflowOlfdis2D(
        olf,
        nrow=nrow,
        ncol=ncol,
        delr=dx,
        delc=dx,
        bottom=botm,
        idomain=idomain,
        xorigin=xorigin,
        yorigin=yorigin,
        angrot=angrot,
    )

    dfw = flopy.mf6.ModflowOlfdfw(
        olf,
        print_flows=True,
        save_flows=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        manningsn=3.5,
        idcxs=None,
    )

    ic = flopy.mf6.ModflowOlfic(olf, strt=1.0)

    # output control
    oc = flopy.mf6.ModflowOlfoc(
        olf,
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

    chd = flopy.mf6.ModflowOlfchd(
        olf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(0, 0, 1.0), (nrow - 1, ncol - 1, 0.5)],
    )

    return


def add_olf_model_disv2d(sim):
    name = "overland"
    olf = flopy.mf6.ModflowOlf(sim, modelname=name, save_flows=True)

    disvkwargs = get_disv_kwargs(
        1,
        nrow,
        ncol,
        dx,
        dx,
        1.0,  # top
        0.0,  # botm
        0.0,  # xoffset
        0.0,  # yoffset
    )

    _ = disvkwargs.pop("top")
    _ = disvkwargs.pop("nlay")
    bottom = disvkwargs.pop("botm").reshape((ncpl,))
    disvkwargs["nodes"] = disvkwargs.pop("ncpl")

    dis = flopy.mf6.ModflowOlfdisv2D(
        olf,
        bottom=bottom,
        idomain=idomain.reshape((ncpl,)),
        xorigin=xorigin,
        yorigin=yorigin,
        angrot=angrot,
        **disvkwargs,
    )

    dfw = flopy.mf6.ModflowOlfdfw(
        olf,
        print_flows=True,
        save_flows=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        manningsn=3.5,
        idcxs=None,
    )

    ic = flopy.mf6.ModflowOlfic(olf, strt=1.0)

    # output control
    oc = flopy.mf6.ModflowOlfoc(
        olf,
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

    chd = flopy.mf6.ModflowOlfchd(
        olf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(0, 1.0), (ncpl - 1, 0.5)],
    )

    return


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    mfsim = test.sims[0]
    olf = mfsim.olf[0]

    stage = olf.output.stage().get_data()

    pmv = flopy.plot.PlotMapView(model=olf)
    pmv.plot_array(stage, masked_values=[3e30])
    pmv.plot_grid()

    fname = test.workspace / "results.png"
    plt.savefig(fname)


def check_grb_dis2d(fpth):
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    assert grb.grid_type == "DIS2D", "grb grid type not DIS2D"
    assert grb.ncells == nrow * ncol, "grb ncells is incorrect"
    assert grb.nrow == nrow, "nrow in grb file is not 10"
    assert grb.ncol == ncol, "ncol in grb file is not 10"
    assert grb.nja == 432, "nja in grb file is not 432"
    assert grb.xorigin == xorigin, "xorigin in grb file is not correct"
    assert grb.yorigin == yorigin, "yorigin in grb file is not correct"
    assert grb.angrot == angrot, "angrot in grb file is not correct"
    assert np.allclose(grb.delr, dx * np.ones(ncol)), "grb delr not correct"
    assert np.allclose(grb.delc, dx * np.ones(nrow)), "grb delc not correct"
    assert np.allclose(grb.bot.reshape((nrow, ncol)), np.zeros((nrow, ncol))), (
        "grb botm not correct"
    )
    assert grb.ia.shape[0] == grb.ncells + 1, "ia in grb file is not correct size"
    assert grb.ja.shape[0] == grb.nja, "ja in grb file is not correct size"
    assert np.allclose(grb.idomain.reshape((nrow, ncol)), idomain), (
        "grb idomain not correct"
    )


def check_grb_disv2d(fpth):
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    assert grb.grid_type == "DISV2D", "grb grid type not DISV2D"
    assert grb.ncells == ncpl, "grb ncells is incorrect"
    assert grb._datadict["NODES"] == 96, "grb nodes is incorrect"
    assert grb.verts.shape == (
        (nrow + 1) * (ncol + 1),
        2,
    ), "vertices shape is incorrect"
    assert grb.nja == 432, "nja in grb file is not 432"
    assert grb.xorigin == xorigin, "xorigin in grb file is not correct"
    assert grb.yorigin == yorigin, "yorigin in grb file is not correct"
    assert grb.angrot == angrot, "angrot in grb file is not correct"
    assert np.allclose(grb.bot.reshape((nrow, ncol)), np.zeros((nrow, ncol))), (
        "grb botm not correct"
    )
    cellx, celly = np.meshgrid(
        np.linspace(dx / 2, ncol * dx - dx / 2, ncol),
        np.linspace(dx * nrow - dx / 2.0, dx / 2, nrow),
    )
    assert np.allclose(grb._datadict["CELLX"], cellx.flatten()), "cellx is not right"
    assert np.allclose(grb._datadict["CELLY"], celly.flatten()), "celly is not right"
    assert grb._datadict["IAVERT"].shape[0] == ncpl + 1, "iavert size not right"
    assert grb._datadict["IAVERT"][-1] - 1 == grb._datadict["JAVERT"].shape[0], (
        "javert size not right"
    )
    assert grb.ia.shape[0] == grb.ncells + 1, "ia in grb file is not correct size"
    assert grb.ja.shape[0] == grb.nja, "ja in grb file is not correct size"
    assert np.allclose(grb.idomain.reshape((ncpl,)), idomain.reshape((ncpl,))), (
        "grb idomain not correct"
    )


def check_output(idx, test):
    print(f"evaluating model for case {idx}...")

    modelname = "overland"
    ws = test.workspace
    mfsim = flopy.mf6.MFSimulation.load(sim_ws=ws)

    # read the binary grid file
    if idx == 0:
        fpth = test.workspace / f"{modelname}.dis2d.grb"
        check_grb_dis2d(fpth)
    elif idx == 1:
        fpth = test.workspace / f"{modelname}.disv2d.grb"
        check_grb_disv2d(fpth)

    # read binary stage file
    fpth = test.workspace / f"{modelname}.stage"
    sobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage = sobj.get_data().reshape((nrow, ncol))
    assert np.allclose(stage[idomain == 0], 3.0e30), (
        "stage should have nodata values where idomain is zero"
    )
    assert stage[idomain == 1].max() == 1.0, "maximum stage should be 1.0"
    assert stage[idomain == 1].min() == 0.5, "minimum stage should be 0.5"


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
