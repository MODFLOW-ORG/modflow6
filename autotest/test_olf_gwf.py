"""

Simple 3-cell OLF model connected to a simple
3-cell GWF model


Zero-based diagram below

   0      1      2      OLF
o------o------o------o
|      |      |      |
|  0   |  1   |  2   |  GWF
|      |      |      |
o------o------o------o


"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "olf-gwf01",
]


def build_models(idx, test):
    sim_ws = test.workspace
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=sim_ws,
        memory_print_option="all",
    )

    tdis = flopy.mf6.ModflowTdis(sim)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="all",
        linear_acceleration="BICGSTAB",
        outer_dvclose=1.0e-12,
        inner_dvclose=1.0e-12,
        outer_maximum=500,
        under_relaxation="simple",
        under_relaxation_gamma=0.1,
    )

    add_olf_model(sim)
    add_gwf_model(sim)

    olfgwf_data = [
        ((0, 0), (0, 0, 0), 1.0, 1.0),
        ((0, 1), (0, 0, 1), 1.0, 1.0),
        ((0, 2), (0, 0, 2), 1.0, 1.0),
    ]
    olfgwf = flopy.mf6.ModflowOlfgwf(
        sim,
        print_input=True,
        print_flows=True,
        exgtype="OLF6-GWF6",
        nexg=len(olfgwf_data),
        exgmnamea="olfmodel",
        exgmnameb="gwfmodel",
        exchangedata=olfgwf_data,
    )

    return sim, None


def add_olf_model(sim):
    name = "olfmodel"
    olf = flopy.mf6.ModflowOlf(sim, modelname=name, save_flows=True)

    dx = 1000.0
    nrow = 1
    ncol = 3
    total_length = dx * ncol

    botm = np.zeros((nrow, ncol), dtype=float)

    dis = flopy.mf6.ModflowOlfdis2D(
        olf,
        nrow=nrow,
        ncol=ncol,
        delr=dx,
        delc=dx,
        bottom=botm,
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

    sto = flopy.mf6.ModflowOlfsto(olf, save_flows=True)

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

    flw_spd = [
        (0, 0, 100),
    ]
    maxbound = len(flw_spd)
    flw = flopy.mf6.ModflowOlfflw(
        olf,
        maxbound=maxbound,
        print_input=True,
        print_flows=True,
        stress_period_data=flw_spd,
    )

    chd = flopy.mf6.ModflowOlfchd(
        olf,
        maxbound=1,
        print_input=True,
        print_flows=True,
        stress_period_data=[(0, 2, 1.0)],
    )

    return


def add_gwf_model(sim):
    # create gwf model
    name = "gwfmodel"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=3,
        delr=1000.0,
        delc=50.0,
        top=0.0,
        botm=-10.0,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=-5.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=1,
        k=1.0,
    )

    sto = flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1.0e-5, iconvert=1)

    # chd files
    # chd_spd = [(0, 0, 0, 1.0), (0, 0, 2, 1.0)]
    # chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
    #     gwf,
    #     maxbound=len(chd_spd),
    #     stress_period_data=chd_spd,
    #     pname="CHD-1",
    # )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )


def check_output(idx, test):
    print("evaluating model...")

    # assign name
    name = "olfmodel"

    # read the binary grid file
    fpth = test.workspace / f"{name}.dis2d.grb"
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    ia = grb.ia
    ja = grb.ja
    assert ia.shape[0] == grb.ncells + 1, "ia in grb file is not correct size"

    # read stage file
    fpth = test.workspace / f"{name}.stage"
    qobj = flopy.utils.HeadFile(fpth, precision="double", text="STAGE")
    stage = qobj.get_alldata()
    print(f"{stage=}")

    # read the budget file
    fpth = test.workspace / f"{name}.bud"
    budobj = flopy.utils.binaryfile.CellBudgetFile(fpth)
    flowja = budobj.get_data(text="FLOW-JA-FACE")[0].flatten()
    qstorage = budobj.get_data(text="STORAGE")
    qflw = budobj.get_data(text="FLW")
    qchd = budobj.get_data(text="CHD")
    qresidual = np.zeros(grb.nodes)
    nodes = ia.shape[0] - 1
    for n in range(nodes):
        qresidual[n] = flowja[ia[n]]

    print(budobj.list_records())
    print(f"flowja: {flowja}")
    print(f"qstorage: {qstorage}")
    print(f"qchd: {qchd}")
    print(f"qflw: {qflw}")
    print(f"qresidual: {qresidual}")

    assert np.allclose(qresidual, 0.0, atol=1.0e-3), "Flowja residual > 1.e-3"

    return


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
