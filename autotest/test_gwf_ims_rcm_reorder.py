import os

import flopy
import pytest
from flopy.utils.compare import eval_bud_diff
from framework import TestFramework

paktest = "ims"
cases = ["ims_rcm"]
cmp_prefix = "mf6"

# spatial discretization data
nlay, nrow, ncol = 2, 5, 30
delr, delc = 100.0, 100.0
top = 0.0
botm = [-10.0, -20.0]
strt = 0.0
chd_left = 10.0
chd_right = 5.0


def build_model(idx, ws):
    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = [(1.0, 1, 1.0)]

    # build MODFLOW 6 files
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=nper,
        perioddata=tdis_rc,
    )

    if not str(ws).endswith(cmp_prefix):
        reordering_method = "rcm"
    else:
        reordering_method = None

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        reordering_method=reordering_method,
        preconditioner_levels=10,
        preconditioner_drop_tolerance=1e-6,
        outer_dvclose=1e-9,
        outer_maximum=100,
        inner_dvclose=1e-12,
        inner_maximum=100,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf)

    # chd files
    # chd data
    spd = [[(0, i, 0), chd_left] for i in range(nrow)]
    spd += [[(0, i, ncol - 1), chd_right] for i in range(nrow)]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd, pname="chd-1")

    # output control
    hdspth = f"{name}.hds"
    budpth = f"{name}.cbc"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=hdspth,
        budget_filerecord=budpth,
        printrecord=[
            ("BUDGET", "ALL"),
        ],
        saverecord=[
            ("BUDGET", "ALL"),
            ("HEAD", "ALL"),
        ],
    )

    return sim


def build_models(idx, test):
    return build_model(idx, test.workspace), build_model(
        idx, os.path.join(test.workspace, cmp_prefix)
    )


def check_output(idx, test):
    name = test.name
    fpth = os.path.join(test.workspace, f"{name}.dis.grb")
    ia = flopy.mf6.utils.MfGrdFile(fpth).ia

    fpth = os.path.join(test.workspace, f"{name}.cbc")
    b0 = flopy.utils.CellBudgetFile(fpth, precision="double")

    fpth = os.path.join(test.workspace, cmp_prefix, f"{name}.cbc")
    b1 = flopy.utils.CellBudgetFile(fpth, precision="double")

    fpth = os.path.join(test.workspace, f"{name}.cbc.cmp.out")
    eval_bud_diff(fpth, b0, b1, ia=ia)

    # close the budget files
    b0.close()
    b1.close()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
