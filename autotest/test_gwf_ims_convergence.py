import os

import flopy
import pandas as pd
import pytest
from framework import TestFramework

outers_record = "outers.csv"

cases = ["ims_strict", "ims_noinnercnvg"]
xfail = [False, True]
outer_dvclose = [1.0e3, 1.0e-7]  # to activate IMS strict
inner_dvclose = [1.0e-07, 1.0e-100]  # large to force non-convergence on inners
rclose_rec = ["1e-07 strict", "1e-07"]
outer_max = 10
inner_max = 20

# spatial discretization data
nlay, nrow, ncol = 2, 10, 10
delr, delc = 100.0, 100.0
top = 0.0
botm = [-10.0, -20.0]
strt = 0.0
chd_left = 10.0
chd_right = 5.0


def build_models(idx, test):
    ws = test.workspace

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
        time_units="days",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=outers_record,
        outer_dvclose=outer_dvclose[idx],
        outer_maximum=outer_max,
        inner_dvclose=inner_dvclose[idx],
        inner_maximum=inner_max,
        rcloserecord=rclose_rec[idx],
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=1.0)

    # chd files
    # chd data
    spd = [[(0, irow, 0), chd_left] for irow in range(nrow)]
    spd += [[(0, irow, ncol - 1), chd_right] for irow in range(nrow)]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd, pname="chd-1")

    return sim, None


def check_output(idx, test):
    csv_file = os.path.join(test.workspace, outers_record)
    assert os.path.exists(csv_file), "outer iterations are stored in csv"

    df = pd.read_csv(csv_file)

    if idx == 0:  # ims_strict
        assert len(df["solution_outer_dvmax"].values) > 1, (
            "can't converge on first outer with 'strict'"
        )
        assert df["inner_iterations"].values[-1] == 1, (
            "strict: single inner iteration on last outer"
        )

    if idx == 1:  # ims_noinnercnvg
        assert df["total_inner_iterations"].values[-1] == outer_max * inner_max, (
            "test should fail, but after trying for maximum nr. of iters"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=xfail[idx],
    )
    test.run()
