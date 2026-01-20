import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sfr_div"]
inflows = np.array([10, 0, 10, 0, 10])
diversion = np.array([0.5, 0.5, 0.5, 0.5, 0.0])


def build_models(idx, test):
    # static model data
    # temporal discretization
    nper = len(inflows)
    tdis_rc = [(1.0, 1, 1.0)] * nper

    # spatial discretization data
    nlay, nrow, ncol = 1, 1, 1
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = -10.0
    strt = 0.0

    # build MODFLOW 6 files
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    sim.simulation_data.verify_data = False

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="days",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(sim, print_option="ALL")

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=name + ".cbb", saverecord=[["BUDGET", "ALL"]]
    )

    # sfr file
    cellid = (0, 0, 0)
    nreaches = 3
    rlen = 10.0
    rwid = 10.0
    roughness = 0.001
    rbth = 1.0
    rhk = 0.0
    slope = 0.001

    sfrrch_prop = [cellid, rlen, rwid, slope, top, rbth, rhk, roughness]
    packagedata = [
        [0] + sfrrch_prop + [2, 1.0, 1],
        [1] + sfrrch_prop + [1, 0.0, 0],
        [2] + sfrrch_prop + [1, 1.0, 0],
    ]
    connectiondata = [
        [0, -1, -2],
        [1, 0],
        [2, 0],
    ]
    diversiondata = [
        [0, 0, 1, "FRACTION"],
    ]
    perioddata = {
        i: [[0, "inflow", qin], [0, "diversion", 0, qdiv]]
        for i, (qin, qdiv) in enumerate(zip(inflows, diversion))
    }

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        budgetcsv_filerecord=name + ".sfr.csv",
        budget_filerecord=name + ".sfr.cbb",
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=diversiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )

    return sim, None


def check_output(idx, test):
    # check flow for individual reach
    fname = os.path.join(test.workspace, f"{test.name}.sfr.cbb")
    with flopy.utils.CellBudgetFile(fname) as cbb:
        outflows = cbb.get_data(text="EXT-OUTFLOW")

    # check outflow for reach 2 and 3
    assert np.allclose([r.q[1] for r in outflows], -inflows * diversion), (
        "Incorrect outflow for diversion reach"
    )
    assert np.allclose([r.q[2] for r in outflows], -inflows * (1 - diversion)), (
        "Incorrect outflow for outlet reach"
    )

    # load SFR budget CSV and check overall budget
    with open(fname.replace(".cbb", ".csv")) as f:
        header = f.readline().strip().split(",")
        flux = np.loadtxt(f, delimiter=",")

    assert np.allclose(flux[:, header.index("EXT-OUTFLOW_IN")], 0), (
        "External flow IN larger than zero"
    )
    assert np.allclose(flux[:, header.index("PERCENT_DIFFERENCE")], 0), (
        "Large mass balance error in SFR"
    )


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
