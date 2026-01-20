"""
Test to make sure that array based recharge is applied correctly when idomain
is used with -1, 0, and 1 for top layer.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["rch03"]


def build_models(idx, test):
    nlay, nrow, ncol = 2, 4, 5
    perlen = [1.0]
    nper = len(perlen)
    nstp = nper * [1]
    tsmult = nper * [1.0]

    delr = delc = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "rch"

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # set ims csv files
    csv0 = f"{name}.outer.ims.csv"
    csv1 = f"{name}.inner.ims.csv"

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=csv0,
        csv_inner_output_filerecord=csv1,
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    idomain = [
        [
            [0, 0, 0, 0, 0],
            [0, -1, 1, -1, 0],
            [0, -1, 1, -1, 0],
            [0, 0, 0, 0, 0],
        ],
        [
            [1, 1, 1, 1, 1],
            [1, 1, 1, -1, 1],
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1],
        ],
    ]
    idomain = np.array(idomain, dtype=int)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=100.0,
        botm=[50.0, 0.0],
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=100.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=1.0)

    # chd
    chdspd = [[(1, 0, 0), 100.0]]
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

    recharge = np.arange(nrow * ncol).reshape(nrow, ncol) + 1.0
    irch = [
        [1, 0, 0, 0, 0],
        [0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0],
    ]
    rch = flopy.mf6.ModflowGwfrcha(gwf, print_flows=True, recharge=recharge, irch=irch)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, "rch.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="rch")[0]

    for i, (node, node2, q) in enumerate(records):
        print(f"Checking: index={i + 1}; node={node}; node2={node2}; q={q}")

    answer = [21, 27, 8, 32, 13, 34]
    errmsg = f"node must be {answer}.  found {records['node']}"
    assert np.allclose(records["node"], answer), errmsg

    answer = [1, 7, 8, 12, 13, 14]
    errmsg = f"node2 must be {answer}.  found {records['node2']}"
    assert np.allclose(records["node2"], answer), errmsg

    answer = [0.0, 7.0, 8.0, 12.0, 13.0, 14.0]
    errmsg = f"rech q must be {answer}.  found {records['q']}"
    assert np.allclose(records["q"], answer), errmsg

    fpth = os.path.join(test.workspace, "rch.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    heads = hobj.get_alldata()


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
