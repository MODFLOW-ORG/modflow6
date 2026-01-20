import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["aux02"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 10, 10
    nper = 3
    perlen = [1.0, 1.0, 1.0]
    nstp = [10, 10, 10]
    tsmult = [1.0, 1.0, 1.0]

    lenx = 300.0
    delr = delc = lenx / float(nrow)
    strt = 100.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
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
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=90.0,
        botm=0.0,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=1.0, k33=0.01)

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 100.0] + [a for a in range(100)])
    chdlist0.append([(0, nrow - 1, ncol - 1), 95.0] + [a for a in range(100)])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=True,
        auxiliary=[f"aux{i}" for i in range(100)],
        pname="CHD-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.bud",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    # maw budget aux variables
    fpth = os.path.join(test.workspace, "aux02.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="CHD")
    for r in records:
        for a in range(100):
            aname = f"AUX{a}"
            assert np.allclose(r[aname], a)


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
