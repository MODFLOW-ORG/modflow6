"""
Test to make sure that auxmultcol is working when used with a time series
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["auxmult01", "auxmult02"]

wellist = [[(0, 2, 2), "tsq", "tsqfact"], [(0, 2, 2), 1.0000000, "tsqfact"]]


def build_models(idx, test):
    global numstep
    nlay, nrow, ncol = 1, 3, 3
    perlen = [1.0, 1.0, 1.0, 1.0]
    nstp = [10, 1, 1, 1]
    tsmult = [1.0, 1.0, 1.0, 1.0]
    nper = len(perlen)
    lenx = 300.0
    delr = delc = lenx / float(nrow)
    botm = -1.0
    hk = 1.0
    numstep = sum(nstp)

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-3, 1.0

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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
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
        top=0.0,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=0,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 1.0])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        pname="wel",
        print_input=True,
        print_flows=True,
        stress_period_data={0: [wellist[idx]]},
        auxiliary=["auxmult"],
        auxmultname="auxmult",
    )
    # ts_filerecord='well-rates.ts')

    # wel obs
    obs = {
        "wel.obs.csv": [["q", "wel", (0, 2, 2)]],
    }
    welobs = wel.obs.initialize(
        print_input=True,
        continuous=obs,
    )

    # well ts package
    ts_recarray = [
        (0.0, 0.0, 1.0),
        (0.1, 1.0, 1.0),
        (0.2, 0.0, 1.0),
        (0.3, 1.0, 1.0),
        (0.4, 0.0, 1.0),
        (0.5, 1.0, 1.0),
        (0.6, 0.0, 1.0),
        (0.7, 1.0, 1.0),
        (0.8, 0.0, 1.0),
        (0.9, 1.0, 1.0),
        (1.0, 0.0, 1.0),
        (2.0, 1.0, 1.0),
        (3.0, 0.0, 1.0),
        (4.0, 1.0, 1.0),
    ]

    filename = name + ".wel.ts"
    time_series_namerecord = [("tsqfact", "tsq")]
    interpolation_methodrecord = [("linearend", "linearend")]
    wel.ts.initialize(
        filename=filename,
        timeseries=ts_recarray,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
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
    name = cases[idx]
    fpth = os.path.join(test.workspace, f"{name}.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double", verbose=False)
    records = bobj.get_data(text="wel")

    qlist = []
    for kper, r in enumerate(records):
        node, node2, q, qauxmult = r[0]
        qlist.append(q)
    qlist = np.array(qlist)

    answer = np.array(7 * [1.0, 0.0])[:-1]
    msg = f"err {qlist} /= {answer}"
    assert np.allclose(qlist, answer), msg

    # MODFLOW 6 observations
    fpth = os.path.join(test.workspace, "wel.obs.csv")
    try:
        obs = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    rate = obs["Q"]
    obs_answer = [1.0 if x % 2 == 0 else 0.0 for x in range(numstep)]
    msg = f"err {rate} /= {obs_answer}"
    assert np.allclose(rate, obs_answer), msg


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
