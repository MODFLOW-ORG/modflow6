import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["newton01"]
nlay = 2
nrow, ncol = 3, 3
top = 20
botm = [10, 0]
laytyp = 1
hk = 10.0
delr = delc = 1.0
chdloc = [(1, i, j) for i in range(nrow) for j in range(ncol)]
chd = 7.0
strt = chd
rch = 1.0
oname = "head_obs.csv"
obs_recarray = {oname: [("h1", "HEAD", (0, 1, 1)), ("h2", "HEAD", (1, 1, 1))]}


def build_models(idx, test):
    c6 = []
    for loc in chdloc:
        c6.append([loc, chd])
    cd6 = {0: c6}

    nper = 1
    tdis_rc = [(1.0, 1, 1.0)]

    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        inner_dvclose=1e-9,
        outer_dvclose=1e-9,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )

    flopy.mf6.ModflowGwfdis(
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk)

    # gwf observation
    flopy.mf6.ModflowUtlobs(gwf, digits=10, print_input=True, continuous=obs_recarray)

    # chd files
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=cd6)

    # rch files
    flopy.mf6.modflow.ModflowGwfrcha(gwf, recharge={0: rch})

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, oname)
    v = np.genfromtxt(fpth, delimiter=",", names=True)

    msg = f"head in layer 1 != 8. ({v['H1']})"
    assert np.allclose(v["H1"], 8.0), msg

    msg = f"head in layer 2 != 7. ({v['H2']})"
    assert np.allclose(v["H2"], 7.0), msg


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
