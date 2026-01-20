"""
Test the GWT Sorption (RCT) Package by running a simple 2-cell test with
mass injected into the first cell at a rate of 1.0 unit/time.  Compare the
simulated aqueous and sorbate concentrations with a known solution.  Problem
uses 10 time steps to discretize a 1.0 time unit period.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mst02a", "mst02b"]
distcoef = [0.0, 1.0]
nlay, nrow, ncol = 1, 1, 2

# time series answers for the two simulations
ts1 = np.array(
    [
        [0.1, 1.0, 0.0],
        [0.2, 2.0, 0.0],
        [0.3, 3.0, 0.0],
        [0.4, 4.0, 0.0],
        [0.5, 5.0, 0.0],
        [0.6, 6.0, 0.0],
        [0.7, 7.0, 0.0],
        [0.8, 8.0, 0.0],
        [0.9, 9.0, 0.0],
        [1.0, 10.0, 0.0],
    ]
)
ts2 = np.array(
    [
        [0.1, 0.09090909, 0.0],
        [0.2, 0.18181818, 0.0],
        [0.3, 0.27272727, 0.0],
        [0.4, 0.36363636, 0.0],
        [0.5, 0.45454545, 0.0],
        [0.6, 0.54545455, 0.0],
        [0.7, 0.63636364, 0.0],
        [0.8, 0.72727273, 0.0],
        [0.9, 0.81818182, 0.0],
        [1.0, 0.90909091, 0.0],
    ]
)
tsanswers = [ts1, ts2]


def build_models(idx, test):
    nper = 1
    perlen = [1.0]
    nstp = [10]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0]
    hnoflo = 1e30
    hdry = -1e30
    strt = [1.0]
    hk = 1.0
    laytyp = 1
    ss = 0.0
    sy = 0.1

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for id in range(nper):
        tdis_rc.append((perlen[id], nstp[id], tsmult[id]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk)

    # chd files
    chddict = {0: [[(0, 0, 0), 1.0]]}
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chddict, save_flows=False, pname="CHD-1"
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)

    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=sy,
        sorption="linear",
        bulk_density=1.0,
        sorbate_filerecord=f"{gwtname}.mst.csrb",
        distcoef=distcoef[idx],
    )

    # mass loading source
    srcdict = {0: [[(0, 0, 0), 1.0]]}
    src = flopy.mf6.ModflowGwtsrc(
        gwt, stress_period_data=srcdict, save_flows=False, pname="SRC-1"
    )

    # sources
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[[]], filename=f"{gwtname}.ssm")

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.bud",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name

    # Check aqueous concentrations
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        ts = cobj.get_ts([(0, 0, 0), (0, 0, 1)])
    except:
        assert False, f'could not load data from "{fpth}"'
    assert np.allclose(ts, tsanswers[idx]), (
        "simulated concentrations do not match with known solution."
    )

    # Check sorbate concentrations
    fpth = os.path.join(test.workspace, f"{gwtname}.mst.csrb")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="SORBATE")
        ts = cobj.get_ts([(0, 0, 0), (0, 0, 1)])
    except:
        assert False, f'could not load data from "{fpth}"'
    d = distcoef[idx]
    tsa_csrb = tsanswers[idx]
    tsa_csrb[:, 1:] *= d
    assert np.allclose(ts, tsa_csrb), (
        "Sorbate concentrations do not match with known solution.",
        ts,
    )

    # Check budget file
    fpth = os.path.join(test.workspace, f"{gwtname}.bud")
    try:
        bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        ra = bobj.get_data(totim=1.0)
    except:
        assert False, f'could not load data from "{fpth}"'


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
