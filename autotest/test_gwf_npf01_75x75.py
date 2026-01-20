import os
import pathlib as pl

import flopy
import numpy as np
import pytest
from conftest import try_get_target
from framework import TestFramework

cases = ["npf01a_75x75", "npf01b_75x75"]
top = [100.0, 0.0]
laytyp = [1, 0]
ss = [0.0, 1.0e-4]
sy = [0.1, 0.0]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 75, 75
    nper = 3
    perlen = [1.0, 1000.0, 1.0]
    nstp = [1, 10, 1]
    tsmult = [1.0, 1.5, 1.0]
    steady = [True, False, True]
    lenx = 20000.0
    delr = delc = lenx / float(nrow)
    botm = [-100.0]
    strt = 40.0
    hnoflo = 1e30
    hdry = -1e30
    mu = 5.0
    sigma = 1.23
    np.random.seed(seed=9001)
    hk = np.random.lognormal(mu, sigma, (nrow, ncol))

    nc = int((nrow - 1) / 2) + 1
    w = [(0, nc, nc, -1000.0)]
    wd = {1: w}
    ws = [((0, nc, nc), -1000.0)]
    wd6 = {1: ws}

    c = []
    c6 = []
    for i in range(nrow):
        c.append([0, i, 0, 48.0, 48.0])
        c.append([0, i, ncol - 1, 40.0, 40.0])
        c6.append([(0, i, 0), 48.0])
        c6.append([(0, i, ncol - 1), 40.0])
    cd = {0: c}
    cd6 = {0: c6}

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 0.01, 1.0

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
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )

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
        export_array_ascii=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top[idx],
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(
        gwf, export_array_ascii=True, strt=strt, filename=f"{name}.ic"
    )

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        export_array_ascii=True,
        save_flows=False,
        icelltype=laytyp[idx],
        k=hk,
        k33=hk,
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp[idx],
        ss=ss[idx],
        sy=sy[idx],
        steady_state={0: True, 1: False, 2: True},
        transient={1: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=len(c6), stress_period_data=cd6, save_flows=False
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(str(ws)),
        stress_period_data=wd6,
        save_flows=False,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # build MODFLOW-2005 files
    ws = os.path.join(test.workspace, "mf2005")
    mc = flopy.modflow.Modflow(
        name, model_ws=ws, exe_name=try_get_target(test.targets, "mf2005")
    )
    dis = flopy.modflow.ModflowDis(
        mc,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        tsmult=tsmult,
        steady=steady,
        delr=delr,
        delc=delc,
        top=top[idx],
        botm=botm,
    )
    bas = flopy.modflow.ModflowBas(mc, ibound=1, strt=strt, hnoflo=hnoflo)
    lpf = flopy.modflow.ModflowLpf(
        mc,
        laytyp=laytyp[idx],
        hk=hk,
        vka=hk,
        ss=ss[idx],
        sy=sy[idx],
        constantcv=True,
        hdry=hdry,
    )
    chd = flopy.modflow.ModflowChd(mc, stress_period_data=cd)
    wel = flopy.modflow.ModflowWel(mc, stress_period_data=wd)
    oc = flopy.modflow.ModflowOc(mc, stress_period_data=None)
    pcg = flopy.modflow.ModflowPcg(
        mc,
        mxiter=nouter,
        iter1=ninner,
        hclose=hclose,
        rclose=rclose,
        relax=relax,
    )
    return sim, mc


def check_output(idx, test):
    print("evaluating model...")
    ws = test.workspace

    # ensure export array is working properly
    name = cases[idx]
    layered = [
        "dis.botm",
        "dis.idomain",
        "ic.strt",
        "npf.icelltype",
        "npf.k",
        "npf.k33",
    ]
    flist = [
        "dis.botm",
        "dis.delc",
        "dis.delr",
        "dis.idomain",
        "dis.top",
        "ic.strt",
        "npf.icelltype",
        "npf.k",
        "npf.k33",
    ]
    files = [
        (
            pl.Path(ws / f"{name}-{f}.l1.txt")
            if f in layered
            else pl.Path(ws / f"{name}-{f}.txt")
        )
        for f in flist
    ]
    gwf = test.sims[0].gwf[0]
    for i, fpth in enumerate(files):
        assert fpth.is_file(), f"Expected file does not exist: {fpth.name}"
        a = np.loadtxt(fpth)
        array_name = flist[i][flist[i].index(".") + 1 :]
        package_name = flist[i][0 : flist[i].index(".")]
        package = getattr(gwf, package_name)
        b = getattr(package, array_name).array
        assert np.allclose(a, b)
        print(f"compared: {fpth}")
        print(f"a={a}")
        print(f"b={b}")
    return


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
