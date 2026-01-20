import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["moc3d03"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 30, 30
    nper = 1
    perlen = [1000]
    nstp = [1000]
    tsmult = [1.0]
    steady = [True]
    delr = 10.0
    delc = 10.0
    top = 0.0
    delz = 10.0
    botm = [-delz]
    strt = 0.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 3.6 / delz
    laytyp = 0
    diffc = 0.0
    alphal = 10.0
    alphath = 10.0
    alphatv = 10.0
    porosity = 0.20
    # ss = 0.
    # sy = 0.1

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-8, 1e-6, 1.0

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
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
    )

    # chd files
    chdlist = []
    j = ncol - 1
    for i in range(nrow - 1):
        chdlist.append([(0, i, j), 0.0])
    i = nrow - 1
    for j in range(ncol):
        chdlist.append([(0, i, j), 0.0])
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chdlist, save_flows=False, pname="CHD-1"
    )

    # wel files
    wellq = 56.25
    wellist = []
    wellist.append([(0, 0, 0), wellq, 1.0])  # source well

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=wellist,
        save_flows=False,
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

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
    strt = np.zeros((nlay, nrow, ncol))
    ic = flopy.mf6.ModflowGwtic(gwt, strt=strt, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD", filename=f"{gwtname}.adv")

    # dispersion
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        diffc=diffc,
        alh=alphal,
        alv=alphal,
        ath1=alphath,
        atv=alphatv,
        filename=f"{gwtname}.dsp",
    )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
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
    name = test.name
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        times = cobj.get_times()
        t = times[-1]
        csim = cobj.get_data(totim=t)
    except:
        assert False, f'could not load data from "{fpth}"'

    cres = np.array(
        [
            0.9999896,
            0.9999803,
            0.9999563,
            0.99990016,
            0.99977394,
            0.99949988,
            0.99892559,
            0.99776642,
            0.99551783,
            0.9913363,
            0.98390115,
            0.97129499,
            0.95097311,
            0.91991835,
            0.87507215,
            0.81406604,
            0.73614146,
            0.64297874,
            0.53905123,
            0.43118821,
            0.32731015,
            0.23468851,
            0.1583593,
            0.10028465,
            0.05950455,
            0.03306359,
            0.01721311,
            0.00841033,
            0.0038776,
            0.00177687,
        ]
    )

    assert np.allclose(cres, csim.diagonal().ravel()), (
        "simulated concentrations do not match with known solution."
    )


@pytest.mark.slow
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
