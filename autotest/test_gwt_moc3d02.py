import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["moc3d02a", "moc3d02b"]
xt3d = [None, True]


def build_models(idx, test):
    nlay, nrow, ncol = 40, 12, 30
    nper = 1
    perlen = [400]
    nstp = [400]
    tsmult = [1.0]
    steady = [True]
    delr = 3.0
    delc = 0.5
    top = 0.0
    delz = 0.05
    botm = np.arange(-delz, -nlay * delz - delz, -delz)
    strt = 0.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 0.0125 / delz
    laytyp = 0
    diffc = 0.0
    alphal = 0.6
    alphath = 0.03
    alphatv = 0.006
    porosity = 0.25
    # ss = 0.
    # sy = 0.1

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-4, 1e-3, 0.97

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
    for k in range(nlay):
        for i in range(nrow):
            chdlist.append([(k, i, j), 0.0])
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chdlist, save_flows=False, pname="CHD-1"
    )

    # wel files
    wellist = []
    j = 0
    qwell = 0.1 * delz * delc * porosity
    for k in range(nlay):
        for i in range(nrow):
            wellist.append([(k, i, j), qwell, 0.0])
    wellist.append([(0, 0, 7), 1.0e-6, 2.5e6])  # source well

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
    strt[0, 0, 0] = 0.0
    ic = flopy.mf6.ModflowGwtic(gwt, strt=strt, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD", filename=f"{gwtname}.adv")

    # dispersion
    xt3d_off = not xt3d[idx]
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=xt3d_off,
        diffc=diffc,
        alh=alphal,
        alv=alphal,
        ath1=alphath,
        ath2=alphatv,
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
    name = cases[idx]
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
            6.19536860e01,
            6.10977621e01,
            5.94331187e01,
            5.70585587e01,
            5.41016012e01,
            5.06360797e01,
            4.67594881e01,
            4.26131186e01,
            3.83506623e01,
            3.41052863e01,
            2.99834184e01,
            2.60685066e01,
            2.24239296e01,
            1.90924262e01,
            1.60976346e01,
            1.34466507e01,
            1.11331123e01,
            9.14042952e00,
            7.44487655e00,
            6.01834933e00,
            4.83067716e00,
            3.85144623e00,
            3.05134631e00,
            2.40308745e00,
            1.88195336e00,
            1.46606535e00,
            1.13642906e00,
            8.76828593e-01,
            6.73622452e-01,
            5.15484358e-01,
            3.93121351e-01,
            2.98992183e-01,
            2.27041270e-01,
            1.72457391e-01,
            1.31461916e-01,
            1.01128195e-01,
            7.92317943e-02,
            6.41300864e-02,
            5.46692350e-02,
            5.01165648e-02,
        ]
    )

    csim = csim[:, 2, 12]
    # rtol is set larger here because convergence is looser
    assert np.allclose(cres, csim, rtol=1.0e-4), (
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
