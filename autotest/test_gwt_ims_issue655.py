"""
Test the GWT with initial concentration of 0 in the model domain and right-hand
side GWT boundary conditions. Versions 6.2.1 and earlier failed because of a
divide by zero error in IMS. This test confirms the fix implemented as part
of the version 6.2.2 release that addressed Issue 655.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["issue655a", "issue655b"]
newton = [False, True]
laytyp = [1]
ss = [1.0e-10]
sy = [0.1]
nlay, nrow, ncol = 1, 11, 11


def build_models(idx, test):
    nper = 1
    perlen = [1000.0]
    nstp = [5]
    tsmult = [1.05]
    xlen = 8020.0
    ylen = 8020.0
    delr = xlen / float(ncol)
    delc = ylen / float(nrow)
    botm = [0.0]
    strt = 20.0
    hk = 2834.8

    ss = 1.0e-6
    sy = 0.1

    top = 20.0
    laytyp = 0

    # put constant heads all around the box
    chdlist = (
        [(0, i, 0, strt) for i in range(nrow)]
        + [(0, i, ncol - 1, strt) for i in range(nrow)]
        + [(0, 0, i, strt) for i in range(ncol)]
        + [(0, nrow - 1, i, strt) for i in range(ncol)]
    )
    chdspdict = {0: list(set(chdlist))}

    # injection well with rate (1000 GPM) and concentration of 100.
    q = 192528.0
    c = 100
    w = {0: [[(0, int(nrow / 2), int(ncol / 2)), q, c]]}

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    if newton[idx]:
        newtonoptions = "NEWTON"
        linear_acceleration = "bicgstab"
    else:
        newtonoptions = None
        linear_acceleration = "cg"

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions=newtonoptions,
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
        linear_acceleration=linear_acceleration,
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)

    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        # steady_state={0: True, 2: True},
        transient={0: True},
    )

    # chd files
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chdspdict, save_flows=False, pname="CHD-1"
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=w,
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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwtname}.adv")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

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
        saverecord=[("CONCENTRATION", "LAST")],
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
    gwfname = "gwf_" + name

    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        head = hobj.get_alldata().flatten()
    except:
        assert False, f'could not load data from "{fpth}"'

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_alldata().flatten()
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculations
    times = hobj.get_times()

    # vold = 1000 * 0.1
    # massold = 100.0 * vold  # c * v
    # nstp = len(times)
    # dt = 4.0 / nstp
    # for t, h, c in zip(times, head, conc):
    #     q = -25.0
    #     if t > 2.0:
    #         q = 25.0
    #     v = vold + q * dt
    #     conc_calc = massold / v
    #     # print(t, h, c, v, conc_calc)
    #     msg = "simulated conc ({}) /= calc conc ({}) for time {}"
    #     msg = msg.format(c, conc_calc, t)
    #     assert np.allclose(c, conc_calc, atol=0.001), msg
    #     vold = v


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
