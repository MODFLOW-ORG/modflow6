import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["henry01"]


def build_models(idx, test):
    lx = 2.0
    lz = 1.0

    nlay = 10
    nrow = 1
    ncol = 20
    nper = 1
    delr = lx / ncol
    delc = 1.0
    top = 1.0
    delz = lz / nlay
    botm = list(top - np.arange(delz, nlay * delz + delz, delz))

    perlen = [0.5]
    nstp = [500]
    tsmult = [1.0]
    steady = [True]
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-10, 1e-6, 0.97

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
    gwtname = "gwt_" + name

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwfname])

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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=1.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=864.0,
    )

    pd = [(0, 0.7, 0.0, gwtname, "none")]
    buy = flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    def chd_value(k):
        # depth = k * delz + 0.5 * delz
        # hf = top + 0.025 * depth
        hf = top
        return hf

    # chd files
    chdlist1 = []
    for k in range(nlay):
        chdlist1.append([(k, 0, ncol - 1), chd_value(k), 35.0])
    chd1 = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.chd",
    )

    wellist1 = []
    qwell = 5.7024 / nlay
    for k in range(nlay):
        wellist1.append([(k, 0, 0), qwell, 0.0])
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wellist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.wel",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=35.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # dispersion
    diffc = 0.57024
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        diffc=diffc,
        # alh=0., alv=0., ath=0., atv=0.,
        filename=f"{gwtname}.dsp",
    )

    # mass storage and transfer
    porosity = 0.35
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")

    # sources
    sourcerecarray = [
        ("CHD-1", "AUX", "CONCENTRATION"),
        ("WEL-1", "AUX", "CONCENTRATION"),
    ]
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
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # time step 500 and only for the bottom layer.

    # these commented out concs are for equivalent freshwater head formulation
    # cres = [[2.89573493e-05, 1.14446653e-04, 4.51581911e-04, 1.76640638e-03,
    #         6.82342381e-03, 2.59046347e-02, 9.59028163e-02, 3.41705012e-01,
    #         1.14554311e+00, 3.47615268e+00, 8.98354607e+00, 1.51879954e+01,
    #         2.06578219e+01, 2.51309261e+01, 2.85181814e+01, 3.09868580e+01,
    #         3.26940512e+01, 3.37932558e+01, 3.44410870e+01, 3.47819736e+01]]

    cres = [
        [
            6.62568993e-05,
            2.62136709e-04,
            1.03376866e-03,
            4.03149437e-03,
            1.54745356e-02,
            5.80960685e-02,
            2.11109069e-01,
            7.29225328e-01,
            2.32090039e00,
            6.46492264e00,
            1.28547693e01,
            1.86965727e01,
            2.36023472e01,
            2.74207659e01,
            3.02383176e01,
            3.22172463e01,
            3.35101828e01,
            3.42835666e01,
            3.46963233e01,
            3.48890176e01,
        ]
    ]

    cres = np.array(cres)
    assert np.allclose(cres, conc[-1, :, :]), (
        "simulated concentrations do not match with known solution.",
        cres,
        conc[-1, :, :],
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
