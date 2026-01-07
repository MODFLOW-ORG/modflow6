"""
Test first-order decay by running a one-cell model with ten 1-day time steps
with a decay rate of 1.  And a starting concentration of -10.  Results should
remain -10 throughout the simulation.
"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["src02_ha", "src02_hawd"]


def build_models(idx, test):
    nlay, nrow, ncol = 2, 1, 10
    nper = 1
    perlen = [10.0]
    nstp = [100]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = 10.0
    botm = [0.0, -1.0]
    laytyp = 1
    ss = 1e-5
    sy = 1.0
    strt = -0.9
    hk = 10.0
    recharge = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    if idx == 0:
        rewet_record = None
        wetdry = None
        gwf_linaccel = "bicgstab"
        newtonoptions = "newton under_relaxation"
    elif idx == 1:
        rewet_record = [("WETFCT", 1.0, "IWETIT", 1, "IHDWET", 1)]
        wetdry = -0.001  # [0.001, 0.001, 0.001]
        gwf_linaccel = "cg"
        newtonoptions = None

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
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
        newtonoptions=newtonoptions,
        save_flows=True,
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
        linear_acceleration=gwf_linaccel,
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
        idomain=1,
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
        rewet_record=rewet_record,
        wetdry=wetdry,
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=laytyp,
        sy=sy,
        ss=ss,
        transient={0: True},
        steady_state={0: False},
    )

    rch = flopy.mf6.ModflowGwfrch(gwf, stress_period_data=[(0, 0, 0, recharge)])

    drn = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=[(0, 0, ncol - 1, 0.1, 2.0)])

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(
        sim,
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
        save_flows=True,
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

    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=1)

    ssm = flopy.mf6.ModflowGwtssm(gwt)

    srcs = {0: [[(0, 0, 0), 1.00]]}
    src = flopy.mf6.ModflowGwtsrc(
        gwt,
        highest_saturated=True,
        maxbound=len(srcs),
        stress_period_data=srcs,
        save_flows=False,
        pname="SRC-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
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
    name = test.name
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name

    fpth = pl.Path(test.workspace) / f"{gwfname}.cbc"
    gwf_bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    fpth = pl.Path(test.workspace) / f"{gwtname}.cbc"
    gwt_bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    times = gwf_bobj.get_times()
    nodes = np.zeros(len(times), dtype=int)
    nodes_ans = np.zeros(len(times), dtype=int)
    saturation = np.zeros((len(times), 2), dtype=float)
    for idx, time in enumerate(times):
        src = gwt_bobj.get_data(text="src", totim=time)[0]
        nodes[idx] = src["node"][0] - 1
        sat = gwf_bobj.get_data(text="DATA-SAT", totim=time)[0]
        sat_on = sat[0]["sat"]
        saturation[idx] = sat_on
        if sat_on > 0.0:
            nodes_ans[idx] = 0
        else:
            nodes_ans[idx] = 10
    assert np.array_equal(nodes, nodes_ans), "src nodes not based on saturation"


@pytest.mark.developmode  # TODO remove for 6.7.0
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
