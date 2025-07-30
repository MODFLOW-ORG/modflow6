"""
Test the dependent_variable_scale option to make sure modflow fails if the option
is not specified for all gwt/gwe models in the same IMS solution. Tested for for a
one-dimensional model grid of square cells.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["dvscalefail"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 100
    nper = 1
    perlen = [5.0]
    nstp = [200]
    tsmult = [1.0]
    steady = [True]
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

    c = {0: [[(0, 0, 99), 0.0000000]]}
    w = {0: [[(0, 0, 0), 1.0, 1e6]]}

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

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
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    # create iterative model solution and register the gwf model with it
    gwf_ims = flopy.mf6.ModflowIms(
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
    sim.register_ims_package(gwf_ims, [gwf.name])

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
        icelltype=laytyp,
        k=hk,
        k33=hk,
        save_specific_discharge=True,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        pname="CHD-1",
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(w),
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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )

    # create gwt model
    gwt1name = "gwt1_" + name
    gwt1 = flopy.mf6.ModflowGwt(
        sim,
        dependent_variable_scaling=True,
        modelname=gwt1name,
        model_nam_file=f"{gwt1name}.nam",
        save_flows=True,
    )

    gwt2name = "gwt2_" + name
    gwt2 = flopy.mf6.ModflowGwt(
        sim,
        dependent_variable_scaling=False,
        modelname=gwt2name,
        model_nam_file=f"{gwt2name}.nam",
        save_flows=True,
    )

    # create iterative model solution and register the gwf model with it
    gwt_ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="bicgstab",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"gwt_{name}.ims",
    )
    sim.register_ims_package(gwt_ims, gwt1name)
    sim.register_ims_package(gwt_ims, gwt2name)

    for gwt in (gwt1, gwt2):
        gwtname = gwt.name

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
            saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("BUDGET", "LAST")],
        )

        # GWF GWT exchange
        gwfgwt = flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            filename=f"{gwtname}.gwfgwt",
            pname=gwtname,
        )

    return sim, None


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        xfail=True,
    )
    test.run()
