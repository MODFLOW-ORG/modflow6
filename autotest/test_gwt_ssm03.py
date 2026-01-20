"""
Test the SSM FILEINPUT option for specifying source and sink
concentrations.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ssm03"]


def build_models(idx, test):
    nlay, nrow, ncol = 3, 5, 5
    perlen = [5.0]
    nstp = [5]
    tsmult = [1.0]
    nper = len(perlen)
    delr = 1.0
    delc = 1.0
    top = 4.0
    botm = [3.0, 2.0, 1.0]
    strt = 4.0
    hk = 1.0
    laytyp = 0

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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=laytyp,
        k=hk,
        save_specific_discharge=True,
    )

    # chd files
    spd = [[(0, 0, 0), 4.0]]
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="CHD-1",
    )

    # wel files
    spd = [[(0, nrow - 1, ncol - 1), 1.0]]
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="WEL-1",
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
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt)

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    pd = [(0, "concentration", 100.0)]
    spc = flopy.mf6.ModflowUtlspc(
        gwt, perioddata=pd, maxbound=len(pd), filename=f"{gwtname}.wel1.spc"
    )
    sourcerecarray = [()]
    fileinput = [
        ("WEL-1", f"{gwtname}.wel1.spc"),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, print_flows=True, sources=sourcerecarray, fileinput=fileinput
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        f"{gwtname}.obs.csv": [
            ("(1-1-1)", "CONCENTRATION", (0, 0, 0)),
            ("(1-5-5)", "CONCENTRATION", (0, 4, 4)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
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

    # load concentration file
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()

    # load transport budget file
    fpth = os.path.join(test.workspace, f"{gwtname}.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    ssmbudall = bobj.get_data(text="SOURCE-SINK MIX")
    for ssmbud in ssmbudall:
        node, node2, q = ssmbud[0]
        assert node == 25, "node location for well must be 25 (GWT cell 25)"
        assert node2 == 1, "node2 location for well must be 1 (first well)"
        assert q == 100.0, "mass flux for well must be 100."

        node, node2, q = ssmbud[1]
        assert node == 1, "node location for chd must be 1 (first GWT cell)"
        assert node2 == 1, "node2 location for chd must be 1 (first chd)"
        assert q < 0.0, "mass flux for chd must be less than zero"


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
