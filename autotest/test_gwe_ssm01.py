"""
Test the SSM FILEINPUT option for specifying source and sink
temperatures.  Similar to test_gwt_ssm03.py
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ssm01"]

Cpw = 4184.0
rhow = 1000.0


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
    Cps = 703.7
    rhos = 2700.0

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
    gwename = "gwe_" + name
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    gwe.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
    imsgwe = flopy.mf6.ModflowIms(
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
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    dis = flopy.mf6.ModflowGwedis(
        gwe,
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
    ic = flopy.mf6.ModflowGweic(gwe, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGweadv(gwe)

    # mass storage and transfer
    est = flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.1,
        heat_capacity_solid=Cps,
        density_solid=rhos,
    )

    # sources
    pd = [(0, "temperature", 100.0)]
    spc = flopy.mf6.ModflowUtlspc(
        gwe, perioddata=pd, maxbound=len(pd), filename=f"{gwename}.wel1.spc"
    )
    sourcerecarray = [()]
    fileinput = [
        ("WEL-1", f"{gwename}.wel1.spc"),
    ]
    ssm = flopy.mf6.ModflowGwessm(
        gwe, print_flows=True, sources=sourcerecarray, fileinput=fileinput
    )

    # output control
    oc = flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        f"{gwename}.obs.csv": [
            ("(1-1-1)", "TEMPERATURE", (0, 0, 0)),
            ("(1-5-5)", "TEMPERATURE", (0, 4, 4)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwe,
        pname=f"{gwename}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    # GWF GWE exchange
    gwfgwe = flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )

    return sim, None


def check_output(idx, test):
    name = test.name
    gwename = "gwe_" + name

    # load temperature file
    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
    temp = tobj.get_data()

    # load transport budget file
    fpth = os.path.join(test.workspace, f"{gwename}.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    ssmbudall = bobj.get_data(text="SOURCE-SINK MIX")
    for ssmbud in ssmbudall:
        node, node2, q = ssmbud[0]
        assert node == 25, "node location for well must be 25 (GWE cell 25)"
        assert node2 == 1, "node2 location for well must be 1 (first well)"
        msg0 = (
            "energy flux for well must the input temperature (100 deg C) "
            "multiplied the heat capacity and density of water."
        )
        assert q == 100.0 * Cpw * rhow, msg0

        node, node2, q = ssmbud[1]
        assert node == 1, "node location for chd must be 1 (first GWE cell)"
        assert node2 == 1, "node2 location for chd must be 1 (first chd)"
        assert q < 0.0, "energy flux for chd must be less than zero"


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
