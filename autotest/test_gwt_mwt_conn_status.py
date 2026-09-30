"""
Tests the MWT Package with the MAW CONNECTION_STATUS and STATUS period settings.

A single pumping well is screened over all three layers of the model, and only
the deepest layer holds solute, so the well draws mass only through its deepest
connection.

Cases:
  - mwt_cs : the deepest connection is deactivated in period 2, the well is made
             INACTIVE in period 3 and ACTIVE again in period 4. A connection that
             is not simulated, because it or its well is inactive, must exchange
             no mass with the aquifer, and the mass budget must close.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mwt_cs"]

nlay, nrow, ncol, nper = 3, 1, 5, 4
botm = [-10.0, -20.0, -30.0]
jcol = 2
deep = nlay - 1  # the deepest connection
mawrate = -100.0

mawperioddata = {
    0: [[0, "rate", mawrate]],
    1: [[0, "connection_status", deep, "inactive"]],
    2: [[0, "status", "inactive"]],
    3: [[0, "status", "active"]],
}


def build_models(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"
    gwtname = f"gwt-{name}"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(1.0, 1, 1.0)] * nper
    )

    # flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=100.0,
        top=0.0,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=-1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=10.0, k33=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), -1.0], [(0, 0, ncol - 1), -1.0]]
    )
    conndata = [
        [0, k, (k, 0, jcol), 0.0 if k == 0 else botm[k - 1], botm[k], 1.0, 0.1]
        for k in range(nlay)
    ]
    flopy.mf6.ModflowGwfmaw(
        gwf,
        save_flows=True,
        packagedata=[[0, 0.15, botm[-1], -1.0, "THIEM", nlay]],
        connectiondata=conndata,
        perioddata=mawperioddata,
        pname="MAW-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # transport model; only the deepest layer holds solute
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-8,
        inner_dvclose=1e-9,
        linear_acceleration="bicgstab",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=100.0,
        top=0.0,
        botm=botm,
    )
    strt = np.zeros((nlay, nrow, ncol))
    strt[deep] = 100.0
    flopy.mf6.ModflowGwtic(gwt, strt=strt)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.2)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtmwt(
        gwt,
        flow_package_name="MAW-1",
        save_flows=True,
        print_flows=True,
        budget_filerecord=f"{gwtname}.mwt.bud",
        packagedata=[(0, 0.0)],
        pname="MWT-1",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    fpth = os.path.join(test.workspace, f"gwt-{name}.mwt.bud")
    bud = flopy.utils.CellBudgetFile(fpth, precision="double")

    # mass exchanged with the aquifer through each connection, for each period
    mass = []
    for kper in range(nper):
        rec = bud.get_data(kstpkper=(0, kper), text="GWF")[0]
        mass.append(np.array([r["q"] for r in rec]))

    # period 1: the well draws solute through its deepest connection
    assert abs(mass[0][deep]) > 0.0, "the deepest connection must draw solute"
    # period 2: the inactive connection exchanges no mass
    assert mass[1][deep] == 0.0, f"inactive connection exchanged {mass[1][deep]}"
    # period 3: an inactive well exchanges no mass through any connection
    assert np.all(mass[2] == 0.0), f"inactive well exchanged mass: {mass[2]}"
    # period 4: the connection is still inactive after the well is reactivated
    assert mass[3][deep] == 0.0, f"inactive connection exchanged {mass[3][deep]}"

    # the transport mass budget closes in every period
    lst = flopy.utils.Mf6ListBudget(
        os.path.join(test.workspace, f"gwt-{name}.lst"),
        budgetkey="MASS BUDGET FOR ENTIRE MODEL",
    )
    pdiff = lst.get_incremental()["PERCENT_DISCREPANCY"]
    assert np.all(np.abs(pdiff) < 1e-3), f"mass budget does not close: {pdiff}"


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
