"""
Tests wellbore conduction in the MWE Package with the MAW CONNECTION_STATUS
and STATUS period settings.

A single well is screened over all three layers of the model and does not pump,
so it exchanges heat with the aquifer only by conduction through the wetted area
of each connection. The well starts warmer than the aquifer.

Cases:
  - mwe_cs : the deepest connection is deactivated in period 2, the well is
             made INACTIVE in period 3 and ACTIVE again in period 4. A
             connection that is not simulated, because it or its well is
             inactive, must not conduct heat.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mwe_cs"]

nlay, nrow, ncol, nper = 3, 1, 5, 4
botm = [-10.0, -20.0, -30.0]
jcol = 2
deep = nlay - 1  # the deepest connection

# MAW period data; the well does not pump, so it exchanges heat by conduction
mawperioddata = {
    0: [[0, "rate", 0.0]],
    1: [[0, "connection_status", deep, "inactive"]],
    2: [[0, "status", "inactive"]],
    3: [[0, "status", "active"]],
}


def build_models(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"
    gwename = f"gwe-{name}"
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

    # energy transport model; the well starts warmer than the aquifer
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename)
    imsgwe = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-8,
        inner_dvclose=1e-9,
        linear_acceleration="bicgstab",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwename])
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=100.0,
        top=0.0,
        botm=botm,
    )
    flopy.mf6.ModflowGweic(gwe, strt=10.0)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.2,
        heat_capacity_water=4180.0,
        density_water=1000.0,
        heat_capacity_solid=840.0,
        density_solid=2500.0,
    )
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.27)
    flopy.mf6.ModflowGwessm(gwe)
    flopy.mf6.ModflowGwemwe(
        gwe,
        flow_package_name="MAW-1",
        save_flows=True,
        budget_filerecord=f"{gwename}.mwe.bud",
        packagedata=[(0, 40.0, 1.0, 0.1)],
        pname="MWE-1",
    )
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    fpth = os.path.join(test.workspace, f"gwe-{name}.mwe.bud")
    bud = flopy.utils.CellBudgetFile(fpth, precision="double")

    # conduction through each connection, for each stress period
    cond = []
    for kper in range(nper):
        rec = bud.get_data(kstpkper=(0, kper), text="WELLBORE-COND")[0]
        cond.append(np.array([r["q"] for r in rec]))

    # period 1: every connection of the warmer well conducts heat
    assert np.all(np.abs(cond[0]) > 0.0), f"all connections must conduct: {cond[0]}"
    # period 2: the inactive connection does not
    assert cond[1][deep] == 0.0, f"inactive connection conducted {cond[1][deep]}"
    assert abs(cond[1][0]) > 0.0, "active connections must still conduct"
    # period 3: an inactive well conducts no heat through any connection
    assert np.all(cond[2] == 0.0), f"inactive well conducted heat: {cond[2]}"
    # period 4: the reactivated well conducts again, except through the
    # connection that is still inactive
    assert abs(cond[3][0]) > 0.0, "reactivated well must conduct"
    assert cond[3][deep] == 0.0, f"inactive connection conducted {cond[3][deep]}"


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
