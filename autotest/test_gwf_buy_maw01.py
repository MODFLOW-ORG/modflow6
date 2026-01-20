"""
Test the buoyancy package and the variable density flows between maw
and the gwf model.  This model has 4 layers with a single maw.
The model is transient and has heads in the aquifer higher than the initial
stage in the well.  As the model runs, the well and aquifer equalize and
should end up at the same level.  The test ensures that the initial and
final water volumes in the entire system are the same.  There are three
different cases:
 1.  No buoyancy package
 2.  Buoyancy package with maw and aquifer density = 1000.
 3.  Buoyancy package with maw and aquifer density = 1024.5
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["buy_maw_01a"]  # , 'buy_maw_01b', 'buy_maw_01c']
buy_on_list = [False]  # , True, True]
concbuylist = [0.0]  # , 0., 35.]


def build_models(idx, test):
    lx = 7.0
    lz = 4.0
    nlay = 4
    nrow = 1
    ncol = 7
    nper = 1
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = 4.0
    botm = [3.0, 2.0, 1.0, 0.0]

    perlen = [10.0]
    nstp = [50]
    tsmult = [1.0]

    Kh = 1.0
    Kv = 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    nouter, ninner = 700, 10
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

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

    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions=newtonoptions)

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
    strt = 4.0
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=Kh,
        k33=Kv,
    )

    sto = flopy.mf6.ModflowGwfsto(gwf, sy=0.3, ss=0.0, iconvert=1)

    c = concbuylist[idx]
    mawdense = 1000.0 + 0.7 * c
    buy_on = buy_on_list[idx]
    if buy_on:
        pd = [(0, 0.7, 0.0, "none", "none")]
        buy = flopy.mf6.ModflowGwfbuy(
            gwf, packagedata=pd, denseref=1000.0, concentration=c
        )

    mawradius = 0.1
    mawbottom = 0.0
    mawstrt = 3.5
    mawcondeqn = "THIEM"
    mawngwfnodes = nlay
    # <ifno> <radius> <bottom> <strt> <condeqn> <ngwfnodes>
    mawpackagedata = [
        [0, mawradius, mawbottom, mawstrt, mawcondeqn, mawngwfnodes, mawdense]
    ]
    # <ifno> <icon> <cellid(ncelldim)> <scrn_top> <scrn_bot> <hk_skin> <radius_skin>
    mawconnectiondata = [
        [0, icon, (icon, 0, 0), top, mawbottom, -999.0, -999.0] for icon in range(nlay)
    ]
    # <ifno> <mawsetting>
    mawperioddata = [[0, "STATUS", "ACTIVE"]]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        head_filerecord=f"{gwfname}.maw.bin",
        budget_filerecord=f"{gwfname}.maw.bud",
        packagedata=mawpackagedata,
        connectiondata=mawconnectiondata,
        perioddata=mawperioddata,
        auxiliary=["DENSITY"],
        pname="MAW-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    return sim, None


def check_output(idx, test):
    # calculate volume of water and make sure it is conserved
    gwfname = "gwf_" + test.name
    fname = gwfname + ".maw.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.HeadFile(fname, text="HEAD")
    stage = bobj.get_alldata().flatten()

    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    hobj = flopy.utils.HeadFile(fname)
    head = hobj.get_alldata()

    # calculate initial volume of water in well and aquifer
    v0maw = 3.5 * np.pi * 0.1**2
    v0gwf = 4 * 7 * 0.3
    v0 = v0maw + v0gwf
    top = [4.0, 3.0, 2.0, 1.0]
    botm = [3.0, 2.0, 1.0, 0.0]
    nlay = 4
    ncol = 7

    # calculate current volume of water in well and aquifer and compare with
    # initial volume
    for kstp, mawstage in enumerate(stage):
        vgwf = 0
        for k in range(nlay):
            for j in range(ncol):
                tp = min(head[kstp, k, 0, j], top[k])
                dz = tp - botm[k]
                vgwf += 0.3 * max(0.0, dz)
        vmaw = stage[kstp] * np.pi * 0.1**2
        vnow = vmaw + vgwf
        errmsg = f"kstp {kstp}: current volume ({v0}) not equal initial volume ({vnow})"
        assert np.allclose(v0, vnow), errmsg

    # compare the maw-gwf flows in maw budget file with the gwf-maw flows in
    # gwf budget file.  Values should be the same but reversed in sign
    fname = gwfname + ".maw.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    mbud = flopy.utils.CellBudgetFile(fname, precision="double")
    maw_gwf = mbud.get_data(text="GWF")

    fname = gwfname + ".cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    gbud = flopy.utils.CellBudgetFile(fname, precision="double")
    gwf_maw = gbud.get_data(text="MAW")

    assert len(maw_gwf) == len(gwf_maw), "number of budget records not equal"

    for istp, (ra_maw, ra_gwf) in enumerate(zip(maw_gwf, gwf_maw)):
        for i in range(ra_maw.shape[0]):
            qmaw = ra_maw[i]["q"]
            qgwf = ra_gwf[i]["q"]
            msg = f"step {istp} record {i} comparing qmaw with qgwf: {qmaw} {qgwf}"
            print(msg)
            assert np.allclose(qmaw, -qgwf), msg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
