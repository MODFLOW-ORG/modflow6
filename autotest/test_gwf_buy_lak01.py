"""
Test the buoyancy package and the variable density flows between the lake
and the gwf model.  This model has 4 layers and a lake incised within it.
The model is transient and has heads in the aquifer higher than the initial
stage in the lake.  As the model runs, the lake and aquifer equalize and
should end up at the same level.  The test ensures that the initial and
final water volumes in the entire system are the same.  There are three
different cases:
 1.  No buoyancy package
 2.  Buoyancy package with lake and aquifer density = 1000.
 3.  Buoyancy package with lake and aquifer density = 1024.5
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["buy_lak_01a"]  # , 'buy_lak_01b', 'buy_lak_01c']
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

    nouter, ninner = 700, 300
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

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="NEWTON")

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

    idomain = np.full((nlay, nrow, ncol), 1)
    idomain[0, 0, 1:6] = 0
    idomain[1, 0, 2:5] = 0
    idomain[2, 0, 3:4] = 0
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # initial conditions
    strt = np.zeros((nlay, nrow, ncol), dtype=float)
    strt[0, 0, :] = 3.5
    strt[1, 0, :] = 3.0
    strt[1, 0, 1:6] = 2.5
    strt[2, 0, :] = 2.0
    strt[3, 0, :] = 1.0
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
    lake_dense = 1000.0 + 0.7 * c
    buy_on = buy_on_list[idx]
    if buy_on:
        pd = [(0, 0.7, 0.0, "none", "none")]
        buy = flopy.mf6.ModflowGwfbuy(
            gwf, packagedata=pd, denseref=1000.0, concentration=c
        )

    nlakeconn = 11  # note: number of connections for this lake
    # pak_data = [ifno, strt, nlakeconn, dense, boundname]
    pak_data = [(0, 2.25, nlakeconn, lake_dense)]

    connlen = delr / 2.0
    connwidth = delc
    bedleak = "None"
    con_data = [
        # con_data=(ifno,iconn,(cellid),claktype,bedleak,belev,telev,connlen,connwidth )
        (0, 0, (0, 0, 0), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
        (0, 1, (1, 0, 1), "VERTICAL", bedleak, 10, 10, connlen, connwidth),
        (0, 2, (1, 0, 1), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
        (0, 3, (2, 0, 2), "VERTICAL", bedleak, 10, 10, connlen, connwidth),
        (0, 4, (2, 0, 2), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
        (0, 5, (3, 0, 3), "VERTICAL", bedleak, 10, 10, connlen, connwidth),
        (0, 6, (2, 0, 4), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
        (0, 7, (2, 0, 4), "VERTICAL", bedleak, 10, 10, connlen, connwidth),
        (0, 8, (1, 0, 5), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
        (0, 9, (1, 0, 5), "VERTICAL", bedleak, 10, 10, connlen, connwidth),
        (0, 10, (0, 0, 6), "HORIZONTAL", bedleak, 10, 10, connlen, connwidth),
    ]

    # period data
    p_data = [
        (0, "STATUS", "ACTIVE"),
    ]

    # note: for specifying lake number, use fortran indexing!
    fname = f"{gwfname}.lak.obs.csv"
    lak_obs = {
        fname: [
            ("lakestage", "stage", 1),
            ("lakevolume", "volume", 1),
            ("lak1", "lak", 1, 1),
            ("lak2", "lak", 1, 2),
            ("lak3", "lak", 1, 3),
            ("lak4", "lak", 1, 4),
            ("lak5", "lak", 1, 5),
            ("lak6", "lak", 1, 6),
            ("lak7", "lak", 1, 7),
            ("lak8", "lak", 1, 8),
            ("lak9", "lak", 1, 9),
            ("lak10", "lak", 1, 10),
            ("lak11", "lak", 1, 11),
        ],
        # "digits": 10,
    }

    lak = flopy.mf6.modflow.ModflowGwflak(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        stage_filerecord=f"{gwfname}.lak.bin",
        budget_filerecord=f"{gwfname}.lak.bud",
        nlakes=len(pak_data),
        ntables=0,
        packagedata=pak_data,
        pname="LAK-1",
        connectiondata=con_data,
        perioddata=p_data,
        observations=lak_obs,
        auxiliary=["DENSITY"],
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    # calculate volume of water and make sure it is conserved
    gwfname = "gwf_" + test.name
    fname = gwfname + ".lak.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.HeadFile(fname, text="STAGE")
    stage = bobj.get_alldata().flatten()
    # print(stage)

    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    hobj = flopy.utils.HeadFile(fname)
    head = hobj.get_data()
    # print(head)

    # calculate initial water volume
    v0 = 3.5 * 2  # outermost columns
    v0 += 2.5 * 2  # next innermost columns
    v0 += 2.0 * 2  # next innermost columns
    v0 += 1.0 * 1  # middle column
    v0 = v0 * 0.3  # specific yield
    v0 = v0 + (2.25 - 2.0) * 2 + (2.25 - 1.0)
    print(f"initial volume of water in model = {v0}")

    # calculate ending water volume in model
    h = head[0, 0, 0]
    s = stage[-1]
    v = h * 4 + 2.0 * 2 + 1.0 * 1
    v = v * 0.3  # specific yield
    v = v + (s - 2.0) * 2 + (s - 1.0)
    print(f"final volume of water in model = {v}")

    # check to make sure starting water volume same as equalized final volume
    errmsg = f"initial and final water volume not equal: {v0} {v}"
    assert np.allclose(v0, v)

    # todo: add a better check of the lake concentrations
    # assert False


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
