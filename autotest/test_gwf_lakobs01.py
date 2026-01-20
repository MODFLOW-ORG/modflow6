"""
Test for checking lak observation input.  The following observation types:
'lak', 'wetted-area', and 'conductance,' require that ID2 be provided when
ID is an integer corresponding to a lake number and not BOUNDNAME.
See table in LAK Package section of mf6io.pdf for an explanation of ID,
ID2, and Observation Type.
"""

import os

import flopy
import numpy as np

cases = "gwf_lakobs_01a"
gwf = None


def get_idomain(nlay, nrow, ncol, lakend):
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    for k, j in enumerate(lakend):
        idomain[k, 0, 0:j] = 0

    return idomain


def build_model(dir, exe):
    lx = 300.0
    lz = 45.0
    nlay = 45
    nrow = 1
    ncol = 30
    nper = 1
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = 5.0
    botm = [top - (k + 1) * delz for k in range(nlay)]

    perlen = [20.0]
    nstp = [1]
    tsmult = [1.0]

    Kh = 1.0
    Kv = 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

    name = cases

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=exe,
        sim_ws=dir,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = name
    global gwf
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="NEWTON")

    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )

    # number of columns to be a lake for layer 1, 2, , ... len(lakend)
    lakend = [10, 9, 8, 7, 6]
    idomain = get_idomain(nlay, nrow, ncol, lakend)
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
    strt += top
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

    sy = 0.3
    ss = np.zeros((nlay, nrow, ncol), dtype=float)
    # ss[0, :, :] = sy
    idx = np.where(idomain == 0)
    for k, i, j in zip(*idx):
        ss[k + 1, i, j] = 0.0  # sy
    sto = flopy.mf6.ModflowGwfsto(gwf, sy=sy, ss=ss, iconvert=1)

    irch = np.zeros((nrow, ncol), dtype=int)
    lake_vconnect = []
    idx = np.where(idomain == 0)
    for k, i, j in zip(*idx):
        if idomain[k + 1, i, j] == 1:
            lake_vconnect.append((k + 1, i, j))
            irch[i, j] = k + 1
    nlakeconn = len(lake_vconnect)

    # pak_data = [ifno, strt, nlakeconn]
    initial_stage = 0.1
    pak_data = [(0, initial_stage, nlakeconn)]

    bedleak = 100.0  # "None"
    belev = 0.0
    con_data = [
        (0, i, idx, "VERTICAL", bedleak, belev, -99, -99, -99)
        for i, idx in enumerate(lake_vconnect)
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
            ("lak1", "lak", 1),
        ],
        "digits": 10,
    }

    lak = flopy.mf6.modflow.ModflowGwflak(
        gwf,
        surfdep=0.0,
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
    )

    chdspd = [((0, 0, ncol - 1), 5.0)]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=chdspd)

    rech = 0.0001 * np.ones((nrow, ncol), dtype=float)
    # rech[:, 0:20] = 0.
    rch = flopy.mf6.modflow.ModflowGwfrcha(
        gwf, print_flows=True, save_flows=True, recharge=rech, irch=irch
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def test_mf6model(function_tmpdir, targets):
    # build the models
    sim = build_model(str(function_tmpdir), targets["mf6"])

    # write model input
    sim.write_simulation()

    # attempt to run model, should fail
    sim.run_simulation()

    # ensure that the error msg is contained in the mfsim.lst file
    f = open(str(function_tmpdir / "mfsim.lst"), "r")
    lines = f.readlines()
    error_count = 0
    expected_msg = False
    for line in lines:
        if "ID2 (iconn) is missing" in line:
            expected_msg = True
            error_count += 1

    assert error_count == 1, "error count = " + str(error_count) + "but should equal 1"

    # fix the error and attempt to rerun model
    orig_fl = str(function_tmpdir / (cases + ".lak.obs"))
    new_fl = str(function_tmpdir / (cases + ".lak.obs.new"))
    sr = open(orig_fl, "r")
    sw = open(new_fl, "w")

    lines = sr.readlines()
    error_free_line = "  lak1  lak  1  1\n"
    for line in lines:
        if " lak " in line:
            sw.write(error_free_line)
        else:
            sw.write(line)

    sr.close()
    sw.close()

    # delete original and replace with corrected lab obs input
    os.remove(orig_fl)
    os.rename(new_fl, orig_fl)

    # rerun the model, should be no errors
    sim.run_simulation()
