"""
Test the time array series for the recharge package
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "utl06_tas_a",
    "utl06_tas_b",
    "utl06_tas_c",
    "utl06_tas_d",
]

nlay, nrow, ncol = 3, 5, 5
idomain_lay0 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 1, 1, 1],
]
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, :, :] = np.array(idomain_lay0)


def build_models(idx, test):
    perlen = [5.0]
    nstp = [5]
    tsmult = [1.0]
    nper = len(perlen)
    delr = [1.0, 2.0, 3.0, 4.0, 5.0]
    delc = [1.0, 2.0, 3.0, 4.0, 5.0]
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

    sim_name = "sim"

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf"
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
        idomain=idomain,
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
    spd = [[(nlay - 1, nrow - 1, ncol - 1), 4.0]]
    chd = flopy.mf6.modflow.ModflowGwfchd(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="CHD-1",
    )

    # assign recharge rate equal to one-based user node number
    recharge_rate = np.arange(nrow * ncol).reshape((nrow, ncol)) + 1
    if idx == 0:
        irch = None
        recharge = recharge_rate
    elif idx == 1:
        irch = [
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]
        recharge = recharge_rate
    elif idx == 2:
        irch = None
        recharge = "TIMEARRAYSERIES rcharray"
    elif idx == 3:
        irch = [
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]
        recharge = "TIMEARRAYSERIES rcharray"

    # array-based rch files
    rch4 = flopy.mf6.ModflowGwfrcha(
        gwf,
        print_flows=True,
        irch=irch,
        recharge=recharge,
        pname="RCH-4",
        filename=f"{gwfname}.rch4",
    )
    if idx in [2, 3]:
        filename = f"{gwfname}.rch4.tas"
        # for now write the recharge concentration to a dat file because there
        # is a bug in flopy that will not correctly write this array as internal
        tas_array = {0.0: 0.0, perlen[0]: f"{gwfname}.rch4.tas.dat"}
        time_series_namerecord = "rcharray"
        interpolation_methodrecord = "linear"
        rch4.tas.initialize(
            filename=filename,
            tas_array=tas_array,
            time_series_namerecord=time_series_namerecord,
            interpolation_methodrecord=interpolation_methodrecord,
        )
        np.savetxt(
            os.path.join(ws, f"{gwfname}.rch4.tas.dat"), recharge_rate, fmt="%7.1f"
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
    gwfname = "gwf"

    # load concentration file
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        head = hobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # load gwf budget file
    fpth = os.path.join(test.workspace, f"{gwfname}.cbc")
    try:
        bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    except:
        assert False, f'could not load data from "{fpth}"'

    rchbudall = bobj.get_data(text="RCH")
    times = bobj.get_times()

    print(times)
    for itime, totim in enumerate(times):
        print(f"Checking records for time {totim}")

        # Check records for each of the four recharge packages
        rchbud = rchbudall[itime]
        print(rchbud)
        print("  Checking records for recharge package")

        # id1 is the GWF user-based cell number
        print("    Checking id1")
        id1 = rchbud["node"]
        if idx in [0, 2]:
            id1a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
        elif idx in [1, 3]:
            id1a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                38,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
        assert np.allclose(id1, id1a), f"{id1} /= {id1a}"

        # id2 is the recharge package bound number
        print("    Checking id2")
        id2 = rchbud["node2"]
        if idx in [0, 2]:
            id2a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
        elif idx in [1, 3]:
            id2a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
        assert np.allclose(id2, id2a), f"{id2} /= {id2a}"

        print("    Checking q")
        q = rchbud["q"]
        if idx in [2, 3]:
            # time array series
            frac = (totim - 0.5) / 5.0
        else:
            frac = 1.0
        area = np.zeros((5, 5), dtype=float)
        for i in range(5):
            for j in range(5):
                area[i, j] = float((i + 1) * (j + 1))
        area = area.flatten()
        area = area[np.array(id2a) - 1]
        qa = [a * rate * frac for a, rate in zip(area, id2a)]
        assert np.allclose(q, qa), f"{q} /=\n {qa}"


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
