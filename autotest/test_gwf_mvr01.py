import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

name = "gwf_mvr01"
cases = [name]


def build_models(idx, test):
    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = []
    for _ in range(nper):
        tdis_rc.append((1.0, 1, 1.0))

    # spatial discretization data
    nlay, nrow, ncol = 3, 10, 10
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = [-10, -20, -30]
    strt = 0.0

    # calculate hk
    hk = 1.0e-4

    # solver options
    nouter, ninner = 600, 100
    hclose, rclose, relax = 1e-6, 0.1, 1.0
    newtonoptions = "NEWTON"
    imsla = "BICGSTAB"

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions=newtonoptions,
        save_flows=True,
        print_flows=True,
    )

    # create iterative model solution and register the gwf model with it
    csv0 = f"{name}.outer.ims.csv"
    csv1 = f"{name}.inner.ims.csv"
    ims = flopy.mf6.ModflowIms(
        sim,
        csv_outer_output_filerecord=csv0,
        csv_inner_output_filerecord=csv1,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk)

    # chd files
    # chd data
    spd = [
        [(0, 0, 0), 1.0],
        [(0, nrow - 1, ncol - 1), 0.0],
    ]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd, pname="chd-1")

    # drn file
    drn6 = [
        [(0, 1, 2), -1.0, 1.0],
        [(0, 2, 3), -1.0, 1.0],
    ]
    drn = flopy.mf6.modflow.ModflowGwfdrn(
        gwf, mover=True, stress_period_data=drn6, pname="drn-1"
    )

    # sfr file
    packagedata = [
        [
            0,
            (1 - 1, 4 - 1, 1 - 1),
            3.628e001,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            1.0e-4,
            1.0e-1,
            1,
            0.0,
            0,
        ],
        [
            1,
            (1 - 1, 4 - 1, 2 - 1),
            1.061e002,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            1.0e-4,
            1.0e-1,
            2,
            1.0,
            0,
        ],
        [
            2,
            (1 - 1, 4 - 1, 3 - 1),
            6.333e001,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            1.0e-4,
            1.0e-1,
            2,
            1.0,
            0,
        ],
        [
            3,
            (1 - 1, 5 - 1, 3 - 1),
            4.279e001,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            1.0e-4,
            1.0e-1,
            2,
            1.0,
            0,
        ],
        [
            4,
            (1 - 1, 5 - 1, 4 - 1),
            6.532e001,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            1.0e-4,
            1.0e-1,
            1,
            1.0,
            0,
        ],
    ]
    connectiondata = [
        [0, -1],
        [1, 0, -2],
        [2, 1, -3],
        [3, 2, -4],
        [4, 3],
    ]
    perioddata = [
        [0, "status", "active"],
        [1, "status", "active"],
        [2, "status", "active"],
        [3, "status", "active"],
        [4, "status", "active"],
    ]
    cnvgpth = f"{name}.sfr.cnvg.csv"
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        mover=True,
        nreaches=5,
        maximum_depth_change=1.0e-5,
        package_convergence_filerecord=cnvgpth,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )

    packagedata = [
        [0, 1.0, -20.0, 0.0, "SPECIFIED", 2],
    ]
    nmawwells = len(packagedata)
    connectiondata = [
        [1 - 1, 1 - 1, (1 - 1, 5 - 1, 8 - 1), 0.0, -20, 1.0, 1.1],
        [1 - 1, 2 - 1, (2 - 1, 5 - 1, 8 - 1), 0.0, -20, 1.0, 1.1],
    ]
    perioddata = [[0, "FLOWING_WELL", 0.0, 0.0, 0.0], [0, "RATE", 1.0e-3]]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        mover=True,
        nmawwells=nmawwells,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="maw-1",
    )

    packagedata = [(0, 1.0, 11), (1, 0.5, 11)]
    outlets = [
        (0, 0, 1, "manning", 0.001, 0.0, 0.1, 0.001),
        (1, 0, -1, "specified", -999.0, -999.0, -999.0, -999.0),
        (2, 0, -1, "specified", -999.0, -999.0, -999.0, -999.0),
    ]
    nlakes = len(packagedata)
    noutlets = len(outlets)
    connectiondata = [
        (0, 0, (0, 0, 5), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 1, (0, 1, 4), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 2, (1, 1, 5), "vertical", 1.0e-05, -5.0, 0.0, 1.0, 0.0),
        (0, 3, (0, 2, 4), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 4, (0, 3, 5), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 5, (0, 2, 6), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 6, (1, 2, 5), "vertical", 1.0e-05, -5.0, 0.0, 1.0, 0.0),
        (0, 7, (0, 0, 6), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 8, (0, 2, 6), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 9, (0, 1, 7), "horizontal", 1.0e-05, -5.0, 0.0, 100.0, 100.0),
        (0, 10, (1, 1, 6), "vertical", 1.0e-05, -5.0, 0.0, 1.0, 0.0),
        (1, 0, (0, 0, 8), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 1, (0, 1, 7), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 2, (0, 1, 9), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 3, (1, 1, 8), "vertical", 1.0e-05, -1.0, 0.0, 0.0, 0.0),
        (1, 4, (0, 2, 7), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 5, (0, 2, 9), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 6, (1, 2, 8), "vertical", 1.0e-05, -1.0, 0.0, 0.0, 0.0),
        (1, 7, (0, 3, 7), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 8, (0, 4, 8), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 9, (0, 3, 9), "horizontal", 1.0e-05, -1.0, 0.0, 100.0, 100.0),
        (1, 10, (1, 3, 8), "vertical", 1.0e-05, -1.0, 0.0, 0.0, 0.0),
    ]
    lakeperioddata = [
        (1, "status", "active"),
        (1, "rainfall", "0.0"),
        (1, "evaporation", "0.000000000000e+000"),
        (1, "runoff", "0.000000000000e+000"),
        (1, "withdrawal", "0.000000000000e+000"),
    ]
    outletperioddata = [
        (0, "rate", 1.000000000000e000),
        (0, "invert", 1.000000000000e-003),
        (0, "width", 0.000000000000e000),
        (0, "slope", 1.000000000000e-003),
        (0, "rough", 1.000000000000e-001),
        (1, "rate", -0.001),
        (2, "rate", 0.0),
    ]
    perioddata = lakeperioddata + outletperioddata
    cnvgpth = f"{name}.lak.cnvg.csv"
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        mover=True,
        nlakes=nlakes,
        noutlets=noutlets,
        print_stage=True,
        print_flows=True,
        package_convergence_filerecord=cnvgpth,
        packagedata=packagedata,
        connectiondata=connectiondata,
        outlets=outlets,
        perioddata=perioddata,
        pname="lak-1",
    )

    packagedata = [
        (0, (0, 5, 1), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (1, (0, 5, 2), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (2, (0, 5, 3), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (3, (0, 6, 1), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (4, (0, 6, 2), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (5, (0, 6, 3), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (6, (0, 7, 1), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (7, (0, 7, 2), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
        (8, (0, 7, 3), 1, -1, 1.0, 1.0e-05, 0.2, 0.4, 0.3, 3.5),
    ]
    perioddata = [
        [0, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [1, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [2, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [3, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [4, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [5, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [6, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [7, 1.0e-8, 0, 0, 0, 0, 0, 0],
        [8, 1.0e-8, 0, 0, 0, 0, 0, 0],
    ]
    cnvgpth = f"{name}.uzf.cnvg.csv"
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        mover=True,
        package_convergence_filerecord=cnvgpth,
        nuzfcells=len(packagedata),
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )

    packages = [("drn-1",), ("lak-1",), ("maw-1",), ("sfr-1",), ("uzf-1",)]
    perioddata = [
        ("drn-1", 0, "lak-1", 1, "excess", 1.0),
        ("drn-1", 0, "maw-1", 0, "threshold", 2.0),
        ("drn-1", 0, "sfr-1", 2, "upto", 3.0),
        ("drn-1", 1, "lak-1", 1, "excess", 1.0),
        ("drn-1", 1, "maw-1", 0, "threshold", 2.0),
        ("drn-1", 1, "sfr-1", 2, "upto", 3.0),
        ("lak-1", 0, "sfr-1", 0, "factor", 1.0),
        ("lak-1", 1, "sfr-1", 1, "factor", 0.5),
        ("lak-1", 2, "sfr-1", 2, "factor", 0.5),
        ("uzf-1", 0, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 1, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 2, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 3, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 4, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 5, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 6, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 7, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 8, "sfr-1", 0, "factor", 1.0),
    ]
    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=len(perioddata),
        budget_filerecord=f"{name}.mvr.bud",
        maxpackages=len(packages),
        print_flows=True,
        packages=packages,
        perioddata=perioddata,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    # mvr budget terms
    fpth = os.path.join(test.workspace, "gwf_mvr01.mvr.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    times = bobj.get_times()
    records = bobj.get_data(totim=times[-1])
    adt = [("node", "<i4"), ("node2", "<i4"), ("q", "<f8")]
    assert len(records) == 5

    assert records[0].shape == (2,)
    a = np.array([(1, 2, -0.0), (2, 2, -0.0)], dtype=adt)
    assert np.array_equal(records[0], a)

    assert records[1].shape == (2,)
    a = np.array([(1, 1, -0.0), (2, 1, -0.0)], dtype=adt)
    assert np.array_equal(records[1], a)

    assert records[2].shape == (2,)
    a = np.array([(1, 3, -0.00545875), (2, 3, -0.00468419)], dtype=adt)
    assert np.allclose(records[2]["node"], a["node"])
    assert np.allclose(records[2]["node2"], a["node2"])
    assert np.allclose(records[2]["q"], a["q"], atol=0.001), "{}\n{}".format(
        records[2]["q"], a["q"]
    )

    assert records[3].shape == (3,)
    a = np.array([(1, 1, -0.0), (1, 2, -0.0005), (1, 3, -0.0)], dtype=adt)
    assert np.array_equal(records[3], a)

    assert records[4].shape == (9,)
    a = np.array(
        [
            (1, 1, -1.0e-04),
            (2, 1, -1.0e-04),
            (3, 1, -1.0e-04),
            (4, 1, -1.0e-04),
            (5, 1, -1.0e-04),
            (6, 1, -1.0e-04),
            (7, 1, -1.0e-04),
            (8, 1, -1.0e-04),
            (9, 1, -1.0e-04),
        ],
        dtype=adt,
    )
    assert np.allclose(records[4]["node"], a["node"])
    assert np.allclose(records[4]["node2"], a["node2"])
    assert np.allclose(records[4]["q"], a["q"], atol=0.001), "{}\n{}".format(
        records[4]["q"], a["q"]
    )


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
