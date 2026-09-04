import os

import flopy
import numpy as np
import pytest
from flopy.utils.compare import eval_bud_diff
from framework import TestFramework

paktest = "sfr"
cases = ["ts_sfr01", "sfr_pkgdata_ts", "sfr_period_wins", "sfr_pd_ts_wins"]

# ---------------------------------------------------------------------------
# Shared data for SFR auxvar TS cases (cases 1-3)
# ---------------------------------------------------------------------------
_sfr_nper = 3
_sfr_period_data_aux = [(1.0, 1, 1.0)] * _sfr_nper
_sfr_ts_times = [0.0, 1.0, 2.0, 3.0]

# PKGDATA baseline aux values — change each period so any failure to update is visible
_sfr_pkgd_temp = [32.5, 40.0, 50.0]
_sfr_pkgd_conc = [0.1, 0.2, 0.3]

# PERIOD override values — distinct from pkgdata so a priority failure is visible
_sfr_per_temp = [20.0, 25.0, 15.0]
_sfr_per_conc = [0.5, 0.6, 0.7]

# TS arrays with a leading duplicate (LINEAREND gives the right value at period-end)
_sfr_ts_pkgd_temp = [_sfr_pkgd_temp[0]] + _sfr_pkgd_temp
_sfr_ts_pkgd_conc = [_sfr_pkgd_conc[0]] + _sfr_pkgd_conc
_sfr_ts_per_temp = [_sfr_per_temp[0]] + _sfr_per_temp
_sfr_ts_per_conc = [_sfr_per_conc[0]] + _sfr_per_conc


def _sfr_base_sim(ws, name):
    """Minimal 1-layer 1-row 5-col GWF sim for SFR auxvar TS tests."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=_sfr_nper, perioddata=_sfr_period_data_aux)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=5,
        delr=100.0,
        delc=100.0,
        top=10.0,
        botm=[0.0],
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 4), 4.0]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )
    return sim, gwf


def _sfr3_build(gwf, name, packagedata, perioddata, timeseries=None):
    """Attach a 3-reach SFR with aux=[temp,conc] to gwf."""
    connectiondata = [[0, -1], [1, 0, -2], [2, 1]]
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        auxiliary=["temp", "conc"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    if timeseries is not None:
        ts_names, ts_methods, ts_data = timeseries
        sfr.ts.initialize(
            filename=f"{name}.sfr.ts",
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=ts_methods,
        )
    return sfr


def _sfr3_pkgdata(temp, conc):
    """3-reach packagedata with literal temp/conc aux values."""
    rows = []
    for i in range(3):
        ncon = 1 if i in (0, 2) else 2
        rows.append(
            [
                i,
                (0, 0, i + 1),
                100.0,
                5.0,
                1e-3,
                4.0,
                1.0,
                1e-5,
                0.04,
                ncon,
                1.0,
                0,
                temp,
                conc,
            ]
        )
    return rows


def _sfr3_pkgdata_ts(temp_ts, conc_ts):
    """3-reach packagedata with TS name strings for temp/conc aux values."""
    rows = []
    for i in range(3):
        ncon = 1 if i in (0, 2) else 2
        rows.append(
            [
                i,
                (0, 0, i + 1),
                100.0,
                5.0,
                1e-3,
                4.0,
                1.0,
                1e-5,
                0.04,
                ncon,
                1.0,
                0,
                temp_ts,
                conc_ts,
            ]
        )
    return rows


# ---------------------------------------------------------------------------
# Case 1 — PKGDATA AUX TS only
#
# Reference: PERIOD AUXILIARY literal values change each period.
# TS model:  PKGDATA AUX TS drives aux; no PERIOD AUXILIARY override.
# Both produce identical GWF+SFR budgets (auxvars do not affect GWF flow).
# ---------------------------------------------------------------------------


def _get_sfr_pkgdata_ts_ref(ws, name):
    """Reference: PERIOD AUXILIARY literals change each period."""
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata(_sfr_pkgd_temp[0], _sfr_pkgd_conc[0])
    perioddata = {
        0: [
            [0, "inflow", 0.1],
            [0, "AUXILIARY", "temp", _sfr_pkgd_temp[0]],
            [0, "AUXILIARY", "conc", _sfr_pkgd_conc[0]],
        ],
        1: [
            [0, "AUXILIARY", "temp", _sfr_pkgd_temp[1]],
            [0, "AUXILIARY", "conc", _sfr_pkgd_conc[1]],
        ],
        2: [
            [0, "AUXILIARY", "temp", _sfr_pkgd_temp[2]],
            [0, "AUXILIARY", "conc", _sfr_pkgd_conc[2]],
        ],
    }
    _sfr3_build(gwf, name, pkgdata, perioddata)
    return sim


def _get_sfr_pkgdata_ts(ws, name):
    """TS: PKGDATA AUX TS drives aux; no PERIOD AUXILIARY."""
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata_ts("pkgd_temp", "pkgd_conc")
    perioddata = {0: [[0, "inflow", 0.1]]}
    ts_data = list(zip(_sfr_ts_times, _sfr_ts_pkgd_temp, _sfr_ts_pkgd_conc))
    _sfr3_build(
        gwf,
        name,
        pkgdata,
        perioddata,
        timeseries=(["pkgd_temp", "pkgd_conc"], ["linearend", "linearend"], ts_data),
    )
    return sim


# ---------------------------------------------------------------------------
# Case 2 — PERIOD override wins over PKGDATA AUX TS
#
# Reference: PERIOD AUXILIARY literal values (_sfr_per_*); no PKGDATA TS.
# TS model:  PKGDATA AUX TS (_sfr_pkgd_*) + PERIOD AUXILIARY (_sfr_per_*).
# PERIOD must win → same auxvars → identical budgets.
# ---------------------------------------------------------------------------


def _get_sfr_period_wins_ref(ws, name):
    """Reference: PERIOD AUXILIARY literal (per) values; PKGDATA literal."""
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata(_sfr_per_temp[0], _sfr_per_conc[0])
    perioddata = {
        0: [
            [0, "inflow", 0.1],
            [0, "AUXILIARY", "temp", _sfr_per_temp[0]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[0]],
        ],
        1: [
            [0, "AUXILIARY", "temp", _sfr_per_temp[1]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[1]],
        ],
        2: [
            [0, "AUXILIARY", "temp", _sfr_per_temp[2]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[2]],
        ],
    }
    _sfr3_build(gwf, name, pkgdata, perioddata)
    return sim


def _get_sfr_period_wins_ts(ws, name):
    """TS: PKGDATA AUX TS (pkgd values) + PERIOD AUXILIARY override (per values).

    PERIOD must win: auxvar should equal _sfr_per_* each period.
    """
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata_ts("pkgd_temp", "pkgd_conc")
    perioddata = {
        0: [
            [0, "inflow", 0.1],
            [0, "AUXILIARY", "temp", _sfr_per_temp[0]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[0]],
        ],
        1: [
            [0, "AUXILIARY", "temp", _sfr_per_temp[1]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[1]],
        ],
        2: [
            [0, "AUXILIARY", "temp", _sfr_per_temp[2]],
            [0, "AUXILIARY", "conc", _sfr_per_conc[2]],
        ],
    }
    ts_data = list(zip(_sfr_ts_times, _sfr_ts_pkgd_temp, _sfr_ts_pkgd_conc))
    _sfr3_build(
        gwf,
        name,
        pkgdata,
        perioddata,
        timeseries=(["pkgd_temp", "pkgd_conc"], ["linearend", "linearend"], ts_data),
    )
    return sim


# ---------------------------------------------------------------------------
# Case 3 — PERIOD AUX TS wins over PKGDATA AUX TS
#
# Reference: PKGDATA literal + PERIOD AUXILIARY TS (per_ts values).
# TS model:  PKGDATA AUX TS (pkgd values) + PERIOD AUXILIARY TS (per_ts values).
# PERIOD TS must win → same auxvars → identical budgets.
# ---------------------------------------------------------------------------


def _get_sfr_pd_ts_wins_ref(ws, name):
    """Reference: PKGDATA literal + PERIOD AUXILIARY TS (per values)."""
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata(_sfr_per_temp[0], _sfr_per_conc[0])
    perioddata = {
        0: [
            [0, "inflow", 0.1],
            [0, "AUXILIARY", "temp", "per_temp"],
            [0, "AUXILIARY", "conc", "per_conc"],
        ],
    }
    ts_data = list(zip(_sfr_ts_times, _sfr_ts_per_temp, _sfr_ts_per_conc))
    _sfr3_build(
        gwf,
        name,
        pkgdata,
        perioddata,
        timeseries=(["per_temp", "per_conc"], ["linearend", "linearend"], ts_data),
    )
    return sim


def _get_sfr_pd_ts_wins_ts(ws, name):
    """TS: PKGDATA AUX TS (pkgd values) + PERIOD AUXILIARY TS (per values).

    PERIOD TS must win: auxvar should equal _sfr_per_* values.
    Both series are written to a single SFR TS file.
    """
    sim, gwf = _sfr_base_sim(ws, name)
    pkgdata = _sfr3_pkgdata_ts("pkgd_temp", "pkgd_conc")
    perioddata = {
        0: [
            [0, "inflow", 0.1],
            [0, "AUXILIARY", "temp", "per_temp"],
            [0, "AUXILIARY", "conc", "per_conc"],
        ],
    }
    ts_data = [
        (t, pt, pc, tt, tc)
        for t, pt, pc, tt, tc in zip(
            _sfr_ts_times,
            _sfr_ts_pkgd_temp,
            _sfr_ts_pkgd_conc,
            _sfr_ts_per_temp,
            _sfr_ts_per_conc,
        )
    ]
    _sfr3_build(
        gwf,
        name,
        pkgdata,
        perioddata,
        timeseries=(
            ["pkgd_temp", "pkgd_conc", "per_temp", "per_conc"],
            ["linearend"] * 4,
            ts_data,
        ),
    )
    return sim


def get_model(ws, name, timeseries=False):
    # static model data
    # temporal discretization
    nper = 3
    tdis_rc = [(1.0, 1, 1.0), (1.0, 10, 1.0), (1.0, 10, 1.0)]
    ts_times = np.arange(0.0, float(nper) + 1.0, 1.0, dtype=float)

    auxnames = ["temp", "conc"]
    temp, conc = 32.5, 0.1

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
        sim_name=name,
        memory_print_option="all",
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    # set ims csv files
    csv0 = f"{name}.outer.ims.csv"
    csv1 = f"{name}.inner.ims.csv"

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=csv0,
        csv_inner_output_filerecord=csv1,
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

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions=newtonoptions,
        print_input=True,
        save_flows=True,
        print_flows=True,
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
            2,
            0.0,
            1,
            temp,
            conc,
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
            3,
            1.0,
            1,
            temp,
            conc,
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
            4,
            1.0,
            2,
            temp,
            conc,
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
            3,
            1.0,
            1,
            temp,
            conc,
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
            temp,
            conc,
        ],
        [
            5,
            (1 - 1, 4 - 1, 1 - 1),
            10.0,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            0.0,
            1.0e-1,
            1,
            0.0,
            0,
            temp,
            conc,
        ],
        [
            6,
            (1 - 1, 4 - 1, 2 - 1),
            10.0,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            0.0,
            1.0e-1,
            1,
            0.0,
            0,
            temp,
            conc,
        ],
        [
            7,
            (1 - 1, 4 - 1, 3 - 1),
            10.0,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            0.0,
            1.0e-1,
            1,
            0.0,
            0,
            temp,
            conc,
        ],
        [
            8,
            (1 - 1, 4 - 1, 3 - 1),
            10.0,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            0.0,
            1.0e-1,
            1,
            0.0,
            0,
            temp,
            conc,
        ],
        [
            9,
            (1 - 1, 5 - 1, 4 - 1),
            10.0,
            1.0,
            1.0e-003,
            0.0,
            1.0,
            0.0,
            1.0e-1,
            1,
            0.0,
            0,
            temp,
            conc,
        ],
    ]
    connectiondata = [
        [0, -1, -5],
        [1, 0, -2, -6],
        [2, -3, -7, -8, 1],
        [3, -4, -9, 2],
        [4, 3],
        [5, 0],
        [6, 1],
        [7, 2],
        [8, 2],
        [9, 3],
    ]
    cprior = "upto"
    divdata = [
        [0, 0, 5, cprior],
        [1, 0, 6, cprior],
        [2, 1, 7, cprior],
        [2, 0, 8, cprior],
        [3, 0, 9, cprior],
    ]
    inflow, divflow, divflow2, upstream_fraction = 1.0, 0.05, 0.04, 0.0
    ts_names = ["inflow", "divflow", "ustrf"] + auxnames
    perioddata = [
        [0, "status", "active"],
        [1, "status", "active"],
        [2, "status", "active"],
        [3, "status", "active"],
        [4, "status", "active"],
        [0, "diversion", 0, divflow],
        [1, "diversion", 0, divflow],
        [2, "diversion", 0, divflow2],
        [3, "diversion", 0, divflow],
    ]
    if timeseries:
        perioddata.append([0, "inflow", "inflow"])
        perioddata.append([2, "diversion", 1, "divflow"])
        perioddata.append([0, "AUXILIARY", "conc", "conc"])
        perioddata.append([2, "AUXILIARY", "temp", "temp"])
        perioddata.append([5, "upstream_fraction", "ustrf"])
        perioddata.append([7, "upstream_fraction", "ustrf"])
        perioddata.append([9, "upstream_fraction", "ustrf"])
        ts_methods = ["linearend"] * len(ts_names)
        ts_data = []
        for t in ts_times:
            ts_data.append((t, inflow, divflow, upstream_fraction, temp, conc))
    else:
        perioddata.append([0, "inflow", inflow])
        perioddata.append([2, "diversion", 1, divflow])

    budpth = f"{name}.{paktest}.cbc"
    cnvgpth = f"{name}.sfr.cnvg.csv"
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        maximum_picard_iterations=1,
        auxiliary=auxnames,
        print_input=True,
        budget_filerecord=budpth,
        mover=True,
        nreaches=len(packagedata),
        maximum_depth_change=1.0e-5,
        package_convergence_filerecord=cnvgpth,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=divdata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    if timeseries:
        fname = f"{name}.sfr.ts"
        sfr.ts.initialize(
            filename=fname,
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=ts_methods,
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
        print_head=True,
        mover=True,
        nmawwells=nmawwells,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="maw-1",
    )

    packagedata = [(0, 1.0, 11), (1, 0.5, 11)]
    outlets = [(0, 0, 1, "manning", 0.001, 0.0, 0.1, 0.001)]
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
    perioddata = [
        (1, "status", "active"),
        (1, "rainfall", "0.0"),
        (1, "evaporation", "0.000000000000e+000"),
        (1, "runoff", "0.000000000000e+000"),
        (1, "withdrawal", "0.000000000000e+000"),
        (0, "rate", "1.000000000000e+000"),
        (0, "invert", "1.000000000000e-003"),
        (0, "width", "0.000000000000e+000"),
        (0, "slope", "1.000000000000e-003"),
        (0, "rough", "1.000000000000e-001"),
    ]
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
        ("uzf-1", 0, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 1, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 2, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 3, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 4, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 5, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 6, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 7, "sfr-1", 0, "factor", 1.0),
        ("uzf-1", 8, "sfr-1", 0, "factor", 1.0),
        ("sfr-1", 2, "sfr-1", 3, "factor", 0.5),
        ("sfr-1", 6, "sfr-1", 4, "factor", 0.5),
        ("sfr-1", 8, "sfr-1", 4, "factor", 0.5),
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
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "LAST"), ("HEAD", "LAST")],
    )

    return sim


def build_models(idx, test):
    name = cases[idx]
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")

    if idx == 0:
        return get_model(ws0, name), get_model(ws1, name, timeseries=True)
    elif idx == 1:
        return _get_sfr_pkgdata_ts_ref(ws0, name), _get_sfr_pkgdata_ts(ws1, name)
    elif idx == 2:
        return _get_sfr_period_wins_ref(ws0, name), _get_sfr_period_wins_ts(ws1, name)
    else:
        return _get_sfr_pd_ts_wins_ref(ws0, name), _get_sfr_pd_ts_wins_ts(ws1, name)


def _check_sfr_aux(cobj, expected_temp, expected_conc):
    """Assert SFR's AUXILIARY budget term's per-reach temp/conc equal
    expected_temp/expected_conc (shape (nper, nreaches)) each period."""
    aux = cobj.get_data(text="AUXILIARY")
    for iper in range(_sfr_nper):
        temp = aux[iper]["TEMP"]
        conc = aux[iper]["CONC"]
        assert np.allclose(temp, expected_temp[iper]), (
            f"period {iper}: SFR AUX temp {temp} != expected {expected_temp[iper]}"
        )
        assert np.allclose(conc, expected_conc[iper]), (
            f"period {iper}: SFR AUX conc {conc} != expected {expected_conc[iper]}"
        )


def check_result(idx, test):
    name = os.path.basename(test.name)
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")

    # GWF budget comparison (all cases)
    ia = (
        flopy.mf6.utils.MfGrdFile(os.path.join(ws0, f"{name}.dis.grb"))._datadict["IA"]
        - 1
    )
    cobj0 = flopy.utils.CellBudgetFile(
        os.path.join(ws0, f"{name}.cbc"), precision="double"
    )
    cobj1 = flopy.utils.CellBudgetFile(
        os.path.join(ws1, f"{name}.cbc"), precision="double"
    )
    eval_bud_diff(os.path.join(ws0, f"{name}.cbc.cmp.out"), cobj0, cobj1, ia)

    # SFR package budget comparison (all cases)
    sfr_cbc = f"{name}.{paktest}.cbc"
    cobj0 = flopy.utils.CellBudgetFile(os.path.join(ws0, sfr_cbc), precision="double")
    cobj1 = flopy.utils.CellBudgetFile(os.path.join(ws1, sfr_cbc), precision="double")
    eval_bud_diff(os.path.join(ws0, f"{name}.{paktest}.cbc.cmp.out"), cobj0, cobj1)

    if idx == 0:
        # Spot checks specific to the ts_sfr01 geometry
        v0 = cobj0.get_data(totim=1.0, text="FLOW-JA-FACE")[0]
        q = [v0["q"][i] for i, node in enumerate(v0["node"]) if node > 5]
        v0 = np.array(q)
        check = np.ones(v0.shape, dtype=float) * 5e-2
        check[-2] = 4e-2
        assert np.allclose(v0, check), "FLOW-JA-FACE failed"

        v0 = cobj0.get_data(totim=1.0, text="EXT-OUTFLOW")[0]
        v0 = v0["q"][4:]
        check = np.array([-0.80871, -5e-2, -2.5e-2, -5e-2, -2.0e-2, -5e-2])
        assert np.allclose(v0, check), "EXT-OUTFLOW failed"

        v0 = cobj0.get_data(totim=1.0, text="FROM-MVR")[0]
        v0 = v0["q"][4:]
        check = np.array([4.5e-2, 0.0, 0.0, 0.0, 0.0, 0.0])
        assert np.allclose(v0, check), "FROM-MVR failed"

        v0 = cobj0.get_data(totim=1.0, text="TO-MVR")[0]
        v0 = v0["q"][4:]
        check = np.array([0.0, 0.0, -2.5e-2, 0.0, -2.0e-2, 0.0])
        assert np.allclose(v0, check), "TO-MVR failed"
    elif idx == 1:
        # PKGDATA AUX TS drives all reaches; no PERIOD override
        expected_temp = [[t] * 3 for t in _sfr_pkgd_temp]
        expected_conc = [[c] * 3 for c in _sfr_pkgd_conc]
        _check_sfr_aux(cobj1, expected_temp, expected_conc)
    else:
        # idx 2/3: PERIOD AUXILIARY (literal or TS) overrides reach 0 only;
        # reaches 1/2 stay driven by the PKGDATA AUX TS
        expected_temp = [[t, p, p] for t, p in zip(_sfr_per_temp, _sfr_pkgd_temp)]
        expected_conc = [[c, p, p] for c, p in zip(_sfr_per_conc, _sfr_pkgd_conc)]
        _check_sfr_aux(cobj1, expected_temp, expected_conc)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_result(idx, t),
        targets=targets,
    )
    test.run()
