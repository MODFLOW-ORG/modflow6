"""PERIOD RATE/WELL_HEAD/AUXILIARY time-series persistence for GWF-MAW,
across three related scenarios:

  continue  (maw_ratecont, maw_wellheadcont, maw_auxcont):
      a value linked to a TS in period 0 keeps tracking that TS in later
      periods whose own PERIOD block doesn't reissue it.
  switch    (maw_rateswitch, maw_whswitch):
      a single well's setting is a literal, then a TS, then a new
      literal, then unreissued -- checks that a literal override clears
      the earlier TS link.
  partial   (maw_ratepartial, maw_whpartial):
      two wells; well 0 switches to a TS partway through, well 1's
      setting is set once in period 0 and never reissued again -- checks
      a partial PERIOD reissue only touches the well it names.

Model geometry: 1-layer, 1-row, single-well cases use 3 columns (CHD at
    columns 0 and 2); the two-well "partial" case uses 4 columns (CHD at
    columns 0 and 3, wells in the middle two). All steady-state, 1
    timestep per period.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "maw"
cases = [
    "maw_ratecont",
    "maw_wellheadcont",
    "maw_auxcont",
    "maw_rateswitch",
    "maw_whswitch",
    "maw_ratepartial",
    "maw_whpartial",
]

radius = 0.1
well_bot = -10.0
strt = 90.0
top = 200.0
botm = [0.0]
k11 = 1.0
chd_head = 100.0
delr = delc = 100.0

# -- continue group --
CONT_NPER = 3
CONT_TS_TIMES = [0.0, 1.0, 2.0, 3.0]
rate_vals = [-1.0e-3, -2.0e-3, -3.0e-3]
well_head_vals = [80.0, 60.0, 40.0]
aux_vals = [12.0, 24.0, 36.0]
cont_ts_rate = [rate_vals[0]] + rate_vals
cont_ts_head = [well_head_vals[0]] + well_head_vals
cont_ts_aux = [aux_vals[0]] + aux_vals

# -- switch group --
SWITCH_NPER = 4
SWITCH_TS_TIMES = [0.0, 1.0, 2.0, 3.0, 4.0]
# period 0: literal; period 1: TS-linked; period 2: literal override;
# period 3: not reissued (must persist period 2's literal)
rate_lit0 = -1.0e-3
rate_lit2 = -5.0e-3
rate_ts_val = -7.0e-3  # distinct from rate_lit2 so link-clearing is visible
head_lit0 = 70.0
head_lit2 = 85.0
head_ts_val = 50.0

# -- partial group --
PARTIAL_NPER = 3
PARTIAL_TS_TIMES = [0.0, 1.0, 2.0, 3.0]
# well 0: literal in period 0, TS-linked from period 1 on
# well 1: literal in period 0 only, never reissued again
rate_well0_lit = -1.0e-3
rate_well1_lit = -2.0e-3
rate_ts_vals = [-6.0e-3, -9.0e-3]  # period 1, period 2 (well 0 only)
head_well0_lit = 70.0
head_well1_lit = 50.0
head_ts_vals = [55.0, 65.0]  # period 1, period 2 (well 0 only)


def _base_sim(ws, name, nper, ncol):
    """Return a minimal GWF sim ready for a MAW package to be attached."""
    period_data = [(1.0, 1, 1.0)] * nper
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=period_data)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=k11)
    flopy.mf6.ModflowGwfic(gwf, strt=chd_head)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), chd_head], [(0, 0, ncol - 1), chd_head]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )
    return sim, gwf


# ===== continue group =====


def _get_rate_continue(ws, name):
    """RATE linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name, CONT_NPER, 3)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "active"), (0, "rate", "rate_ts")],
            1: [(0, "status", "active")],
            2: [(0, "status", "active")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_rate", "RATE", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(CONT_TS_TIMES, cont_ts_rate)),
        time_series_namerecord=["rate_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_wellhead_continue(ws, name):
    """WELL_HEAD linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name, CONT_NPER, 3)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "constant"), (0, "well_head", "head_ts")],
            1: [(0, "status", "constant")],
            2: [(0, "status", "constant")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_head", "HEAD", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(CONT_TS_TIMES, cont_ts_head)),
        time_series_namerecord=["head_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_aux_continue(ws, name):
    """AUXILIARY linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name, CONT_NPER, 3)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["conc"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1, 0.0)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "active"), (0, "AUXILIARY", "conc", "conc_ts")],
            1: [(0, "status", "active")],
            2: [(0, "status", "active")],
        },
        pname="maw-1",
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(CONT_TS_TIMES, cont_ts_aux)),
        time_series_namerecord=["conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _check_continue(idx, test):
    name = cases[idx]
    if idx == 0:
        fname = test.workspace / f"{name}.{paktest}.obs.csv"
        tc = np.genfromtxt(fname, names=True, delimiter=",")
        obs = tc["WELL_RATE"]
        expected = rate_vals
        label = "RATE"
    elif idx == 1:
        fname = test.workspace / f"{name}.{paktest}.obs.csv"
        tc = np.genfromtxt(fname, names=True, delimiter=",")
        obs = tc["WELL_HEAD"]
        expected = well_head_vals
        label = "WELL_HEAD"
    else:
        fname = test.workspace / f"{name}.{paktest}.cbc"
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        auxcol = data[0].dtype.names[-1]
        obs = np.array([rec[auxcol][0] for rec in data])
        expected = aux_vals
        label = "AUXILIARY"

    assert len(obs) == CONT_NPER, (
        f"Expected {CONT_NPER} {label} records, got {len(obs)}"
    )
    assert np.allclose(obs, expected), (
        f"{label} expected to track its time series per period {expected} "
        f"(periods 1/2's PERIOD block doesn't reissue it), got {obs}."
    )


# ===== switch group =====


def _get_rate_switch(ws, name):
    sim, gwf = _base_sim(ws, name, SWITCH_NPER, 3)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "active"), (0, "rate", rate_lit0)],
            1: [(0, "status", "active"), (0, "rate", "rate_ts")],
            2: [(0, "status", "active"), (0, "rate", rate_lit2)],
            3: [(0, "status", "active")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_rate", "RATE", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(
            zip(
                SWITCH_TS_TIMES,
                [rate_lit0, rate_lit0, rate_ts_val, rate_ts_val, rate_ts_val],
            )
        ),
        time_series_namerecord=["rate_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_wellhead_switch(ws, name):
    sim, gwf = _base_sim(ws, name, SWITCH_NPER, 3)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "constant"), (0, "well_head", head_lit0)],
            1: [(0, "status", "constant"), (0, "well_head", "head_ts")],
            2: [(0, "status", "constant"), (0, "well_head", head_lit2)],
            3: [(0, "status", "constant")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_head", "HEAD", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(
            zip(
                SWITCH_TS_TIMES,
                [head_lit0, head_lit0, head_ts_val, head_ts_val, head_ts_val],
            )
        ),
        time_series_namerecord=["head_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _check_switch(local_idx, test):
    name = cases[3 + local_idx]
    fname = test.workspace / f"{name}.{paktest}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")

    if local_idx == 0:
        obs = tc["WELL_RATE"]
        expected = [rate_lit0, rate_ts_val, rate_lit2, rate_lit2]
        label = "RATE"
    else:
        obs = tc["WELL_HEAD"]
        expected = [head_lit0, head_ts_val, head_lit2, head_lit2]
        label = "WELL_HEAD"

    assert len(obs) == SWITCH_NPER, (
        f"Expected {SWITCH_NPER} {label} records, got {len(obs)}"
    )
    assert np.allclose(obs, expected), (
        f"{label} expected {expected} across periods 0-3, got {obs}."
    )


# ===== partial group =====


def _get_rate_partial(ws, name):
    """well 0's RATE switches to a TS in period 1; well 1's RATE is only
    ever set once, in period 0."""
    sim, gwf = _base_sim(ws, name, PARTIAL_NPER, 4)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=2,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[
            (0, radius, well_bot, strt, "THIEM", 1),
            (1, radius, well_bot, strt, "THIEM", 1),
        ],
        connectiondata=[
            (0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0),
            (1, 0, (0, 0, 2), top, well_bot, 0.0, 0.0),
        ],
        perioddata={
            0: [
                (0, "status", "active"),
                (1, "status", "active"),
                (0, "rate", rate_well0_lit),
                (1, "rate", rate_well1_lit),
            ],
            1: [
                (0, "status", "active"),
                (1, "status", "active"),
                (0, "rate", "rate_ts"),
            ],
            2: [
                (0, "status", "active"),
                (1, "status", "active"),
            ],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={
            f"{name}.{paktest}.obs.csv": [
                ("w0_rate", "RATE", (0,)),
                ("w1_rate", "RATE", (1,)),
            ]
        },
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(
            zip(PARTIAL_TS_TIMES, [rate_well0_lit, rate_well0_lit] + rate_ts_vals)
        ),
        time_series_namerecord=["rate_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_wellhead_partial(ws, name):
    """well 0's WELL_HEAD switches to a TS in period 1; well 1's
    WELL_HEAD is only ever set once, in period 0."""
    sim, gwf = _base_sim(ws, name, PARTIAL_NPER, 4)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=2,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[
            (0, radius, well_bot, strt, "THIEM", 1),
            (1, radius, well_bot, strt, "THIEM", 1),
        ],
        connectiondata=[
            (0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0),
            (1, 0, (0, 0, 2), top, well_bot, 0.0, 0.0),
        ],
        perioddata={
            0: [
                (0, "status", "constant"),
                (1, "status", "constant"),
                (0, "well_head", head_well0_lit),
                (1, "well_head", head_well1_lit),
            ],
            1: [
                (0, "status", "constant"),
                (1, "status", "constant"),
                (0, "well_head", "head_ts"),
            ],
            2: [
                (0, "status", "constant"),
                (1, "status", "constant"),
            ],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={
            f"{name}.{paktest}.obs.csv": [
                ("w0_head", "HEAD", (0,)),
                ("w1_head", "HEAD", (1,)),
            ]
        },
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(
            zip(PARTIAL_TS_TIMES, [head_well0_lit, head_well0_lit] + head_ts_vals)
        ),
        time_series_namerecord=["head_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _check_partial(local_idx, test):
    name = cases[5 + local_idx]
    fname = test.workspace / f"{name}.{paktest}.obs.csv"
    tc = np.genfromtxt(fname, names=True, delimiter=",")

    if local_idx == 0:
        w0_expected = [rate_well0_lit] + rate_ts_vals
        w1_expected = [rate_well1_lit] * PARTIAL_NPER
        w0_obs, w1_obs = tc["W0_RATE"], tc["W1_RATE"]
        label = "RATE"
    else:
        w0_expected = [head_well0_lit] + head_ts_vals
        w1_expected = [head_well1_lit] * PARTIAL_NPER
        w0_obs, w1_obs = tc["W0_HEAD"], tc["W1_HEAD"]
        label = "WELL_HEAD"

    assert len(w0_obs) == PARTIAL_NPER, (
        f"Expected {PARTIAL_NPER} well-0 {label} records, got {len(w0_obs)}"
    )
    assert np.allclose(w0_obs, w0_expected), (
        f"well 0 {label} expected {w0_expected}, got {w0_obs}."
    )
    assert np.allclose(w1_obs, w1_expected), (
        f"well 1 {label} expected {w1_expected} (never reissued after "
        f"period 0), got {w1_obs}."
    )


# ===== dispatch =====

_builders = [
    _get_rate_continue,
    _get_wellhead_continue,
    _get_aux_continue,
    _get_rate_switch,
    _get_wellhead_switch,
    _get_rate_partial,
    _get_wellhead_partial,
]


def build_models(idx, test):
    name = cases[idx]
    return _builders[idx](test.workspace, name)


def check_output(idx, test):
    if idx <= 2:
        _check_continue(idx, test)
    elif idx <= 4:
        _check_switch(idx - 3, test)
    else:
        _check_partial(idx - 5, test)


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
