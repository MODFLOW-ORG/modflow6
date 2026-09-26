"""LAK PACKAGEDATA AUX time-series linkage, mirroring
test_gwf_maw_pkgdata_ts.py (LAK's STRT has no time_series flag, so there
is no case for that):
  Case 0: PACKAGEDATA AUX time series matches equivalent literal PERIOD
          AUXILIARY values
  Case 1: a literal PERIOD AUXILIARY value overrides an active PACKAGEDATA
          AUX time series and persists through periods that don't reissue it
  Case 2: a PERIOD AUXILIARY time series overrides an active PACKAGEDATA
          AUX time series
  Case 3: a PACKAGEDATA AUX time series re-evaluates to a different value
          at each of several timesteps within a single stress period
  Case 4: a PERIOD AUXILIARY time series set once keeps tracking without
          being reissued (inter-period)
  Case 5: a PERIOD AUXILIARY time series re-evaluates at each of several
          timesteps within one period (intra-period)
  Case 6: a PERIOD AUXILIARY time series is overridden by a literal in a
          later period, which then holds without reissue
  Case 7: a PERIOD AUXILIARY time series with a non-matching AUXNAME is
          ignored, not applied
  Case 8: reissuing one of two AUX names with a literal does not affect
          the other AUX name's own time-series tracking
  Case 9: reissuing STAGE with a literal does not affect a separate,
          concurrently time-series-linked AUXILIARY name

Model geometry: 1-layer, 1-row, 3-column; CHD at flanking columns; single
    lake with one connection; cases 0-2, 4, 6-9 use 3 steady-state
    stress periods, 1 step each; cases 3 and 5 use 1 period with 3 steps.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.compare import eval_bud_diff
from framework import TestFramework

paktest = "lak"
cases = [
    "lak_aux_pd_ts",
    "lak_pd_oride",
    "lak_pd_ts_oride",
    "lak_substep_ts",
    "lak_auxcont",
    "lak_pd_substep",
    "lak_auxswitch",
    "lak_aux_nomatch",
    "lak_aux_dedup",
    "lak_cross_dedup",
]

# Case 3 -- one period, 3 timesteps, distinct TS value at each step end
substep_vals = [5.0, 15.0, 25.0]
substep_ts_times = [0.0, 1.0, 2.0, 3.0]
substep_ts_vals = [substep_vals[0]] + substep_vals

nper = 3
perlen = 1.0
period_data = [(perlen, 1, 1.0)] * nper

nlay, nrow, ncol = 1, 1, 3
delr = delc = 100.0
top = 200.0
botm = [0.0]
k11 = 1.0
chd_head = 100.0

strt = 100.0

# Per-period values -- distinct so any failure to update is visible
conc_vals = [10.0, 20.0, 30.0]

# PERIOD AUXILIARY override values -- distinct from conc_vals so a failure
# to apply the override is distinguishable from the PACKAGEDATA TS path.
period_conc_vals = [40.0, 50.0, 60.0]

# Case 1's literal override and its own competing PACKAGEDATA AUX TS --
# distinct from every other value so a leak between them is unambiguous.
override_conc = 20.0
pkgd_conc_vals = [100.0, 200.0, 300.0]

# TS times at period-end times so LINEAREND gives exact per-period values.
ts_times = [0.0, 1.0, 2.0, 3.0]
ts_conc = [conc_vals[0]] + conc_vals
ts_period_conc = [period_conc_vals[0]] + period_conc_vals
ts_pkgd_conc = [pkgd_conc_vals[0]] + pkgd_conc_vals

# Case 4 -- PERIOD AUXILIARY TS set once, not reissued in periods 1/2
cont_conc_vals = [70.0, 80.0, 90.0]
ts_cont_conc = [cont_conc_vals[0]] + cont_conc_vals

# Case 5 -- PERIOD AUXILIARY TS, one period, distinct value per step
pd_substep_vals = [7.0, 17.0, 27.0]
pd_substep_ts_vals = [pd_substep_vals[0]] + pd_substep_vals

# Case 6 -- PERIOD AUXILIARY TS in period 0, literal override in period 1
switch_ts_val = 11.0
switch_literal = 55.0
switch_ts_vals = [switch_ts_val] * len(ts_times)
switch_expected = [switch_ts_val, switch_literal, switch_literal]

# Case 7 -- PERIOD AUXILIARY TS with a non-matching AUXNAME
unmatched_conc = 99.0
unmatched_ts_val = 50.0

# Case 8 -- two AUX names, reissue one with a literal
dedup_temp_ts = [10.0, 10.0, 20.0, 30.0]
dedup_temp_literal = 55.0
dedup_temp_expected = [10.0, dedup_temp_literal, dedup_temp_literal]
dedup_conc_ts = [100.0, 100.0, 200.0, 300.0]
dedup_conc_expected = [100.0, 200.0, 300.0]

# Case 9 -- STAGE and AUXILIARY both TS-linked; reissue STAGE with a literal
cross_stage_ts = [5.0, 5.0, 15.0, 25.0]
cross_stage_literal = 77.0
cross_stage_expected = [5.0, cross_stage_literal, cross_stage_literal]
cross_conc_ts = [50.0, 50.0, 150.0, 250.0]
cross_conc_expected = [50.0, 150.0, 250.0]


def _base_sim(ws, name):
    """Return a minimal GWF sim ready for a LAK package to be attached."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=period_data)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
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


_connectiondata = [
    (0, 0, (0, 0, 1), "vertical", 0.0, 0.0, top, 0.0, 0.0),
]


# Case 0 -- PACKAGEDATA AUX TS vs equivalent literal PERIOD AUXILIARY
def _get_aux_period(ws, name):
    """Reference: AUX changes via PERIOD AUXILIARY each period; STRT fixed."""
    sim, gwf = _base_sim(ws, name)
    lak_spd = {
        0: [(0, "status", "constant"), (0, "AUXILIARY", "concentration", conc_vals[0])],
        1: [(0, "AUXILIARY", "concentration", conc_vals[1])],
        2: [(0, "AUXILIARY", "concentration", conc_vals[2])],
    }
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, conc_vals[0])],
        connectiondata=_connectiondata,
        perioddata=lak_spd,
        pname="lak-1",
    )
    return sim


def _get_aux_pkgdata_ts(ws, name):
    """TS: AUX changes via PACKAGEDATA AUX TS; STRT fixed (same as reference)."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, "conc_ts")],
        connectiondata=_connectiondata,
        perioddata={0: [(0, "status", "constant")]},
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_conc)),
        time_series_namerecord=["conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 1 -- literal PERIOD AUXILIARY override, issued once, over an active
# PACKAGEDATA AUX TS
def _get_pd_literal_override(ws, name):
    """A literal PERIOD AUXILIARY override in period 0 persists through periods 1/2."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, strt, 1, "pkgd_conc")],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "concentration", override_conc),
            ],
            1: [(0, "status", "constant")],
            2: [(0, "status", "constant")],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_pkgd_conc)),
        time_series_namerecord=["pkgd_conc"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 2 -- PERIOD AUX TS overrides PACKAGEDATA AUX TS
def _get_period_ts_override_ref(ws, name):
    """Reference: PERIOD AUXILIARY TS sets the value; no PACKAGEDATA AUX TS."""
    sim, gwf = _base_sim(ws, name)
    lak_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", "period_conc_ts"),
        ],
        1: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
        2: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
    }
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, period_conc_vals[0])],
        connectiondata=_connectiondata,
        perioddata=lak_spd,
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_period_conc)),
        time_series_namerecord=["period_conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_period_ts_override_ts(ws, name):
    """TS: PACKAGEDATA AUX TS and PERIOD AUXILIARY TS both active; PERIOD TS wins."""
    sim, gwf = _base_sim(ws, name)
    lak_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", "period_conc_ts"),
        ],
        1: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
        2: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
    }
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, "conc_ts")],
        connectiondata=_connectiondata,
        perioddata=lak_spd,
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=[
            (ts_times[i], ts_conc[i], ts_period_conc[i]) for i in range(len(ts_times))
        ],
        time_series_namerecord=["conc_ts", "period_conc_ts"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    return sim


def _get_aux_substep_ts(ws, name):
    """One period, 3 steps; PACKAGEDATA AUX TS must re-evaluate every step."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(3.0, 3, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
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
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, "substep_ts")],
        connectiondata=_connectiondata,
        perioddata={0: [(0, "status", "constant")]},
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(substep_ts_times, substep_ts_vals)),
        time_series_namerecord=["substep_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 4 -- PERIOD AUXILIARY TS set once, continues without reissue
def _get_aux_continue(ws, name):
    """AUXILIARY TS set once; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, cont_conc_vals[0])],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "concentration", "cont_conc_ts"),
            ],
            1: [(0, "status", "constant")],
            2: [(0, "status", "constant")],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_cont_conc)),
        time_series_namerecord=["cont_conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 5 -- PERIOD AUXILIARY TS, one period, 3 steps
def _get_aux_period_substep(ws, name):
    """PERIOD AUXILIARY TS must re-evaluate every step, not just period start."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(3.0, 3, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
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
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, pd_substep_vals[0])],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "concentration", "pd_substep_ts"),
            ],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(substep_ts_times, pd_substep_ts_vals)),
        time_series_namerecord=["pd_substep_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 6 -- PERIOD AUXILIARY TS in period 0, literal override in period 1
def _get_aux_switch(ws, name):
    """AUXILIARY: TS in period 0, literal override in period 1, held in period 2."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, switch_ts_val)],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "concentration", "switch_ts"),
            ],
            1: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "concentration", switch_literal),
            ],
            2: [(0, "status", "constant")],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, switch_ts_vals)),
        time_series_namerecord=["switch_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 7 -- PERIOD AUXILIARY TS with a non-matching AUXNAME
def _get_aux_unmatched_name(ws, name):
    """A non-matching AUXNAME must be ignored, not applied."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, unmatched_conc)],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "wrongname", "unmatched_ts"),
            ],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, [unmatched_ts_val] * len(ts_times))),
        time_series_namerecord=["unmatched_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 8 -- two AUX names TS-linked; reissue one with a literal
def _get_aux_dedup(ws, name):
    """TEMP reissued as a literal in period 1; CONC keeps tracking its own TS."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["temp", "conc"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, dedup_temp_ts[0], dedup_conc_ts[0])],
        connectiondata=_connectiondata,
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "temp", "dedup_temp_ts"),
                (0, "AUXILIARY", "conc", "dedup_conc_ts"),
            ],
            1: [
                (0, "status", "constant"),
                (0, "AUXILIARY", "temp", dedup_temp_literal),
            ],
            2: [(0, "status", "constant")],
        },
        pname="lak-1",
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=[
            (t, temp, conc)
            for t, temp, conc in zip(ts_times, dedup_temp_ts, dedup_conc_ts)
        ],
        time_series_namerecord=["dedup_temp_ts", "dedup_conc_ts"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    return sim


# Case 9 -- STAGE and AUXILIARY both TS-linked; reissue STAGE with a literal
def _get_cross_dedup(ws, name):
    """STAGE reissued as a literal in period 1; AUXILIARY keeps tracking its own TS."""
    sim, gwf = _base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, strt, 1, cross_conc_ts[0])],
        connectiondata=[
            (0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, top, 10.0, 10.0),
        ],
        perioddata={
            0: [
                (0, "status", "constant"),
                (0, "stage", "cross_stage_ts"),
                (0, "AUXILIARY", "concentration", "cross_conc_ts"),
            ],
            1: [
                (0, "status", "constant"),
                (0, "stage", cross_stage_literal),
            ],
            2: [(0, "status", "constant")],
        },
        pname="lak-1",
    )
    lak.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("stg1", "stage", (0,))]},
    )
    lak.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=[
            (t, stage, conc)
            for t, stage, conc in zip(ts_times, cross_stage_ts, cross_conc_ts)
        ],
        time_series_namerecord=["cross_stage_ts", "cross_conc_ts"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")
    if idx == 0:
        # Case 0: baseline uses PERIOD AUXILIARY; TS uses PACKAGEDATA AUX TS.
        return _get_aux_period(ws0, name), _get_aux_pkgdata_ts(ws1, name)
    elif idx == 1:
        # Case 1: literal override, not compared against a reference; the
        # AUXILIARY budget term (in mf6/) is checked directly.
        return _get_pd_literal_override(ws0, name), _get_pd_literal_override(ws1, name)
    elif idx == 2:
        # Case 2: PERIOD AUX TS must win over PACKAGEDATA AUX TS.
        return (
            _get_period_ts_override_ref(ws0, name),
            _get_period_ts_override_ts(ws1, name),
        )
    elif idx == 3:
        # Case 3: not compared against a reference; checked directly.
        return _get_aux_substep_ts(ws0, name), _get_aux_substep_ts(ws1, name)
    elif idx == 4:
        # Case 4: not compared against a reference; checked directly.
        return _get_aux_continue(ws0, name), _get_aux_continue(ws1, name)
    elif idx == 5:
        # Case 5: not compared against a reference; checked directly.
        return _get_aux_period_substep(ws0, name), _get_aux_period_substep(ws1, name)
    elif idx == 6:
        # Case 6: not compared against a reference; checked directly.
        return _get_aux_switch(ws0, name), _get_aux_switch(ws1, name)
    elif idx == 7:
        # Case 7: not compared against a reference; checked directly.
        return _get_aux_unmatched_name(ws0, name), _get_aux_unmatched_name(ws1, name)
    elif idx == 8:
        # Case 8: not compared against a reference; checked directly.
        return _get_aux_dedup(ws0, name), _get_aux_dedup(ws1, name)
    else:
        # Case 9: not compared against a reference; checked directly.
        return _get_cross_dedup(ws0, name), _get_cross_dedup(ws1, name)


def check_output(idx, test):
    name = os.path.basename(test.name)
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")

    if idx == 1:
        # -- Case 1: the literal override holds through periods 1/2, on
        #    both connections --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        expected = [override_conc] * nper
        assert np.allclose(aux_obs, expected), (
            f"AUXILIARY expected to hold the literal override {expected} "
            f"through periods 1/2, which don't reissue it, got {aux_obs}."
        )
    elif idx == 3:
        # -- Case 3: each of the 3 steps in the one period must show a
        #    different, TS-correct value, not the period's initial value --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == 3, f"Expected 3 AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        assert np.allclose(aux_obs, substep_vals), (
            f"AUXILIARY expected {substep_vals} at each of the 3 steps, got {aux_obs}."
        )
    elif idx == 4:
        # -- Case 4: AUXILIARY TS must keep tracking through periods 1/2,
        #    which don't reissue it --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        assert np.allclose(aux_obs, cont_conc_vals), (
            f"AUXILIARY expected to keep tracking its TS {cont_conc_vals} "
            f"through periods 1/2, got {aux_obs}."
        )
    elif idx == 5:
        # -- Case 5: AUXILIARY TS must re-evaluate at each of the 3 steps --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == 3, f"Expected 3 AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        assert np.allclose(aux_obs, pd_substep_vals), (
            f"AUXILIARY expected {pd_substep_vals} at each of the 3 steps, "
            f"got {aux_obs}."
        )
    elif idx == 6:
        # -- Case 6: literal override in period 1 must win and hold into
        #    period 2, not revert to the period-0 TS --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        assert np.allclose(aux_obs, switch_expected), (
            f"AUXILIARY expected {switch_expected} (TS, then literal, then "
            f"held), got {aux_obs}."
        )
    elif idx == 7:
        # -- Case 7: non-matching AUXNAME must be ignored --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        aux_obs = np.array([rec[auxcol][0] for rec in data])
        expected = [unmatched_conc] * nper
        assert np.allclose(aux_obs, expected), (
            f"AUXILIARY expected to stay at {expected} (unmatched AUXNAME "
            f"ignored), got {aux_obs}."
        )
    elif idx == 8:
        # -- Case 8: TEMP reissued as a literal must not disturb CONC's
        #    own time series --
        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        temp_obs = np.array([rec["TEMP"][0] for rec in data])
        conc_obs = np.array([rec["CONC"][0] for rec in data])
        assert np.allclose(temp_obs, dedup_temp_expected), (
            f"TEMP expected {dedup_temp_expected}, got {temp_obs}."
        )
        assert np.allclose(conc_obs, dedup_conc_expected), (
            f"CONC expected {dedup_conc_expected} (own time series, "
            f"unaffected by TEMP's reissue), got {conc_obs}."
        )
    elif idx == 9:
        # -- Case 9: STAGE reissued as a literal must not disturb
        #    CONCENTRATION's own time series --
        obs_path = os.path.join(ws1, f"{name}.{paktest}.obs.csv")
        obs = np.genfromtxt(obs_path, delimiter=",", names=True)
        stage_obs = obs["STG1"]
        assert np.allclose(stage_obs, cross_stage_expected), (
            f"STAGE expected {cross_stage_expected}, got {stage_obs}."
        )

        fname = os.path.join(ws1, f"{name}.{paktest}.cbc")
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        assert len(data) == nper, f"Expected {nper} AUXILIARY records, got {len(data)}"
        auxcol = data[0].dtype.names[-1]
        conc_obs = np.array([rec[auxcol][0] for rec in data])
        assert np.allclose(conc_obs, cross_conc_expected), (
            f"CONCENTRATION expected {cross_conc_expected} (own time "
            f"series, unaffected by STAGE's reissue), got {conc_obs}."
        )
    else:
        # Cases 0 and 2: GWF and LAK budgets must be identical between
        # reference and TS model. For case 2 this asserts that PERIOD TS
        # overrides PACKAGEDATA AUX TS.
        ia = (
            flopy.mf6.utils.MfGrdFile(os.path.join(ws0, f"{name}.dis.grb"))._datadict[
                "IA"
            ]
            - 1
        )

        eval_bud_diff(
            os.path.join(ws0, f"{name}.cbc.cmp.out"),
            flopy.utils.CellBudgetFile(
                os.path.join(ws0, f"{name}.cbc"), precision="double"
            ),
            flopy.utils.CellBudgetFile(
                os.path.join(ws1, f"{name}.cbc"), precision="double"
            ),
            ia,
        )
        eval_bud_diff(
            os.path.join(ws0, f"{name}.{paktest}.cbc.cmp.out"),
            flopy.utils.CellBudgetFile(
                os.path.join(ws0, f"{name}.{paktest}.cbc"), precision="double"
            ),
            flopy.utils.CellBudgetFile(
                os.path.join(ws1, f"{name}.{paktest}.cbc"), precision="double"
            ),
        )


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
