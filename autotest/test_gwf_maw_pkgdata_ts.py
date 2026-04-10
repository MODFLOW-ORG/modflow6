"""test_gwf_maw_pkgdata_ts.py

Focused tests for MAW PACKAGEDATA IDM enabled time-series linkage for two
distinct TS mechanisms:
  Case 0: STRT package data timeseries variable
  Case 1: AUX package data timeseries variable
  Case 2: verify PERIOD AUX data array override with PACKAGEDATA AUX TS
  Case 3: verify PERIOD AUX TS data array override with PACKAGEDATA AUX TS

Model geometry: 1-layer, 1-row, 3-column; CHD at flanking columns; single
    MAW well in the centre column; 3 steady-state stress periods, 1 step each.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.compare import eval_bud_diff
from framework import TestFramework

paktest = "maw"
cases = ["maw_strt_pd_ts", "maw_aux_pd_ts", "maw_pd_oride", "maw_pd_ts_oride"]

# Shared geometry and temporal discretisation
nper = 3
perlen = 1.0
period_data = [(perlen, 1, 1.0)] * nper  # steady-state, 1 timestep each

nlay, nrow, ncol = 1, 1, 3
delr = delc = 100.0
top = 200.0
botm = [0.0]
k11 = 1.0
chd_head = 100.0

# MAW geometry
radius = 0.1
well_bot = -10.0

# Per-period values — chosen to be distinct so any failure to update is visible
strt_vals = [90.0, 95.0, 85.0]
conc_vals = [10.0, 20.0, 30.0]

# PERIOD AUXILIARY override values — different from conc_vals so that a failure
# to apply the PERIOD override is clearly distinguishable from the PKGDATA TS path.
period_conc_vals = [40.0, 50.0, 60.0]

# TS times: entries at period-end times so LINEAREND gives exact per-period values.
# tsmanager%ad() is called with totim = end of current timestep (1.0, 2.0, 3.0).
ts_times = [0.0, 1.0, 2.0, 3.0]
ts_strt = [strt_vals[0]] + strt_vals  # [90, 90, 95, 85]
ts_conc = [conc_vals[0]] + conc_vals  # [10, 10, 20, 30]
ts_period_conc = [period_conc_vals[0]] + period_conc_vals  # [40, 40, 50, 60]


def _base_sim(ws, name):
    """Return a minimal GWF sim ready for a MAW package to be attached."""
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
        stress_period_data=[[(0, 0, 0), chd_head], [(0, 0, 2), chd_head]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )
    return sim, gwf


# Case 0 — STRT PACKAGEDATA TS
def _get_strt_ts(ws, name):
    """STRT as PACKAGEDATA TS; STATUS CONSTANT; HEAD obs recorded."""
    sim, gwf = _base_sim(ws, name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, "strt_ts", "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={0: [(0, "status", "constant")]},
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        digits=10,
        print_input=True,
        continuous={"maw_head.csv": [("well_head", "HEAD", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_strt)),
        time_series_namerecord=["strt_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 1 — PACKAGEDATA AUX TS
def _get_aux_period(ws, name):
    """Reference: AUX changes via PERIOD AUXILIARY each period; STRT fixed."""
    sim, gwf = _base_sim(ws, name)
    maw_spd = {
        0: [(0, "status", "constant"), (0, "AUXILIARY", "concentration", conc_vals[0])],
        1: [(0, "AUXILIARY", "concentration", conc_vals[1])],
        2: [(0, "AUXILIARY", "concentration", conc_vals[2])],
    }
    flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, radius, well_bot, strt_vals[0], "THIEM", 1, conc_vals[0])],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata=maw_spd,
        pname="maw-1",
    )
    return sim


def _get_aux_pkgdata_ts(ws, name):
    """TS: AUX changes via PACKAGEDATA AUX TS; STRT fixed (same as reference)."""
    sim, gwf = _base_sim(ws, name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, radius, well_bot, strt_vals[0], "THIEM", 1, "conc_ts")],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={0: [(0, "status", "constant")]},
        pname="maw-1",
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_conc)),
        time_series_namerecord=["conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 2 — PERIOD override PACKAGEDATA AUX TS
# Reference model: PERIOD AUXILIARY with literal period_conc_vals; no PKGDATA TS.
# TS model:        PACKAGEDATA AUX TS (conc_vals) + same PERIOD AUXILIARY override.
# If PERIOD correctly overrides, both models produce identical budget files.
def _get_pd_ovrride_ref(ws, name):
    """Reference: PERIOD AUXILIARY literal values; no PACKAGEDATA AUX TS."""
    sim, gwf = _base_sim(ws, name)
    maw_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", period_conc_vals[0]),
        ],
        1: [(0, "AUXILIARY", "concentration", period_conc_vals[1])],
        2: [(0, "AUXILIARY", "concentration", period_conc_vals[2])],
    }
    flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[
            (0, radius, well_bot, strt_vals[0], "THIEM", 1, period_conc_vals[0])
        ],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata=maw_spd,
        pname="maw-1",
    )
    return sim


def _get_pd_ovrride_ts(ws, name):
    """TS: PKGDATA AUX TS (conc_vals) + PERIOD AUXILIARY override (period_conc_vals).

    PERIOD must win: auxvar should equal period_conc_vals each period, not conc_vals.
    """
    sim, gwf = _base_sim(ws, name)
    maw_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", period_conc_vals[0]),
        ],
        1: [(0, "AUXILIARY", "concentration", period_conc_vals[1])],
        2: [(0, "AUXILIARY", "concentration", period_conc_vals[2])],
    }
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, radius, well_bot, strt_vals[0], "THIEM", 1, "conc_ts")],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata=maw_spd,
        pname="maw-1",
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_conc)),
        time_series_namerecord=["conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


# Case 3 — PERIOD AUX TS override PACKAGEDATA AUX TS
# Reference model: PERIOD AUXILIARY TS (period_conc_ts) drives auxvar; no PKGDATA TS.
# TS model:        PACKAGEDATA AUX TS (conc_ts, conc_vals) +
#                  PERIOD AUXILIARY TS (period_conc_ts, period_conc_vals).
# If PERIOD TS correctly overrides, both models produce identical budget files.
def _get_period_ts_ovrride_ref(ws, name):
    """Reference: PERIOD AUXILIARY TS drives auxvar; no PACKAGEDATA AUX TS."""
    sim, gwf = _base_sim(ws, name)
    maw_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", "period_conc_ts"),
        ],
        1: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
        2: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
    }
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[
            (0, radius, well_bot, strt_vals[0], "THIEM", 1, period_conc_vals[0])
        ],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata=maw_spd,
        pname="maw-1",
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_period_conc)),
        time_series_namerecord=["period_conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_period_ts_ovrride_ts(ws, name):
    """TS: PKGDATA AUX TS (conc_ts) + PERIOD AUXILIARY TS (period_conc_ts).

    PERIOD TS must win: auxvar should equal period_conc_vals, not conc_vals.
    Both TS series live in the same MAW TS file.
    """
    sim, gwf = _base_sim(ws, name)
    maw_spd = {
        0: [
            (0, "status", "constant"),
            (0, "AUXILIARY", "concentration", "period_conc_ts"),
        ],
        1: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
        2: [(0, "AUXILIARY", "concentration", "period_conc_ts")],
    }
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["concentration"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        packagedata=[(0, radius, well_bot, strt_vals[0], "THIEM", 1, "conc_ts")],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata=maw_spd,
        pname="maw-1",
    )
    # Both TS series in a single file: (time, conc_ts_val, period_conc_ts_val)
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=[
            (ts_times[i], ts_conc[i], ts_period_conc[i]) for i in range(len(ts_times))
        ],
        time_series_namerecord=["conc_ts", "period_conc_ts"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    return sim


# Build the model
def build_models(idx, test):
    name = cases[idx]
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")
    if idx == 0:
        # Case 0: baseline is a fixed-STRT model (not compared); TS model
        # is in mf6/ sub-directory where HEAD obs are checked.
        return _get_strt_ts(ws0, name), _get_strt_ts(ws1, name)
    elif idx == 1:
        # Case 1: baseline uses PERIOD AUXILIARY; TS uses PACKAGEDATA AUX TS.
        return _get_aux_period(ws0, name), _get_aux_pkgdata_ts(ws1, name)
    elif idx == 2:
        # Case 2: PERIOD override (literal) must win over PACKAGEDATA AUX TS.
        return _get_pd_ovrride_ref(ws0, name), _get_pd_ovrride_ts(ws1, name)
    else:
        # Case 3: PERIOD AUX TS must win over PACKAGEDATA AUX TS.
        return (
            _get_period_ts_ovrride_ref(ws0, name),
            _get_period_ts_ovrride_ts(ws1, name),
        )


# Check output
def check_output(idx, test):
    name = os.path.basename(test.name)
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")

    if idx == 0:
        # -- Case 0: assert that MAW HEAD obs (in mf6/) match STRT TS values --
        obs_path = os.path.join(ws1, "maw_head.csv")
        obs = np.genfromtxt(obs_path, delimiter=",", names=True, skip_header=0)
        # One output time per stress period; values should equal strt_vals exactly
        # (LINEAREND at period-end totim gives the exact TS entry value)
        heads = obs["WELL_HEAD"]
        assert len(heads) == nper, f"expected {nper} obs records, got {len(heads)}"
        assert np.allclose(heads, strt_vals, atol=1e-6), (
            f"STRT TS: well_head {heads} does not match expected {strt_vals}"
        )
    else:
        # Cases 1-3: GWF and MAW budgets must be identical between reference and
        # TS model.  For cases 2 and 3 this asserts that PERIOD (or PERIOD TS)
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
