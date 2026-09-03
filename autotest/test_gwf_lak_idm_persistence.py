"""LAK's cross-period value persistence and precedence.

- OUTLETS geometry (INVERT, WIDTH, ROUGH, SLOPE) must be time-series
  capable and re-synced each period, with no PERIOD override involved.
  Reference model bakes in the period-2 geometry from the start; TS model
  uses different geometry in period 1 and switches to the same period-2
  geometry via a time series. Their period-2 outlet flow must match.
- A PERIOD-block override of that same outlet geometry must persist into
  later stress periods that don't repeat it, taking precedence over the
  OUTLETS block's own separately time-series-linked baseline geometry for
  that same outlet -- per the PERIOD block's documented persistence rule
  ("will continue to apply for subsequent stress periods ... until another
  PERIOD block is encountered").
- PERIOD STAGE must override the PACKAGEDATA STRT baseline for a CONSTANT
  lake.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-outlets-ts", "lak-outlet-period-precedence", "lak-strt-stage"]

# -- outlets-ts case data
ts_stage = 10.0
invert2, width2, rough2, slope2 = 2.0, 5.0, 0.03, 0.001  # period-2 geometry
invert1, width1, rough1, slope1 = 8.0, 1.0, 0.1, 0.02  # period-1 (TS model only)

outlets_ts_times = [0.0, 1.0, 2.0]
invert_ts_vals = [invert1, invert1, invert2]
width_ts_vals = [width1, width1, width2]
rough_ts_vals = [rough1, rough1, rough2]
slope_ts_vals = [slope1, slope1, slope2]

# -- period-precedence case data
prec_stage = 100.0
prec_width = 5.0
invert_a, invert_b, invert_c, invert_d = 90.0, 95.0, 92.0, 98.0

# governs periods 0, 1, 2 respectively (OUTLETS block's own baseline,
# absent any PERIOD override); invert_b is never expected to take effect,
# since period 1's PERIOD-block override should preempt it.
prec_ts_times = [0.0, 1.0, 2.0, 3.0]
prec_invert_ts = [invert_a, invert_a, invert_b, invert_d]

# -- strt-stage case data
strt = 100.0
stage_override = 80.0


def _outlets_ts_base_sim(ws, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 2), 100.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim, gwf


def _outlets_ts_ref(ws, name):
    sim, gwf = _outlets_ts_base_sim(ws, name)
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        noutlets=1,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, ts_stage, 1)],
        connectiondata=[(0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)],
        outlets=[(0, 0, -1, "manning", invert2, width2, rough2, slope2)],
        perioddata={0: [(0, "status", "constant"), (0, "stage", ts_stage)]},
        pname="lak-1",
        observations={"lak_outlet.csv": [("out1", "outlet", 1)]},
    )
    return sim


def _outlets_ts_ts(ws, name):
    sim, gwf = _outlets_ts_base_sim(ws, name)
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        noutlets=1,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, ts_stage, 1)],
        connectiondata=[(0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)],
        outlets=[
            (0, 0, -1, "manning", "invert_ts", "width_ts", "rough_ts", "slope_ts")
        ],
        perioddata={0: [(0, "status", "constant"), (0, "stage", ts_stage)]},
        pname="lak-1",
        observations={"lak_outlet.csv": [("out1", "outlet", 1)]},
    )
    lak.ts.initialize(
        filename=f"{name}.lak.ts",
        timeseries=list(
            zip(
                outlets_ts_times,
                invert_ts_vals,
                width_ts_vals,
                rough_ts_vals,
                slope_ts_vals,
            )
        ),
        time_series_namerecord=["invert_ts", "width_ts", "rough_ts", "slope_ts"],
        interpolation_methodrecord=["linearend"] * 4,
    )
    return sim


def build_outlets_ts(test):
    name = cases[0]
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")
    return _outlets_ts_ref(ws0, name), _outlets_ts_ts(ws1, name)


def build_period_precedence(test):
    name = cases[1]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(1.0, 1, 1.0)] * 3)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="lak-prec", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 2), 100.0]]
    )
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        noutlets=1,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, prec_stage, 1)],
        connectiondata=[(0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)],
        outlets=[(0, 0, -1, "weir", "invert_ts", prec_width, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "constant"), (0, "stage", prec_stage)],
            # period 1 override only -- period 2 deliberately has no PERIOD
            # block at all, so it must inherit this override, not the
            # OUTLETS block's own invert_ts value for period 2 (invert_d).
            1: [(0, "invert", invert_c)],
        },
        pname="lak-1",
        observations={"lak_outlet.csv": [("out1", "outlet", 1)]},
    )
    lak.ts.initialize(
        filename=f"{name}.lak.ts",
        timeseries=list(zip(prec_ts_times, prec_invert_ts)),
        time_series_namerecord=["invert_ts"],
        interpolation_methodrecord=["linearend"],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def build_strt_stage(test):
    name = cases[2]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=200.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 0, 2), 100.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        budget_filerecord=f"{name}.lak.cbc",
        packagedata=[(0, strt, 1)],
        connectiondata=[(0, 0, (0, 0, 1), "horizontal", 0.0, 0.0, 150.0, 10.0, 10.0)],
        perioddata={
            0: [(0, "status", "inactive")],
            1: [(0, "status", "constant"), (0, "stage", stage_override)],
        },
        pname="lak-1",
        observations={"lak_stage.csv": [("stg1", "stage", 1)]},
    )
    return sim


def build_models(idx, test):
    return [build_outlets_ts, build_period_precedence, build_strt_stage][idx](test)


def check_outlets_ts(test):
    ws0 = test.workspace
    ws1 = os.path.join(test.workspace, "mf6")

    obs_ref = np.genfromtxt(
        os.path.join(ws0, "lak_outlet.csv"), delimiter=",", names=True
    )
    obs_ts = np.genfromtxt(
        os.path.join(ws1, "lak_outlet.csv"), delimiter=",", names=True
    )
    flow_ref = float(obs_ref["OUT1"][1])
    flow_ts = float(obs_ts["OUT1"][1])
    assert np.isclose(flow_ref, flow_ts, rtol=1e-6), (
        f"period-2 outlet flow should match once TS geometry resyncs to the "
        f"reference values; ref={flow_ref}, ts={flow_ts}"
    )


def check_period_precedence(test):
    obs = np.genfromtxt(test.workspace / "lak_outlet.csv", delimiter=",", names=True)
    flow0, flow1, flow2 = (float(obs["OUT1"][i]) for i in range(3))
    assert not np.isclose(flow0, flow1), (
        "period 2's INVERT override should produce a different outlet flow "
        f"than period 1's OUTLETS-block baseline; got flow0={flow0}, "
        f"flow1={flow1}"
    )
    assert np.isclose(flow1, flow2, rtol=1e-6), (
        "period 2's INVERT override should persist into period 3 (no new "
        "PERIOD block there), not be overwritten by the OUTLETS block's "
        f"own resync to its next TS value; got flow1={flow1}, flow2={flow2}"
    )


def check_strt_stage(test):
    obs_path = os.path.join(test.workspace, "lak_stage.csv")
    obs = np.genfromtxt(obs_path, delimiter=",", names=True, skip_header=0)
    stage = float(obs["STG1"][1])
    assert stage == stage_override, (
        f"CONSTANT lake stage should equal the PERIOD STAGE value "
        f"({stage_override}), not the PACKAGEDATA STRT baseline ({strt}); "
        f"got {stage}"
    )


def check_output(idx, test):
    [check_outlets_ts, check_period_precedence, check_strt_stage][idx](test)


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
