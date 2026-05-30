"""
Tests AUTO_FLOW_REDUCE (AFR) behavior in the WEL package.

Scenario: single isolated unconfined cell, 10 m tall, drained by a well with
AFR enabled.  Because cell thickness (10 m) differs from the AFR unit length
(1 m), fraction mode and length mode produce clearly distinct turnoff thresholds:

  Fraction mode: tp = botm + flowred * (top - botm) = 0 + 0.5 * 10 = 5.0 m
  Length  mode:  tp = botm + flowred * 1.0           = 0 + 0.5       = 0.5 m

Each case runs Newton (main workspace) and non-Newton (mf6/ subdirectory).
The TestFramework's automatic head comparison verifies Newton ≈ non-Newton.
check_output() additionally verifies that the two AFR modes produce heads on
opposite sides of a 1.0 m threshold after 4 days of pumping.

Analytical estimates (smoothstep saturation, Euler forward, dt=0.1 d):
  fraction: full pumping for ~2 d (h: 10→5), then reduced → h(4d) ≈ 1.7 m
  length:   full pumping for ~3.8 d (h: 10→0.5), then near-zero → h(4d) ≈ 0.1 m
The 1.0 m threshold separates them with > 0.7 m margin on each side.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Two AFR-mode cases; each case compares Newton vs non-Newton.
cases = ["wel02-frac", "wel02-len"]

# --- Grid: single isolated cell, 10 m thick ---
nlay, nrow, ncol = 1, 1, 1
top = 10.0
botm = [0.0]
delr = delc = 1.0

# --- Hydraulic properties ---
hk = 1.0
sy = 0.1
strt = 10.0

# --- WEL parameters ---
wellq = -0.25  # m³/d, pumping
flowred = 0.5
# Turnoff thresholds (gwf-wel.f90 formulas):
_tp_frac = botm[0] + flowred * (top - botm[0])  # 5.0 m  (fraction mode)
_tp_len = botm[0] + flowred * 1.0  # 0.5 m  (length  mode)
# Head threshold that fraction mode exceeds and length mode falls below.
_h_threshold = 1.0  # m

# --- Time: 4 days, 40 equal steps of 0.1 d ---
nper = 1
tdis_rc = [(4.0, 40, 1.0)]

# --- Solver ---
nouter, ninner = 500, 300
hclose = 1e-9
rcloserecord = 1e-3  # m3/d; framework uses 5x this for comparison
relax = 1.0


def _build_sim(ws, name, *, use_length, use_newton):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rcloserecord,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="NEWTON UNDER_RELAXATION" if use_newton else None,
        save_flows=True,
    )
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk, k33=hk)
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=sy,
        transient={0: True},
    )
    flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        auto_flow_reduce=flowred,
        flow_reduction_length=True if use_length else None,
        stress_period_data={0: [[0, 0, 0, wellq]]},
        afrcsv_filerecord=f"{name}.afr.csv",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    use_length = name.endswith("-len")
    sim = _build_sim(test.workspace, name, use_length=use_length, use_newton=True)
    mc = _build_sim(
        test.workspace / "mf6", name, use_length=use_length, use_newton=False
    )
    return sim, mc


def check_output(idx, test):
    name = cases[idx]

    # Read all-time heads from Newton (main) and non-Newton (mf6/) runs.
    h_nwt = flopy.utils.HeadFile(test.workspace / f"{name}.hds").get_alldata()
    h_std = flopy.utils.HeadFile(test.workspace / "mf6" / f"{name}.hds").get_alldata()

    # Newton and non-Newton must agree at every time step.
    # 1 mm tolerance is orders of magnitude above cross-compiler floating-point
    # differences for a solution converged to hclose=1e-9.
    max_diff = float(np.max(np.abs(h_nwt - h_std)))
    assert max_diff < 1e-3, (
        f"{name}: max Newton vs non-Newton head difference = {max_diff:.2e} m "
        f"(must be < 1e-3 m)"
    )

    # After 4 days, fraction mode (tp={_tp_frac} m) retains significantly more
    # water than length mode (tp={_tp_len} m) because pumping is curtailed much
    # earlier.  Both sides of the 1.0 m threshold have > 0.7 m of margin.
    h_final = float(h_nwt[-1, 0, 0, 0])
    if name.endswith("-frac"):
        assert h_final > _h_threshold, (
            f"{name}: fraction-mode final head {h_final:.4f} m must be "
            f"> {_h_threshold} m (tp_frac={_tp_frac} m reduces pumping early)"
        )
    else:
        assert h_final < _h_threshold, (
            f"{name}: length-mode final head {h_final:.4f} m must be "
            f"< {_h_threshold} m (tp_len={_tp_len} m reduces pumping late)"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare="mf6",
    )
    test.run()
