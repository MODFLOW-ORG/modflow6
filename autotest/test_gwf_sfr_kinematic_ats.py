"""
Tests for the SFR kinematic-wave ATS_COURANT option.

Single-reach Ponce (1989) Ex. 9-3 routing problem.

Cases:
  - sfr-kwats01 : ATS_COURANT with ATS active; dt0 overshoots so the Courant
                  constraint drives sub-steps down. Verifies multiple sub-steps
                  and that the Courant table is written.
  - sfr-kwats02 : ATS_COURANT without ATS in TDIS; verifies the warning.

Ponce, V. M. (1989). Engineering Hydrology, Principles and Practices.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ("sfr-kwats01", "sfr-kwats02")

# Ponce (1989) Ex 9-3 hydrograph (non-zero flows only; 10 intervals)
_dt = 3600.0  # s, original hydrograph time step
_flows = np.array([30.0, 60.0, 90.0, 120.0, 150.0, 120.0, 90.0, 60.0, 30.0, 10.0])
_nsteps = len(_flows)  # 10 hydrograph intervals

# SFR channel parameters (shared between depth helper and build_models)
_B = 10.0
_dx = 7200.0
_slope = 1.0 / _dx
_roughness = 0.03574737676661647


def _q_from_depth(d):
    """Manning's equation for a rectangular channel."""
    if d <= 0.0:
        return 0.0
    A = _B * d
    P = _B + 2.0 * d
    return (1.0 / _roughness) * A * (A / P) ** (2.0 / 3.0) * _slope**0.5


def _depth_from_q(Q):
    """Invert Manning's equation by bisection."""
    if Q <= 0.0:
        return 0.0
    hi = 1.0
    while _q_from_depth(hi) < Q:
        hi *= 2.0
    lo = 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _q_from_depth(mid) < Q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Shared GWF geometry builder
# ---------------------------------------------------------------------------


def _build_gwf(sim, name, ncol=6):
    """Build a minimal GWF model inside *sim* and return it."""
    nlay, nrow = 1, 1
    top = 20.0
    botm = 0.0

    flopy.mf6.ModflowIms(sim, outer_dvclose=1.0e-8, inner_dvclose=1.0e-9)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=_dx,
        delc=_dx,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0e-4)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            (0, 0, 0, 11.0),
            (0, 0, ncol - 1, 9.0),
        ],
    )
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return gwf


# ---------------------------------------------------------------------------
# Case 0: ATS_COURANT + ATS active
# ---------------------------------------------------------------------------


def _build_ats_active(test):
    name = cases[0]
    perlen = _nsteps * _dt  # 36000 s

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
        version="mf6",
        exe_name="mf6",
    )
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=1,
        perioddata=[(perlen, 1, 1.0)],
    )

    # ATS: dt0 = 4x the Courant-optimal step so the first step overshoots and
    # the Courant constraint actively reduces subsequent steps to ~ _dt.
    dt0 = 4.0 * _dt  # 14400 s
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, dt0, 30.0, dt0, 2.0, 5.0)],
        filename=f"{name}.ats",
    )

    gwf = _build_gwf(sim, name)

    rtp = 20.0
    d0 = _depth_from_q(_flows[0])
    stage0 = rtp + d0
    pak_data = [(0, -1, -1, -1, _dx, _B, _slope, rtp, 1.0, 0.0, _roughness, 0, 0.0, 0)]

    # Inflow varies via a time series with stepwise interpolation.
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_flows=True,
        storage=True,
        ats_courant=1.0,
        maximum_depth_change=1.0e-9,
        nreaches=1,
        packagedata=pak_data,
        connectiondata=[(0,)],
        perioddata={0: [(0, "inflow", "inflow-ts")]},
        initialstages=[(0, stage0)],
        pname="sfr-1",
    )
    ts_times = [i * _dt for i in range(_nsteps)] + [perlen]
    ts_vals = _flows.tolist() + [0.0]
    sfr.ts.initialize(
        filename=f"{name}.sfr.ts",
        timeseries=list(zip(ts_times, ts_vals)),
        time_series_namerecord=[("inflow-ts",)],
        interpolation_methodrecord=[("stepwise",)],
    )

    obs_fname = f"{name}.sfr.obs"
    sfr.obs.initialize(
        filename=obs_fname,
        print_input=True,
        continuous={
            f"{obs_fname}.csv": [
                ("inflow", "ext-inflow", (0,)),
                ("outflow", "ext-outflow", (0,)),
            ]
        },
    )

    return sim


def _check_ats_active(test):
    name = cases[0]

    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    assert len(obs) > 0, "no observations written"

    # With dt0 = 4 * _dt and ATS_COURANT = 1, the first step takes 14400 s
    # and subsequent steps are driven down to ~ _dt = 3600 s.  The total
    # number of steps in the 36000 s period should be well above 1.
    assert len(obs) >= 5, (
        f"expected >= 5 time steps from ATS sub-stepping, got {len(obs)}"
    )

    totim = obs["totim"]
    dt_obs = np.diff(totim, prepend=0.0)
    assert dt_obs[0] > _dt, (
        f"first step ({dt_obs[0]:.0f} s) should exceed original dt ({_dt:.0f} s)"
    )
    assert dt_obs[1] < dt_obs[0], (
        f"second step ({dt_obs[1]:.0f} s) should be smaller than "
        + f"first ({dt_obs[0]:.0f} s)"
    )

    listing = (test.workspace / f"{name}.lst").read_text()
    assert "COURANT NUMBER FOR EACH REACH" in listing, (
        "Courant number table not found in listing file"
    )


# ---------------------------------------------------------------------------
# Case 1: ATS_COURANT without ATS in TDIS (warning test)
# ---------------------------------------------------------------------------


def _build_no_ats(test):
    """Build a model with ATS_COURANT but no ATS package in TDIS."""
    name = cases[1]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
        version="mf6",
        exe_name="mf6",
    )
    # Ten fixed time steps; no ATS block.
    flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=_nsteps,
        perioddata=[(_dt, 1, 1.0)] * _nsteps,
    )

    gwf = _build_gwf(sim, name)

    rtp = 20.0
    pak_data = [(0, -1, -1, -1, _dx, _B, _slope, rtp, 1.0, 0.0, _roughness, 0, 0.0, 0)]
    sfr_perioddata = {i: [(0, "inflow", float(q))] for i, q in enumerate(_flows)}

    flopy.mf6.ModflowGwfsfr(
        gwf,
        print_flows=True,
        storage=True,
        ats_courant=1.0,
        maximum_depth_change=1.0e-9,
        nreaches=1,
        packagedata=pak_data,
        connectiondata=[(0,)],
        perioddata=sfr_perioddata,
        pname="sfr-1",
    )

    return sim


def _check_no_ats(test):
    """Verify that the expected warning appears in mfsim.lst."""
    mfsim_lst = (test.workspace / "mfsim.lst").read_text()
    warn_tag = "ATS_COURANT IS SPECIFIED IN THE SFR OPTIONS BLOCK"
    assert warn_tag in mfsim_lst, (
        f"Expected warning '{warn_tag}' not found in mfsim.lst"
    )


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def build_models(idx, test):
    if idx == 1:
        return _build_no_ats(test)
    return _build_ats_active(test)


def check_output(idx, test):
    if idx == 1:
        _check_no_ats(test)
        return
    _check_ats_active(test)


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
