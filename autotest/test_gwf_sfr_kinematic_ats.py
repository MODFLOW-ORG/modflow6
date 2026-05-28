"""
Test the ATS_COURANT option in the SFR Package for adaptive time stepping
using the kinematic wave approximation.  The problem is the Ponce (1989)
Example 9-3 routing problem used in test_gwf_sfr_kinematic01.py.

The simulation uses a single stress period of 36000 s (10 x 3600 s) so that
ATS can sub-step freely within the period.  Inflow varies via a time series
with stepwise interpolation matching the original 3600 s hydrograph intervals.
ATS is started with dt0 = 14400 s (4x the Courant-optimal step) so that the
first step is deliberately over-large; ATS_COURANT then submits Courant-based
time step requests (targeting Cr = 1) for all subsequent steps, driving them
down to roughly the Courant-optimal interval.

The test verifies that:
- multiple sub-steps are taken within the single period (not just 1), and
- the listing file contains the end-of-simulation Courant number table.

Ponce, V. M. (1989). Engineering Hydrology, Principles and Practices.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ("sfr-kwats01",)

# Ponce (1989) Ex 9-3 hydrograph (non-zero flows only; 10 intervals)
_dt = 3600.0  # s, original hydrograph time step
_flows = np.array([30.0, 60.0, 90.0, 120.0, 150.0, 120.0, 90.0, 60.0, 30.0, 10.0])
_nsteps = len(_flows)  # 10 hydrograph intervals


def build_models(idx, test):
    name = cases[idx]

    # Single stress period long enough to hold all hydrograph intervals.
    # One long period lets ATS take multiple sub-steps instead of being capped
    # by a short perlen on every new period.
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
    dtmin = 30.0
    dtmax = dt0
    dtadj = 2.0
    dtfailadj = 5.0
    # iperats is 0-indexed in flopy (it adds 1 when writing), so 0 maps to
    # stress period 1 in the MODFLOW 6 input.
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, dt0, dtmin, dtmax, dtadj, dtfailadj)],
        filename=f"{name}.ats",
    )

    # spatial discretization
    nlay, nrow, ncol = 1, 1, 6
    dx = 7200.0  # m
    top = 20.0
    botm = 0.0
    h_left = 11.0  # CHD head at column 1 (left)
    h_right = 9.0  # CHD head at column 6 (right)
    h_init = 10.0  # initial head between CHD values

    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-9,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=dx,
        delc=dx,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0e-4)
    flopy.mf6.ModflowGwfic(gwf, strt=h_init)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    # CHDs in columns 1 and 6 drive left-to-right flow through the aquifer.
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            (0, 0, 0, h_left),
            (0, 0, ncol - 1, h_right),
        ],
    )
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    # SFR: single reach, Manning roughness tuned to give Cr ~ 1 at dt = _dt
    slope = 1.0 / dx
    roughness = 0.03574737676661647
    pak_data = [(0, -1, -1, -1, dx, 10.0, slope, top, 1.0, 0.0, roughness, 0, 0.0, 0)]

    # Inflow is set via a time series so it varies within the single long period.
    # With stepwise interpolation each value holds for one 3600 s interval.
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
        pname="sfr-1",
    )
    # Terminal sentinel at perlen ensures the time series covers the full
    # simulation interval (MODFLOW 6 requires this for stepwise interpolation).
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


def check_output(idx, test):
    name = cases[idx]

    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    assert len(obs) > 0, "no observations written"

    # With dt0 = 4 * _dt and ATS_COURANT = 1, the first step takes 14400 s
    # and subsequent steps are driven down to ~ _dt = 3600 s.  The total
    # number of steps in the 36000 s period should be well above 1.
    assert len(obs) >= 5, (
        f"expected >= 5 time steps from ATS sub-stepping, got {len(obs)}"
    )

    # First step should be approximately dt0; subsequent steps should be smaller.
    totim = obs["totim"]
    dt_obs = np.diff(totim, prepend=0.0)
    assert dt_obs[0] > _dt, (
        f"first step ({dt_obs[0]:.0f} s) should exceed original dt ({_dt:.0f} s)"
    )
    assert dt_obs[1] < dt_obs[0], (
        f"second step ({dt_obs[1]:.0f} s) should be smaller than "
        + f"first ({dt_obs[0]:.0f} s)"
    )

    # Listing file must contain the end-of-simulation Courant number table.
    listing = (test.workspace / f"{name}.lst").read_text()
    assert "COURANT NUMBER FOR EACH REACH" in listing, (
        "Courant number table not found in listing file"
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
