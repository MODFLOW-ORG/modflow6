"""
Unit tests for the SFR TVD (van Leer) flux limiter.

Single-reach Ponce (1989) Ex. 9-3 channel. Verifies exact SFR budget closure,
non-negative bounded outflow, and that the Courant table is written.

Cases:
  - sfr-tvd-cr1  : TVD, dt = 3600 s (Cr ~ 1, correction near-zero).
  - sfr-tvd-base : standard transient baseline, no TVD.
  - sfr-tvd-cr05 : TVD, dt = 1800 s (Cr ~ 0.5, van Leer active).

Ponce, V. M. (1989). Engineering Hydrology, Principles and Practices.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"

# case 0: TVD with dt = 3600 s (Cr ~ 1)
# case 1: standard transient with dt = 3600 s (baseline, no TVD)
# case 2: TVD with dt = 1800 s (Cr ~ 0.5, van Leer correction active)
cases = ("sfr-tvd-cr1", "sfr-tvd-base", "sfr-tvd-cr05")
_use_tvd = (True, False, True)
_dts = (3600.0, 3600.0, 1800.0)

# Ponce (1989) Ex 9-3 hydrograph (m3/s, 3600 s intervals)
_dt_hyd = 3600.0
_flows = np.array([30.0, 60.0, 90.0, 120.0, 150.0, 120.0, 90.0, 60.0, 30.0, 10.0])
_nsteps_hyd = len(_flows)

# Channel parameters
_dx = 7200.0
_B = 10.0
_slope = 1.0 / _dx
# Roughness giving Cr = 1 at peak inflow using MODFLOW 6's simplified P=B
# formula (rh = depth, not B*d/(B+2*d)).
# From Cr=1: K^(3/5) * Q_peak^(2/5) = (3/5)*B*(dx/dt),  K = (1/n)*B*sqrt(S)
_roughness = (
    _B
    * _slope**0.5
    / (((3.0 / 5) * _B * _dx / _dt_hyd / max(_flows) ** (2.0 / 5)) ** (5.0 / 3))
)


def _q_from_depth(d):
    """Manning's flow using MODFLOW 6's P=B formula (rh = depth)."""
    if d <= 0.0:
        return 0.0
    return (1.0 / _roughness) * _B * d ** (5.0 / 3.0) * _slope**0.5


def _depth_from_q(Q):
    """Invert Manning's equation by bisection for a rectangular channel."""
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


def _build_one(name, dt, use_tvd, test):
    nsteps_per_period = max(1, round(_dt_hyd / dt))
    nper = _nsteps_hyd
    perioddata = [(_dt_hyd, nsteps_per_period, 1.0) for _ in range(nper)]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
        version="mf6",
        exe_name="mf6",
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=nper,
        perioddata=perioddata,
    )
    flopy.mf6.ModflowIms(sim, outer_dvclose=1.0e-8, inner_dvclose=1.0e-9)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=1,
        nrow=1,
        ncol=6,
        delr=_dx,
        delc=_dx,
        top=20.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0e-4)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    rtp = 20.0
    d0 = _depth_from_q(_flows[0])
    pak_data = [(0, -1, -1, -1, _dx, _B, _slope, rtp, 1.0, 0.0, _roughness, 0, 0.0, 0)]

    sfr_kwargs = dict(
        print_flows=True,
        storage=True,
        nreaches=1,
        packagedata=pak_data,
        connectiondata=[(0,)],
        perioddata={i: [(0, "inflow", float(q))] for i, q in enumerate(_flows)},
        initialstages=[(0, rtp + d0)],
        pname="sfr-1",
    )
    if use_tvd:
        sfr_kwargs["ats_courant"] = 1.0
        sfr_kwargs["maximum_depth_change"] = 1.0e-9

    sfr = flopy.mf6.ModflowGwfsfr(gwf, **sfr_kwargs)

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


def build_models(idx, test):
    return _build_one(cases[idx], _dts[idx], _use_tvd[idx], test)


def check_output(idx, test):
    name = cases[idx]
    use_tvd = _use_tvd[idx]

    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    assert len(obs) > 0, "no observations written"

    outflow = obs["OUTFLOW"]  # negative sign convention
    inflow = obs["INFLOW"]  # positive sign convention

    # non-negative outflow: TVD should not produce spurious positive values
    assert np.all(outflow <= 0.0), (
        f"outflow contains positive values: {outflow[outflow > 0]}"
    )

    # outflow magnitude bounded by inflow peak (no overshooting)
    assert np.all(np.abs(outflow) <= np.max(np.abs(inflow)) + 1.0), (
        "outflow magnitude exceeds inflow peak by more than 1 m3/s"
    )

    if use_tvd:
        listing = (test.workspace / f"{name}.lst").read_text()

        # Courant number table is written for ATS_COURANT runs
        assert "COURANT NUMBER FOR EACH REACH" in listing, (
            "Courant number table not found in listing file"
        )

        # TVD budget closes exactly by construction (volume balance).
        # Parse only the SFR-1 package budget sections; the GWF model budget
        # reports 200% when both TOTAL IN and TOTAL OUT are zero (0/0 artifact).
        sfr_sections = re.split(r"SFR-1 BUDGET FOR ENTIRE MODEL", listing)[1:]
        assert len(sfr_sections) > 0, "no SFR-1 budget sections found in listing"
        for section in sfr_sections:
            match = re.search(r"PERCENT DISCREPANCY\s*=\s*([0-9.Ee+\-]+)", section)
            if match:
                val = match.group(1)
                assert abs(float(val)) < 0.01, (
                    f"TVD SFR budget error {val}% exceeds 0.01% tolerance"
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
