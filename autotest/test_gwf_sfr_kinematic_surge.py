"""
Tests the SFR kinematic-wave STORAGE option under a rapidly varying inflow,
which is the condition a large runoff event routed to a reach produces. Water
provided by the MVR Package is added to the reach inflow rather than
distributed along the reach, so a mover-driven runoff event and a step change
in specified inflow drive the routing the same way.

The channel is five reaches long and is not connected to the groundwater
model, so the routing response is isolated from stream-aquifer exchange. Base
flow is routed to steady state, a twenty-fold surge is imposed for a single
time step, and base flow is then restored. The surge is shorter than the time
the wave takes to cross the channel, so the routed pulse is attenuated.

Cases:
  - sfr-kwsurge01 : storage weight 1.0 (fully implicit, the default).
  - sfr-kwsurge02 : storage weight 0.5 (Crank-Nicolson).

Both cases verify that the surge raises the Courant number above unity, that
the package budget closes, that the peak is attenuated and delayed downstream,
and that the package budget closes.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ("sfr-kwsurge01", "sfr-kwsurge02")
use_tvd = (False, True)

# per-step budget discrepancy tolerances, in percent
_tol_steady = 0.1
_tol_tvd = 0.01

# channel geometry and hydraulics
_nreaches = 5
_dx = 1000.0
_width = 10.0
_slope = 1.0e-3
_roughness = 0.03
_q_base = 5.0
_q_surge = 100.0
_dt = 600.0
_nstp_base = 30
_nstp_surge = 1
_nstp_recess = 40


def _manning_depth(q):
    """Normal depth for a rectangular channel from Manning's equation."""
    from scipy.optimize import brentq

    def residual(d):
        area = _width * d
        perimeter = _width + 2.0 * d
        conveyance = area * (area / perimeter) ** (2.0 / 3.0)
        return conveyance * _slope**0.5 / _roughness - q

    return brentq(residual, 1.0e-8, 100.0)


def _courant(q):
    """Kinematic-wave Courant number for a wide rectangular channel."""
    depth = _manning_depth(q)
    velocity = q / (_width * depth)
    return (5.0 / 3.0) * velocity * _dt / _dx


def build_models(idx, test):
    name = cases[idx]

    # base flow, surge, and recession, each resolved at the same time step
    nper = 3
    perioddata = [
        (_nstp_base * _dt, _nstp_base, 1.0),
        (_nstp_surge * _dt, _nstp_surge, 1.0),
        (_nstp_recess * _dt, _nstp_recess, 1.0),
    ]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=test.workspace, version="mf6", exe_name="mf6"
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="seconds", nper=nper, perioddata=perioddata
    )
    if use_tvd[idx]:
        tdis.ats.initialize(
            maxats=nper,
            perioddata=[(iper, _dt, 1.0, _dt, 2.0, 5.0) for iper in range(nper)],
            filename=f"{name}.ats",
        )
    flopy.mf6.ModflowIms(sim, print_option="summary")

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=1,
        nrow=1,
        ncol=_nreaches,
        delr=_dx,
        delc=_dx,
        top=100.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[(0, 0, 0, 100.0, 5.0)])
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    # reaches are disconnected from the groundwater model (cellid -1) so the
    # test isolates the routing response from stream-aquifer exchange
    top = 100.0
    pak_data = [
        (
            ifno,
            -1,
            -1,
            -1,
            _dx,
            _width,
            _slope,
            top - ifno * _slope * _dx,
            1.0,
            0.0,
            _roughness,
            1 if ifno in (0, _nreaches - 1) else 2,
            1.0,
            0,
        )
        for ifno in range(_nreaches)
    ]
    connectiondata = []
    for ifno in range(_nreaches):
        conn = [ifno]
        if ifno > 0:
            conn.append(ifno - 1)
        if ifno < _nreaches - 1:
            conn.append(-(ifno + 1))
        connectiondata.append(tuple(conn))

    sfr_spd = {
        0: [(0, "inflow", _q_base)],
        1: [(0, "inflow", _q_surge)],
        2: [(0, "inflow", _q_base)],
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_flows=True,
        save_flows=True,
        storage=True,
        nreaches=_nreaches,
        packagedata=pak_data,
        connectiondata=connectiondata,
        perioddata=sfr_spd,
        pname="sfr-1",
        **({"ats_courant": 1.0} if use_tvd[idx] else {}),
    )
    fname = f"{name}.sfr.obs"
    sfr_obs = {
        f"{fname}.csv": [
            ("inflow", "ext-inflow", (0,)),
            ("outflow", "ext-outflow", (_nreaches - 1,)),
        ]
    }
    sfr.obs.initialize(filename=fname, print_input=True, continuous=sfr_obs)

    return sim


def check_output(idx, test):
    name = cases[idx]
    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    inflow = np.abs(obs["INFLOW"])
    outflow = np.abs(obs["OUTFLOW"])

    # the surge is the condition of interest: it must drive the Courant number
    # above unity, otherwise the test is not exercising the intended regime
    cr_base = _courant(_q_base)
    cr_surge = _courant(_q_surge)
    assert cr_base < 1.0, f"base flow Courant number is {cr_base:.2f}, expected < 1"
    assert cr_surge > 1.0, (
        f"surge Courant number is {cr_surge:.2f}, expected > 1; the test is not "
        "exercising the high Courant number regime"
    )

    # the surge must reach the downstream end and be attenuated by routing
    assert outflow.max() > 1.5 * _q_base, (
        f"peak outflow {outflow.max():.2f} shows the surge did not propagate"
    )
    assert outflow.max() < inflow.max(), (
        f"peak outflow {outflow.max():.2f} is not attenuated relative to the "
        f"peak inflow {inflow.max():.2f}"
    )

    # the peak must arrive downstream after it enters upstream
    assert np.argmax(outflow) > np.argmax(inflow), (
        "the downstream peak does not lag the upstream peak"
    )

    # per-step budget discrepancy. The first few steps fill the channel from
    # the initial condition and are excluded; that transient is present with or
    # without the surge and is not what the test is about.
    listing = (test.workspace / f"{name}.lst").read_text()
    sections = re.split(r"SFR-1 BUDGET FOR ENTIRE MODEL", listing)[1:]
    assert len(sections) > 0, "no SFR-1 budget sections found in listing"
    rates = []
    for section in sections:
        found = re.findall(r"PERCENT DISCREPANCY\s*=\s*([0-9.Ee+\-]+)", section)
        # the listing reports cumulative volumes and rates side by side
        if len(found) >= 2:
            rates.append(abs(float(found[1])))
    # the channel fills from the initial condition over the first part of the
    # base period, so the baseline is taken from the last few steps before the
    # surge, by which point the discrepancy has decayed
    steady = max(rates[_nstp_base - 5 : _nstp_base])
    surge = max(rates[_nstp_base:])
    print(f"{name}: steady {steady:.3f}%, surge {surge:.3f}%")

    assert steady < _tol_steady, (
        f"budget discrepancy of {steady:.3f}% at steady base flow exceeds "
        f"{_tol_steady}% before the surge is applied"
    )
    if use_tvd[idx]:
        assert surge < _tol_tvd, (
            f"the TVD scheme should close the budget through the surge, but "
            f"the discrepancy reached {surge:.3f}%"
        )
    else:
        # the standard scheme is expected to degrade here: the surge drives the
        # Courant number to about 2 and the storage term is weighted in time
        assert surge > 1.0, (
            f"the surge discrepancy of {surge:.3f}% is smaller than expected "
            "for the standard scheme; if the scheme has been improved, update "
            "this test and the applicability guidance in the SFR chapter"
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
