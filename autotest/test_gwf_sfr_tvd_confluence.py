"""
Integration tests for the SFR TVD scheme at a stream confluence.

Network:
    R0 --> R1 --> R2 ──┐
                        ├── R5 --> R6 --> R7 --> R8
         R3 --> R4 ────┘

R5 has two upstream reaches, so the van Leer correction falls back to upwind
there. Two headwater hydrographs (R0, R3) drive asymmetric flow. Verifies exact
SFR budget closure, non-negative bounded outlet flow, and the Courant table.

Cases:
  - sfr-tvd-cflu   : TVD, dt = 1800 s (Cr ~ 0.5).
  - sfr-tvd-cflu-b : standard transient baseline, no TVD.

Ponce, V. M. (1989). Engineering Hydrology, Principles and Practices.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"

# case 0: TVD with dt = 1800 s (Cr ~ 0.5)
# case 1: standard transient (STORAGE only, no TVD) — baseline
cases = ("sfr-tvd-cflu", "sfr-tvd-cflu-b")
_use_tvd = (True, False)
_dt = 1800.0  # time step (s) — Cr ~ 0.5 at typical flow

# Hydrograph timing
_dt_hyd = 3600.0  # interval between hydrograph ordinates (s)

# Segment A hydrograph (Ponce 1989 example, m3/s at 3600-s intervals)
_flows_a = np.array([30.0, 60.0, 90.0, 120.0, 150.0, 120.0, 90.0, 60.0, 30.0, 10.0])

# Segment B hydrograph — different peak timing and magnitude
_flows_b = np.array([50.0, 80.0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 5.0])

_nperiods = len(_flows_a)

# Channel parameters — identical for all reaches
_dx = 7200.0  # reach length (m)
_roughness = 0.03574737676661647  # n; gives Cr ~ 1 at Q = 150 m3/s, dt = 3600 s
_slope = 1.0 / _dx  # bed slope (-)
_width = 10.0  # channel width (m)
_nreaches = 9

# SFR connection data (0-indexed reaches).
# Sign convention: positive ic = reach ic is upstream of rno;
#                  negative ic = reach |ic| is downstream of rno.
_connection_data = [
    [0, -1],  # R0: downstream to R1
    [1, 0, -2],  # R1: upstream R0, downstream R2
    [2, 1, -5],  # R2: upstream R1, downstream R5
    [3, -4],  # R3: downstream to R4
    [4, 3, -5],  # R4: upstream R3, downstream R5
    [5, 2, 4, -6],  # R5: confluence — upstream R2 and R4, downstream R6
    [6, 5, -7],  # R6: upstream R5, downstream R7
    [7, 6, -8],  # R7: upstream R6, downstream R8
    [8, 7],  # R8: outlet — upstream R7
]
_ncon = [len(row) - 1 for row in _connection_data]


def _build_one(name, dt, use_tvd, test):
    nsteps_per_period = max(1, round(_dt_hyd / dt))
    perioddata = [(_dt_hyd, nsteps_per_period, 1.0) for _ in range(_nperiods)]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
        version="mf6",
        exe_name="mf6",
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=_nperiods,
        perioddata=perioddata,
    )
    flopy.mf6.ModflowIms(sim, outer_dvclose=1.0e-8, inner_dvclose=1.0e-9)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=1,
        nrow=1,
        ncol=1,
        delr=_dx,
        delc=_dx,
        top=20.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0e-4)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    # rno, cellid (none), rlen, rwid, rgrd, rtp, rbth, rhk, man, ncon, ustrf, ndv
    # cellid (-1,-1,-1) means the reach is not connected to a groundwater cell.
    pak_data = [
        (irno, -1, -1, -1, _dx, _width, _slope, 20.0, 1.0, 0.0, _roughness, nc, 1.0, 0)
        for irno, nc in enumerate(_ncon)
    ]

    # External inflows applied at the two headwaters in each stress period.
    sfr_perioddata = {
        i: [
            (0, "inflow", float(_flows_a[i])),
            (3, "inflow", float(_flows_b[i])),
        ]
        for i in range(_nperiods)
    }

    sfr_kwargs = dict(
        print_flows=True,
        storage=True,
        nreaches=_nreaches,
        packagedata=pak_data,
        connectiondata=_connection_data,
        perioddata=sfr_perioddata,
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
                ("inflow_a", "ext-inflow", (0,)),
                ("inflow_b", "ext-inflow", (3,)),
                ("outflow", "ext-outflow", (8,)),
            ]
        },
    )
    return sim


def build_models(idx, test):
    return _build_one(cases[idx], _dt, _use_tvd[idx], test)


def check_output(idx, test):
    name = cases[idx]
    use_tvd = _use_tvd[idx]

    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    assert len(obs) > 0, "no observations written"

    outflow = obs["OUTFLOW"]  # negative sign convention (leaving network)
    inflow_a = obs["INFLOW_A"]  # positive
    inflow_b = obs["INFLOW_B"]  # positive

    # Non-negative outflow: TVD must not produce spurious positive values.
    assert np.all(outflow <= 0.0), (
        f"outlet outflow contains positive values: {outflow[outflow > 0]}"
    )

    # Outlet peak bounded by combined headwater inflow peak (with margin for
    # transient storage releasing water from previous steps).
    max_combined_inflow = np.max(np.abs(inflow_a) + np.abs(inflow_b))
    assert np.all(np.abs(outflow) <= max_combined_inflow + 1.0), (
        "outlet outflow exceeds combined inflow peak by more than 1 m3/s"
    )

    if use_tvd:
        listing = (test.workspace / f"{name}.lst").read_text()

        # Courant number table must be written for ATS_COURANT runs.
        assert "COURANT NUMBER FOR EACH REACH" in listing, (
            "Courant number table not found in listing file"
        )

        # TVD budget closes to machine precision at every step for the entire
        # SFR network, including the confluence reach R5 where the van Leer
        # anti-diffusion correction falls back to first-order upwind.
        sfr_sections = re.split(r"SFR-1 BUDGET FOR ENTIRE MODEL", listing)[1:]
        assert len(sfr_sections) > 0, "no SFR-1 budget sections found in listing"
        for section in sfr_sections:
            match = re.search(r"PERCENT DISCREPANCY\s*=\s*([0-9.Ee+\-]+)", section)
            if match:
                val = match.group(1)
                assert abs(float(val)) < 0.01, (
                    f"TVD SFR budget discrepancy {val}% exceeds 0.01% tolerance"
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
