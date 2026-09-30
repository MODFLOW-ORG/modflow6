"""PERIOD-block time-series persistence for GWF-UZF, across eight related
scenarios:

  finf_persist   (uzf-finf-persist): a TS-linked FINF value continues to
      be applied in periods that don't reissue it, updating every
      timestep rather than only once per period.
  et_persist     (uzf-et-persist): TS-linked PET/EXTDP/EXTWC values
      continue to be applied in periods that don't reissue them, under
      SIMULATE_ET.
  ts_switch      (uzf-ts-switch): reissuing a cell's FINF in a later
      PERIOD block -- with a different time series, then a literal value
      -- fully replaces the prior time-series link.
  period_partial_ts (uzf-per-part-ts): a TS-linked FINF value on a cell
      keeps tracking its time series in a later period whose PERIOD
      block reissues a different cell but doesn't mention this one.
  period_reorder (uzf-per-reord): PERIOD FINF values are applied to the
      correct cell regardless of the row order they're given in, and a
      later period specifying only a subset of cells leaves the others
      unchanged.
  aux_persist    (uzfauxper): a PERIOD AUX value persists into a later
      period with no PERIOD block at all, the same as other
      advanced-package PERIOD settings.
  aux_dedup      (uzf-aux-dedup): reissuing one AUX column's setting on
      an active cell does not affect a different AUX column on the same
      row -- each column's own time-series link is tracked and cleared
      independently.
  cross_field_dedup (uzf-cross-field): reissuing one field (FINF) on an
      active, SIMULATE_ET cell does not affect a different field's
      (PET/EXTDP/EXTWC) own time-series tracking on the same row.

The first five share a 2-layer, 1-row grid (column count varies with
cell count); aux_persist/aux_dedup use their own single-layer setup.
"""

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = [
    "uzf-finf-persist",
    "uzf-et-persist",
    "uzf-ts-switch",
    "uzf-per-part-ts",
    "uzf-per-reord",
    "uzfauxper",
    "uzf-aux-dedup",
    "uzf-cross-field",
]

PERIOD_AUX1 = 77.0
AUX1_LITERAL = 50.0


def _base_sim(ws, name, nper, tdis_rc, ncol):
    """Return a minimal transient GWF sim ready for a UZF package,
    2-layer/1-row grid with a CHD at layer 2 for each column."""
    nlay, nrow = 2, 1
    top = 10.0
    botm = [0.0, -10.0]

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(
        sim, outer_dvclose=1e-6, inner_dvclose=1e-6, linear_acceleration="BICGSTAB"
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=0.0, sy=0.1, steady_state={0: False}, transient={0: True}
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(1, 0, j), 1.0] for j in range(ncol)]
    )
    return sim, gwf


# ===== finf_persist =====


def _get_finf_persist(ws, name):
    # 3 periods; period 3 has two timesteps to confirm the value updates
    # every timestep, not just at period start
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 2, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 1)

    packagedata = [(0, (0, 0, 0), 1, 0, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5)]
    perioddata = {
        0: [(0, "finf1", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA)],
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        print_flows=True,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    ts_data = [
        (0.0, 0.05),
        (1.0, 0.05),
        (2.0, 0.02),
        (3.0, 0.08),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["finf1"],
        interpolation_methodrecord=["linearend"],
    )
    uzf_obs = {f"{name}.uzf.obs.csv": [("infil1", "infiltration", (0,))]}
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_finf_persist(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[0]}.uzf.obs.csv", delimiter=",", names=True
    )
    infil = obs["INFIL1"]

    # 3 periods x [1, 1, 2] timesteps -> 4 records at t=1.0, 2.0, 2.5, 3.0
    assert len(infil) == 4, f"expected 4 obs records, got {len(infil)}"

    expected = [0.05, 0.02, 0.05, 0.08]
    assert np.allclose(infil, expected), (
        f"infiltration {infil} does not match the time series values "
        f"{expected} expected at t=1.0, 2.0, 2.5, 3.0"
    )


# ===== et_persist =====


def _get_et_persist(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 1)

    packagedata = [(0, (0, 0, 0), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.2, 3.5)]
    perioddata = {
        0: [
            (0, 0.1, "pet1", "extdp1", "extwc1", DNODATA, DNODATA, DNODATA),
        ],
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        simulate_et=True,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    ts_data = [
        (0.0, 0.05, 2.0, 0.15),
        (1.0, 0.05, 2.0, 0.15),
        (2.0, 0.02, 4.0, 0.12),
        (3.0, 0.08, 1.0, 0.18),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["pet1", "extdp1", "extwc1"],
        interpolation_methodrecord=["linearend"] * 3,
    )
    uzf_obs = {f"{name}.uzf.obs.csv": [("uzet1", "uzet", (0,))]}
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_et_persist(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[1]}.uzf.obs.csv", delimiter=",", names=True
    )
    uzet = obs["UZET1"]

    assert len(uzet) == 3, f"expected 3 obs records, got {len(uzet)}"
    # UZET should differ at every period end (t=1, 2, 3), even though
    # PERIOD only reissues PET/EXTDP/EXTWC in period 1
    assert not np.isclose(uzet[0], uzet[1]), (
        f"UZET unchanged between periods 1 and 2: {uzet[0]:.6g} == {uzet[1]:.6g}"
    )
    assert not np.isclose(uzet[1], uzet[2]), (
        f"UZET unchanged between periods 2 and 3: {uzet[1]:.6g} == {uzet[2]:.6g}"
    )


# ===== ts_switch =====


def _get_ts_switch(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 1)

    packagedata = [(0, (0, 0, 0), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5)]
    perioddata = {
        # period 1: FINF linked to finfX
        0: [(0, "finfX", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA)],
        # period 2: reissued, now linked to a different time series, finfY
        1: [(0, "finfY", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA)],
        # period 3: reissued again, now a literal value -- no time series
        2: [(0, 0.15, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA)],
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    # finfX is only referenced in period 1; its values after t=1 (0.99)
    # are deliberately distinct from anything else in this test.
    ts_data = [
        (0.0, 0.05, 0.20),
        (1.0, 0.05, 0.20),
        (2.0, 0.99, 0.20),
        (3.0, 0.99, 0.20),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["finfX", "finfY"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    uzf_obs = {f"{name}.uzf.obs.csv": [("infil1", "infiltration", (0,))]}
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_ts_switch(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[2]}.uzf.obs.csv", delimiter=",", names=True
    )
    infil = obs["INFIL1"]

    assert len(infil) == 3, f"expected 3 obs records, got {len(infil)}"

    # period 1 (t=1): finfX
    assert np.isclose(infil[0], 0.05), infil

    # period 2 (t=2): reissued with finfY, replacing the finfX link
    assert np.isclose(infil[1], 0.20), (
        f"expected finfY's value 0.20 at period 2, got {infil[1]}; a leftover "
        "link to finfX would show up as 0.99 instead"
    )

    # period 3 (t=3): reissued with a literal value, replacing the
    # finfY link
    assert np.isclose(infil[2], 0.15), (
        f"expected the literal value 0.15 at period 3, got {infil[2]}; a "
        "leftover link to finfY would show up as 0.20 instead"
    )


# ===== period_partial_ts =====


def _get_period_partial_ts(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 2)

    packagedata = [
        (i, (0, 0, i), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5) for i in range(2)
    ]
    perioddata = {
        # period 1: both cells given TS-linked FINF
        0: [
            (0, "finfA", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
            (1, "finfB", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
        ],
        # period 2: reissued, but only cell 1 (B) is mentioned -- cell 0
        # (A) is not in this block at all, though the block itself exists
        1: [
            (1, "finfB", DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
        ],
        # period 3: no PERIOD block at all
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    ts_data = [
        (0.0, 0.05, 0.05),
        (1.0, 0.05, 0.05),
        (2.0, 0.02, 0.20),
        (3.0, 0.08, 0.30),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["finfA", "finfB"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    uzf_obs = {
        f"{name}.uzf.obs.csv": [
            ("infila", "infiltration", (0,)),
            ("infilb", "infiltration", (1,)),
        ]
    }
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_period_partial_ts(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[3]}.uzf.obs.csv", delimiter=",", names=True
    )
    infila = obs["INFILA"]
    infilb = obs["INFILB"]

    assert len(infila) == 3, f"expected 3 obs records, got {len(infila)}"

    # cell B is reissued in period 2, then carried with no block in period
    # 3 -- its time series must be followed throughout
    assert np.allclose(infilb, [0.05, 0.20, 0.30]), infilb

    # cell A is never mentioned in periods 2 or 3, but its time series
    # must still be followed throughout, the same as cell B
    assert np.allclose(infila, [0.05, 0.02, 0.08]), (
        f"cell A's time series was not followed: got {infila}, expected "
        "[0.05, 0.02, 0.08]; a later period's PERIOD block reissuing a "
        "different cell must not affect a cell it doesn't mention"
    )


# ===== period_reorder =====


def _get_period_reorder(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 2, tdis_rc, 3)

    packagedata = [
        (i, (0, 0, i), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5) for i in range(3)
    ]
    perioddata = {
        # period 1: all 3 cells given, listed out of order (2, 0, 1), each
        # a distinct literal value
        0: [
            (2, 0.03, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
            (0, 0.01, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
            (1, 0.02, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
        ],
        # period 2: only cell 1 reissued; cells 0 and 2 should keep their
        # period 1 values unchanged
        1: [
            (1, 0.05, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA),
        ],
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    uzf_obs = {
        f"{name}.uzf.obs.csv": [
            ("infil0", "infiltration", (0,)),
            ("infil1", "infiltration", (1,)),
            ("infil2", "infiltration", (2,)),
        ]
    }
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_period_reorder(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[4]}.uzf.obs.csv", delimiter=",", names=True
    )

    assert len(obs["time"]) == 2, f"expected 2 obs records, got {len(obs['time'])}"

    # period 1: order-independent mapping to the correct cell
    assert np.isclose(obs["INFIL0"][0], 0.01), obs["INFIL0"][0]
    assert np.isclose(obs["INFIL1"][0], 0.02), obs["INFIL1"][0]
    assert np.isclose(obs["INFIL2"][0], 0.03), obs["INFIL2"][0]

    # period 2: only cell 1 reissued; cells 0 and 2 persist unchanged
    assert np.isclose(obs["INFIL0"][1], 0.01), (
        f"cell 0 changed to {obs['INFIL0'][1]} despite not being reissued in period 2"
    )
    assert np.isclose(obs["INFIL1"][1], 0.05), obs["INFIL1"][1]
    assert np.isclose(obs["INFIL2"][1], 0.03), (
        f"cell 2 changed to {obs['INFIL2'][1]} despite not being reissued in period 2"
    )


# ===== aux_persist =====


def _get_aux_persist(ws, name):
    nlay, nrow, ncol = 1, 1, 3
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = [-10.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    nuzfcells = ncol
    kv = 1e-4

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=2, perioddata=[(1.0, 5, 1.0)] * 2
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        outer_dvclose=1e-6,
        outer_maximum=600,
        inner_maximum=100,
        inner_dvclose=1e-6,
        rcloserecord=0.1,
        linear_acceleration="BICGSTAB",
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
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=-5.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kv, k33=kv)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-6,
        sy=0.35,
        steady_state={0: False},
        transient={0: True},
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), -5.0], [(0, 0, ncol - 1), -5.0]],
        pname="CHD-1",
    )

    uzf_pkdat = [
        (iuzno, (0, 0, iuzno), 1, -1, 1.0, kv, 0.2, 0.4, 0.3, 3.5)
        for iuzno in range(nuzfcells)
    ]
    # period 1 only -- no PERIOD block for period 2, so FINF/AUX1 should
    # both carry over unchanged
    uzf_spd = {
        0: [
            [iuzno, 1e-6, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, PERIOD_AUX1]
            for iuzno in range(nuzfcells)
        ],
    }
    flopy.mf6.ModflowGwfuzf(
        gwf,
        save_flows=True,
        budget_filerecord=f"{name}.uzf.bud",
        nuzfcells=nuzfcells,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        pname="UZF-1",
        auxiliary=["aux1"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )

    return sim


def _check_aux_persist(test):
    bud_fname = test.workspace / f"{test.name}.uzf.bud"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")

    aux_records = bobj.get_data(text="AUXILIARY")
    assert len(aux_records) == 10, (
        f"Expected 10 AUXILIARY budget records (5/period), got {len(aux_records)}"
    )

    aux1_vals = np.array([rec["AUX1"][0] for rec in aux_records])
    period1 = aux1_vals[:5]
    period2 = aux1_vals[5:]

    assert np.allclose(period1, PERIOD_AUX1), (
        f"Period 1 AUX1 expected {PERIOD_AUX1}, got {period1}"
    )
    assert np.allclose(period2, PERIOD_AUX1), (
        f"Period 2 AUX1 expected {PERIOD_AUX1} (carried over from period 1; "
        f"period 2 has no PERIOD block), got {period2}"
    )


# ===== aux_dedup =====


def _get_aux_dedup(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 1)

    packagedata = [(0, (0, 0, 0), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.1, 3.5)]
    perioddata = {
        # period 1: both AUX columns TS-linked
        0: [
            (
                0,
                0.05,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                "auxA",
                "auxB",
            ),
        ],
        # period 2: row reissued, AUX1 gets a fresh literal; AUX2 is
        # DNODATA and must keep tracking its own time series
        1: [
            (
                0,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                DNODATA,
                AUX1_LITERAL,
                DNODATA,
            ),
        ],
        # period 3: no PERIOD block at all
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        auxiliary=["aux1", "aux2"],
        budget_filerecord=f"{name}.uzf.bud",
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    # auxA's values after t=1 are distinct from AUX1_LITERAL so a leftover
    # link would be obvious; auxB is unrelated and followed every period.
    ts_data = [
        (0.0, 10.0, 100.0),
        (1.0, 10.0, 100.0),
        (2.0, 15.0, 200.0),
        (3.0, 25.0, 300.0),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["auxA", "auxB"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim


def _check_aux_dedup(test):
    bud_fname = test.workspace / f"{cases[6]}.uzf.bud"
    bobj = flopy.utils.CellBudgetFile(bud_fname, precision="double")
    aux_records = bobj.get_data(text="AUXILIARY")

    assert len(aux_records) == 3, (
        f"expected 3 AUXILIARY budget records (1/period), got {len(aux_records)}"
    )

    aux1 = np.array([rec["AUX1"][0] for rec in aux_records])
    aux2 = np.array([rec["AUX2"][0] for rec in aux_records])

    assert np.allclose(aux1, [10.0, 50.0, 50.0]), (
        f"AUX1 expected [10.0, 50.0, 50.0] (TS in period 1, literal "
        f"override from period 2 persisting into period 3), got {aux1}"
    )
    assert np.allclose(aux2, [100.0, 200.0, 300.0]), (
        f"AUX2 expected [100.0, 200.0, 300.0] (its own time series "
        f"followed throughout, unaffected by AUX1's reissue in period "
        f"2), got {aux2}"
    )


# ===== cross_field_dedup =====


def _get_cross_field_dedup(ws, name):
    tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 1, 1.0)]
    sim, gwf = _base_sim(ws, name, 3, tdis_rc, 1)

    packagedata = [(0, (0, 0, 0), 1, -1, 1.0, 1.0, 0.1, 0.3, 0.2, 3.5)]
    perioddata = {
        # period 1: FINF and all three ET fields TS-linked
        0: [(0, "finfZ", "petZ", "extdpZ", "extwcZ", DNODATA, DNODATA, DNODATA)],
        # period 2: row reissued, FINF switched to a literal; PET/EXTDP/EXTWC
        # stay DNODATA and must keep tracking their own time series
        1: [(0, 0.3, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA, DNODATA)],
        # period 3: no PERIOD block at all
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        simulate_et=True,
        ntrailwaves=7,
        nwavesets=40,
        packagedata=packagedata,
        perioddata=perioddata,
        pname="uzf-1",
    )
    # finfZ's values after t=1 (0.9, 0.9) are distinct from the period-2
    # literal (0.3) so a leftover link would be obvious.
    ts_data = [
        (0.0, 0.05, 0.05, 2.0, 0.15),
        (1.0, 0.05, 0.05, 2.0, 0.15),
        (2.0, 0.9, 0.02, 4.0, 0.12),
        (3.0, 0.9, 0.08, 1.0, 0.18),
    ]
    uzf.ts.initialize(
        filename=f"{name}.uzf.ts",
        timeseries=ts_data,
        time_series_namerecord=["finfZ", "petZ", "extdpZ", "extwcZ"],
        interpolation_methodrecord=["linearend"] * 4,
    )
    uzf_obs = {
        f"{name}.uzf.obs.csv": [
            ("infil1", "infiltration", (0,)),
            ("uzet1", "uzet", (0,)),
        ]
    }
    uzf.obs.initialize(filename=f"{name}.uzf.obs", continuous=uzf_obs)

    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_cross_field_dedup(test):
    obs = np.genfromtxt(
        test.workspace / f"{cases[7]}.uzf.obs.csv", delimiter=",", names=True
    )
    infil = obs["INFIL1"]
    uzet = obs["UZET1"]

    assert len(infil) == 3, f"expected 3 obs records, got {len(infil)}"

    assert np.allclose(infil, [0.05, 0.3, 0.3]), (
        f"FINF expected [0.05, 0.3, 0.3] (TS in period 1, literal override "
        f"from period 2 persisting into period 3), got {infil}"
    )
    # UZET must keep differing every period, since PET/EXTDP/EXTWC's own
    # time series is followed throughout, unaffected by FINF's reissue
    assert not np.isclose(uzet[0], uzet[1]), (
        f"UZET unchanged between periods 1 and 2: {uzet[0]:.6g} == {uzet[1]:.6g}"
    )
    assert not np.isclose(uzet[1], uzet[2]), (
        f"UZET unchanged between periods 2 and 3: {uzet[1]:.6g} == {uzet[2]:.6g}"
    )


# ===== dispatch =====

_builders = [
    _get_finf_persist,
    _get_et_persist,
    _get_ts_switch,
    _get_period_partial_ts,
    _get_period_reorder,
    _get_aux_persist,
    _get_aux_dedup,
    _get_cross_field_dedup,
]

_checkers = [
    _check_finf_persist,
    _check_et_persist,
    _check_ts_switch,
    _check_period_partial_ts,
    _check_period_reorder,
    _check_aux_persist,
    _check_aux_dedup,
    _check_cross_field_dedup,
]


def build_models(idx, test):
    name = cases[idx]
    return _builders[idx](test.workspace, name)


def check_output(idx, test):
    _checkers[idx](test)


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
