"""
Tests for the LAK IMPLICIT option, which solves the lake stage as an unknown in
the groundwater flow matrix instead of by the legacy substitution solver. Each
case is built and run twice -- once with the legacy solver and once with the
IMPLICIT option -- and the converged heads and lake stages must agree.

The cases focus on the dry-lake handling:
  * "drylake"  -- a table (embedded) lake pulled toward its bottom by strong
                  evaporation in steady state, exercising the near-bottom loss
                  limiting and the small diagonal term that keeps the lake row
                  solvable as conductance and area drop to zero.
  * "drain"    -- a transient table lake drained to (and held at) its bottom by
                  strong evaporation with the aquifer drawn down below the lake,
                  exercising lak_nur (which keeps the stage at or above the lake
                  bottom under NEWTON UNDER_RELAXATION) and the small diagonal
                  term that keeps the lake row solvable.
  * "fill"     -- a transient lake filled by rainfall with mixed vertical and
                  horizontal connections crossing the wet/dry transition (the
                  lakebotfill scenario), exercising the smoothed seepage
                  coupling (sRamp) for both connection types.

A separate test ("perched") covers a steady-state lake completely disconnected
from the aquifer: the smoothed stage-coupling keeps the perched leakage on the
lake-stage diagonal, so the implicit formulation now converges for this case
and matches the legacy substitution solver.
"""

import os

import flopy
import numpy as np
import pytest

cases = ["drylake", "drain", "fill"]


def _enable_implicit(lak_file, force_fallback=False):
    # add the IMPLICIT keyword to the LAK OPTIONS block (flopy does not yet
    # expose it as a keyword, so it is written directly to the input file).
    # force_fallback also adds the hidden DEV_FORCE_FALLBACK option, which routes
    # every active lake through the substitution fallback path.
    lines = open(lak_file).read().splitlines()
    out = []
    done = False
    for ln in lines:
        out.append(ln)
        if ln.strip().lower() == "begin options" and not done:
            out.append("  IMPLICIT")
            if force_fallback:
                out.append("  DEV_FORCE_FALLBACK")
            done = True
    assert done, f"could not find OPTIONS block in {lak_file}"
    open(lak_file, "w").write("\n".join(out) + "\n")


def _build_drylake(ws, exe):
    name = "lk"
    nlay, nrow, ncol = 1, 15, 15
    delr = delc = 100.0
    top, botm = 50.0, [-50.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-8,
        outer_maximum=500,
        inner_dvclose=1e-10,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=20.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=5.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, steady_state={0: True})
    chd = [[(0, i, 0), 24.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 18.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    # single embedded vertical connection with a bathymetry table
    conn = [[0, 0, (0, 7, 7), "embeddedv", 10.0, 0.0, 0.0, 1.0, 0.0]]
    pkd = [[0, 26.0, 1]]
    # strong evaporation pulls the lake stage down toward the lake bottom (20)
    perdata = [[0, "RAINFALL", 1.0e-3], [0, "EVAPORATION", 5.0e-2]]
    lak_tab = [
        [20.0, 0.0, 0.0, 0.0],
        [22.0, 4.0e3, 4.0e3, 4.0e3],
        [24.0, 1.2e4, 6.0e3, 6.0e3],
        [26.0, 2.4e4, 8.0e3, 8.0e3],
        [30.0, 5.6e4, 1.0e4, 1.0e4],
        [40.0, 1.56e5, 1.0e4, 1.0e4],
    ]
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        ntables=1,
        tables=[0, f"{name}.laktab"],
        surfdep=0.5,
    )
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(lak_tab),
        ncol=4,
        table=lak_tab,
        filename=f"{name}.laktab",
        pname="laktab",
        parent_file=lak,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


def _build_drain(ws, exe):
    # transient table lake that drains toward (and is held at) its bottom by
    # strong evaporation with the aquifer drawn down below the lake -- exercises
    # lak_nur keeps the stage at or above the lake bottom under NEWTON UR
    name = "lk"
    nlay, nrow, ncol = 1, 15, 15
    delr = delc = 100.0
    top, botm = 50.0, [-50.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(2000.0, 20, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-7,
        outer_maximum=500,
        inner_dvclose=1e-9,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=22.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=5.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    # aquifer drawn down below the lake bottom (20) so the lake drains out
    chd = [[(0, i, 0), 16.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 14.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    conn = [[0, 0, (0, 7, 7), "embeddedv", 10.0, 0.0, 0.0, 1.0, 0.0]]
    pkd = [[0, 30.0, 1]]
    perdata = [[0, "RAINFALL", 1.0e-4], [0, "EVAPORATION", 8.0e-2]]
    lak_tab = [
        [20.0, 0.0, 0.0, 0.0],
        [22.0, 4.0e3, 4.0e3, 4.0e3],
        [24.0, 1.2e4, 6.0e3, 6.0e3],
        [26.0, 2.4e4, 8.0e3, 8.0e3],
        [30.0, 5.6e4, 1.0e4, 1.0e4],
        [40.0, 1.56e5, 1.0e4, 1.0e4],
    ]
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        ntables=1,
        tables=[0, f"{name}.laktab"],
        surfdep=0.5,
    )
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(lak_tab),
        ncol=4,
        table=lak_tab,
        filename=f"{name}.laktab",
        pname="laktab",
        parent_file=lak,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


def _build_perched(ws, exe):
    # a steady-state lake that is completely disconnected from the aquifer: the
    # water table (set by CHD) is held well below the lakebed, so every lake-GWF
    # connection is perched. the smoothed stage-coupling keeps the perched
    # leakage on the lake-stage diagonal, so the lake-stage row stays solvable
    # and the IMPLICIT formulation converges to the same stage the legacy
    # substitution solver finds from the lake volume-stage relation.
    name = "lk"
    nlay, nrow, ncol = 1, 15, 15
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-7,
        outer_maximum=500,
        inner_dvclose=1e-9,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=-40.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, steady_state={0: True})
    # water table held well below the flat lakebed (cell top = 0)
    chd = [[(0, i, 0), -35.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), -45.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    lr, lc = range(6, 9), range(6, 9)
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
    # balanced rainfall and evaporation -- the lake budget is otherwise empty
    perdata = [[0, "RAINFALL", 1.0e-3], [0, "EVAPORATION", 1.0e-3]]
    flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        surfdep=0.5,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


def _build_fill(ws, exe):
    # a lake filled by rainfall and connected to the aquifer by a vertical
    # (bottom) connection and four horizontal (side) connections -- the
    # lakebotfill scenario, with mixed connection types crossing the wet/dry
    # transition. Built under NEWTON (no REWET cell drying), where the legacy and
    # implicit formulations converge to the same solution.
    name = "lk"
    nlay, nrow, ncol = 2, 5, 5
    delr = delc = 100.0
    top, botm = 10.0, [0.0, -50.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(50.0, 10, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-7,
        outer_maximum=500,
        inner_dvclose=1e-9,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=5.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    # gentle aquifer gradient straddling the lake bottom (elevation 0)
    chd = [[(0, i, 0), 3.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    # one vertical (bottom) connection and four horizontal (side) connections
    conn = [[0, 0, (1, 2, 2), "vertical", 1.0, 0.0, 0.0, 0.0, 0.0]]
    sides = [(0, 1, 2), (0, 3, 2), (0, 2, 1), (0, 2, 3)]
    for k, (lay, r, c) in enumerate(sides, start=1):
        conn.append([0, k, (lay, r, c), "horizontal", 1.0, 0.0, 10.0, 50.0, 100.0])
    pkd = [[0, 1.0, len(conn)]]
    # rainfall fills the lake from its initial low stage
    perdata = [[0, "RAINFALL", 5.0e-2], [0, "EVAPORATION", 0.0]]
    flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        surfdep=0.5,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


_BUILDERS = {
    "drylake": _build_drylake,
    "drain": _build_drain,
    "fill": _build_fill,
}


def _run(sim):
    success, _ = sim.run_simulation(silent=True)
    return success


def _heads(ws, name):
    return flopy.utils.HeadFile(os.path.join(ws, f"{name}.hds")).get_data()


def _stage(ws, name):
    sf = flopy.utils.HeadFile(os.path.join(ws, f"{name}.lak.stage.bin"), text="STAGE")
    return float(sf.get_data().flatten()[0])


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    builder = _BUILDERS[name]
    sim_l, mname = builder(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = builder(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    ok_l = _run(sim_l)
    ok_i = _run(sim_i)

    # the implicit formulation must converge for these dry-lake cases (a small
    # diagonal term keeps the lake row solvable and lak_nur holds the stage at
    # the lake bottom)
    assert ok_i, f"IMPLICIT option failed to converge for case '{name}'"

    # where the legacy solver also converges, the implicit result must match it
    if ok_l:
        hl = _heads(str(ws_legacy), mname)
        hi = _heads(str(ws_impl), mname)
        hl = np.where(np.abs(hl) < 1e29, hl, np.nan)
        hi = np.where(np.abs(hi) < 1e29, hi, np.nan)
        maxdiff = float(np.nanmax(np.abs(hl - hi)))
        assert maxdiff < 1e-3, f"head mismatch for '{name}': {maxdiff}"

        sl = _stage(str(ws_legacy), mname)
        si = _stage(str(ws_impl), mname)
        assert abs(sl - si) < 1e-3, f"stage mismatch for '{name}': {sl} vs {si}"


@pytest.mark.developmode
def test_perched_disconnected_implicit(function_tmpdir, targets):
    # a steady-state lake completely disconnected from the aquifer (every
    # connection perched). The smoothed stage-coupling keeps the perched leakage
    # on the lake-stage diagonal (sRamp), so the IMPLICIT formulation converges
    # and matches the legacy substitution solver.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, mname = _build_perched(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_perched(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    assert _run(sim_l), "legacy solver failed for the perched lake"
    assert _run(sim_i), "IMPLICIT failed to converge for the perched lake"

    hl = np.where(
        np.abs(_heads(str(ws_legacy), mname)) < 1e29,
        _heads(str(ws_legacy), mname),
        np.nan,
    )
    hi = np.where(
        np.abs(_heads(str(ws_impl), mname)) < 1e29, _heads(str(ws_impl), mname), np.nan
    )
    maxdiff = float(np.nanmax(np.abs(hl - hi)))
    assert maxdiff < 1e-3, f"perched head mismatch: {maxdiff}"


def _build_weak(ws, exe):
    # a perched lake with a very small lakebed leakance and net rainfall inflow
    # and no outlet. The water table (CHD) is held below the lakebed, so the lake
    # sheds its inflow only by downward leakage; with the small leakance the stage
    # must rise far to drive that leakage. This is the weakly connected regime the
    # substitution fallback exists for, so it is used to exercise the fallback
    # assembly (DEV_FORCE_FALLBACK) against the legacy formulation.
    name = "lk"
    nlay, nrow, ncol = 1, 15, 15
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-7,
        outer_maximum=1000,
        inner_dvclose=1e-9,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=-40.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, steady_state={0: True})
    chd = [[(0, i, 0), -35.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), -45.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    lr, lc = range(6, 9), range(6, 9)
    # small lakebed leakance -> a weakly connected lake-stage equation
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0e-5, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
    # net inflow (rainfall, no evaporation) with no outlet
    perdata = [[0, "RAINFALL", 1.0e-3], [0, "EVAPORATION", 0.0]]
    flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        surfdep=0.5,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


@pytest.mark.developmode
def test_fallback_matches_legacy(function_tmpdir, targets):
    # a weakly connected lake solved three ways, all of which must agree:
    #   1. the legacy substitution solver,
    #   2. the IMPLICIT formulation, and
    #   3. the IMPLICIT formulation with every lake forced onto the substitution
    #      fallback (DEV_FORCE_FALLBACK).
    # case 3 routes the lake through the fallback assembly in lak_fc_implicit
    # (solve the stage by substitution, then assemble it like a constant-stage
    # lake), which must reproduce the legacy result. A small synthetic model does
    # not stall the implicit solver on its own, so the fallback path is forced
    # here to give a deterministic regression test of that assembly.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_fb = function_tmpdir / "fallback"
    for ws in (ws_legacy, ws_impl, ws_fb):
        ws.mkdir()

    sim_l, mname = _build_weak(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_weak(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))
    sim_f, _ = _build_weak(str(ws_fb), exe)
    sim_f.write_simulation(silent=True)
    _enable_implicit(str(ws_fb / f"{mname}.lak"), force_fallback=True)

    assert _run(sim_l), "legacy solver failed for the weak lake"
    assert _run(sim_i), "IMPLICIT failed to converge for the weak lake"
    assert _run(sim_f), "IMPLICIT with forced fallback failed for the weak lake"

    hl = _heads(str(ws_legacy), mname)
    for label, ws in (("implicit", ws_impl), ("fallback", ws_fb)):
        hx = _heads(str(ws), mname)
        maxdiff = float(np.nanmax(np.abs(hl - hx)))
        assert maxdiff < 1e-6, f"{label} head mismatch vs legacy: {maxdiff}"
        sl = _stage(str(ws_legacy), mname)
        sx = _stage(str(ws), mname)
        assert abs(sl - sx) < 1e-6, f"{label} stage mismatch: {sl} vs {sx}"
