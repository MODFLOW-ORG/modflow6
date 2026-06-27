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
                  lakebotfill scenario), exercising the lakebed seepage coupling
                  for both connection types.

A separate test ("perched") covers a steady-state lake completely disconnected
from the aquifer: holding the connected-cell head at the lake bottom keeps the
perched leakage on the lake-stage diagonal, so the implicit formulation
converges for this case and matches the legacy substitution solver.

Every implicit run is also checked for budget closure (the model percent
discrepancy must be near zero), which guards the implicit seepage assembly and
its budget reporting (lak_cq) against silently disagreeing.
"""

import os
import re

import flopy
import numpy as np
import pytest

cases = ["drylake", "drain", "fill"]


def _max_discrepancy(ws, name, floor=1e-3):
    # largest absolute percent discrepancy across the budgets reported in the GWF
    # model listing -- the whole-model volume budget and the LAK package
    # sub-budget -- in both the time-step and cumulative columns. Budgets whose
    # total flow is below `floor` are skipped: the percent discrepancy of a
    # near-zero budget (e.g. an essentially empty disconnected lake) is dominated
    # by round-off and is not meaningful.
    text = open(os.path.join(ws, f"{name}.lst")).read()
    worst = None
    for chunk in text.split("BUDGET FOR ENTIRE MODEL")[1:]:
        tin = re.search(r"TOTAL IN\s*=\s*([-+Ee0-9.]+)", chunk)
        disc = re.findall(r"PERCENT DISCREPANCY\s*=\s*([-+Ee0-9.]+)", chunk)
        if tin is None or not disc or abs(float(tin.group(1))) < floor:
            continue
        worst = max([worst or 0.0] + [abs(float(d)) for d in disc])
    return worst


def _assert_budget_closes(ws, name, tol=0.5):
    # the whole-model volume budget and the LAK package budget must both close
    # for the implicit formulation (the package budget covers the stage-driven
    # losses, which the implicit formulation ramps near the lake bottom)
    disc = _max_discrepancy(ws, name)
    assert disc is not None, f"no model budget discrepancy found in {name}.lst"
    assert disc < tol, f"budget discrepancy {disc} exceeds {tol} percent for {ws}"


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

    # the implicit seepage assembly and its budget reporting must close
    _assert_budget_closes(str(ws_impl), mname)

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
    # connection perched). Holding the connected-cell head at the lake bottom
    # keeps the perched leakage on the lake-stage diagonal, so the IMPLICIT
    # formulation converges and matches the legacy substitution solver.
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

    _assert_budget_closes(str(ws_impl), mname)

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


def _build_confined(ws, exe):
    # a lake exchanging with a confined (non-Newton) aquifer. The cells use
    # icelltype=0 and the model has no NEWTON option, so the cells stay fully
    # saturated and there is no REWET drying/rewetting ambiguity; the implicit
    # and legacy formulations must therefore agree. Confirms the implicit lake
    # assembly works without the Newton formulation. The aquifer heads are held
    # above the flat lakebed (cell top = 0) so every connection stays wet.
    name = "lk"
    nlay, nrow, ncol = 1, 11, 11
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
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
    # no newtonoptions -> standard (non-Newton) formulation
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    # confined (icelltype=0): saturated thickness is fixed, cells never dry
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    # gradient with both heads above the lakebed (cell top = 0): connections stay
    # wet, so the lake exchanges freely with the aquifer (not perched)
    chd = [[(0, i, 0), 3.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    lr, lc = range(4, 7), range(4, 7)
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
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


def _build_lkt(ws, exe):
    # a coupled GWF + GWT simulation with an ACTIVE lake (LAK) and lake transport
    # (LKT). The lake starts at concentration 100 and is maintained by a clean
    # (concentration 0) specified inflow while it leaks to the aquifer, so the
    # lake concentration evolves and solute is carried into the aquifer by the
    # lake-aquifer seepage. The lake-aquifer flows that LKT consumes come from the
    # GWF LAK budget, so running this with the legacy and the IMPLICIT LAK and
    # comparing the lake and aquifer concentrations checks that the implicit
    # formulation feeds transport the same flows as the legacy solver.
    gwfname, gwtname = "gwf_lk", "gwt_lk"
    nlay, nrow, ncol = 1, 11, 11
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
    sim = flopy.mf6.MFSimulation(sim_name="lkt", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(100.0, 10, 1.0)])

    # --- GWF model ---
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-8,
        outer_maximum=500,
        inner_dvclose=1e-10,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
        filename=f"{gwfname}.ims",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    sim.register_ims_package(imsgwf, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    # CHD carries a clean (concentration 0) auxiliary so it can be an SSM source
    chd = [[(0, i, 0), 3.0, 0.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0, 0.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chd, auxiliary=["CONCENTRATION"], pname="CHD-1"
    )
    lr, lc = range(4, 7), range(4, 7)
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
    # clean specified inflow maintains the lake while it leaks to the aquifer
    perdata = [[0, "INFLOW", 20.0], [0, "EVAPORATION", 0.0]]
    flopy.mf6.ModflowGwflak(
        gwf,
        save_flows=True,
        print_stage=True,
        stage_filerecord=f"{gwfname}.lak.stage.bin",
        budget_filerecord=f"{gwfname}.lak.cbc",
        nlakes=1,
        noutlets=0,
        packagedata=pkd,
        connectiondata=conn,
        perioddata=perdata,
        pname="LAK-1",
        surfdep=0.5,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # --- GWT model ---
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=1e-8,
        outer_maximum=500,
        inner_dvclose=1e-10,
        inner_maximum=200,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtlkt(
        gwt,
        save_flows=True,
        print_concentration=True,
        concentration_filerecord=f"{gwtname}.lkt.bin",
        budget_filerecord=f"{gwtname}.lkt.cbc",
        packagedata=[[0, 100.0]],
        lakeperioddata=[[0, "STATUS", "ACTIVE"]],
        flow_package_name="LAK-1",
        pname="LKT-1",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return sim, gwfname, gwtname


@pytest.mark.developmode
def test_transport_lkt(function_tmpdir, targets):
    # an active lake coupled to lake transport (LKT). The implicit LAK must feed
    # LKT the same lake-aquifer flows as the legacy solver, so the lake and
    # aquifer concentrations must match between the two formulations.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, gwfname, gwtname = _build_lkt(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _, _ = _build_lkt(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{gwfname}.lak"))

    assert _run(sim_l), "legacy LAK failed in the coupled transport model"
    assert _run(sim_i), "IMPLICIT LAK failed in the coupled transport model"

    _assert_budget_closes(str(ws_impl), gwfname)

    def lake_conc(ws):
        f = flopy.utils.HeadFile(
            os.path.join(ws, f"{gwtname}.lkt.bin"), text="CONCENTRATION"
        )
        return np.array(
            [c.flatten()[0] for c in (f.get_data(totim=t) for t in f.get_times())]
        )

    def aquifer_conc(ws):
        return flopy.utils.HeadFile(
            os.path.join(ws, f"{gwtname}.ucn"), text="CONCENTRATION"
        ).get_data()

    cl, ci = lake_conc(str(ws_legacy)), lake_conc(str(ws_impl))
    assert np.allclose(cl, ci, atol=1e-4), f"lake concentration mismatch: {cl} vs {ci}"
    # the lake must actually exchange solute (concentration evolves from 100)
    assert cl[-1] < 99.0, f"lake concentration did not evolve: {cl[-1]}"

    al, ai = aquifer_conc(str(ws_legacy)), aquifer_conc(str(ws_impl))
    maxdiff = float(np.nanmax(np.abs(al - ai)))
    assert maxdiff < 1e-4, f"aquifer concentration mismatch: {maxdiff}"
    assert float(np.nanmax(ai)) > 1e-3, "no solute reached the aquifer"


def _build_constant_stage(ws, exe):
    # a constant-stage lake (STATUS CONSTANT): the lake stage is held fixed and
    # the aquifer responds to it. The implicit formulation must hold the lake row
    # at the specified stage (unit diagonal, rhs = stage) instead of solving it,
    # and give the same aquifer heads as the legacy solver.
    name = "lk"
    nlay, nrow, ncol = 1, 11, 11
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    chd = [[(0, i, 0), 3.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    lr, lc = range(4, 7), range(4, 7)
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
    # hold the lake at a constant stage of 5.0
    perdata = [[0, "STATUS", "CONSTANT"], [0, "STAGE", 5.0]]
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
def test_constant_stage_lake(function_tmpdir, targets):
    # a constant-stage lake exercises the implicit constant-stage assembly branch
    # (hold the lake row at the specified stage). The aquifer heads and the held
    # stage must match the legacy solver, and the budget must close.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, mname = _build_constant_stage(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_constant_stage(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    assert _run(sim_l), "legacy solver failed for the constant-stage lake"
    assert _run(sim_i), "IMPLICIT failed to converge for the constant-stage lake"

    _assert_budget_closes(str(ws_impl), mname)

    hl = _heads(str(ws_legacy), mname)
    hi = _heads(str(ws_impl), mname)
    maxdiff = float(np.nanmax(np.abs(hl - hi)))
    assert maxdiff < 1e-3, f"constant-stage head mismatch: {maxdiff}"
    # the held stage must be reproduced exactly
    assert abs(_stage(str(ws_impl), mname) - 5.0) < 1e-6, "stage not held constant"


def _build_twolake(ws, exe):
    # two lakes connected by an outlet: an upper lake (lake 0) fed by strong
    # rainfall spills, once its stage exceeds the outlet invert, through a
    # MANNING outlet into a lower lake (lake 1), which sheds the water by leaking
    # to the aquifer. Exercises nlakes>1, noutlets>0, and the lake-to-lake outlet
    # routing (simoutrate), which the implicit formulation lags one outer
    # iteration -- a path the single-lake/no-outlet cases never reach.
    name = "lk"
    nlay, nrow, ncol = 1, 11, 11
    delr = delc = 100.0
    top, botm = 0.0, [-100.0]
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
    # confined cells stay saturated; no REWET ambiguity
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    # aquifer heads above the flat lakebed (cell top = 0): both lakes stay wet
    chd = [[(0, i, 0), 3.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    # upper lake (0): small lakebed leakance, fed by a specified (area-
    # independent) inflow so its stage rises above the outlet invert and it
    # spills; lower lake (1): larger leakance so it passes the received water to
    # the aquifer. A specified INFLOW (not RAINFALL) avoids the empty-lake
    # steady state that area-scaled rainfall would also admit.
    up = [(r, c) for r in range(2, 4) for c in range(4, 7)]
    lo = [(r, c) for r in range(7, 9) for c in range(4, 7)]
    conn = []
    for k, (r, c) in enumerate(up):
        conn.append([0, k, (0, r, c), "VERTICAL", 1.0e-2, 0.0, 0.0, 0.0, 0.0])
    for k, (r, c) in enumerate(lo):
        conn.append([1, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0])
    pkd = [[0, 3.0, len(up)], [1, 2.0, len(lo)]]
    # MANNING outlet from lake 0 into lake 1 (invert above the lakebed, below the
    # upper lake's inflow-fed stage so it spills)
    outlets = [[0, 0, 1, "MANNING", 1.0, 10.0, 0.03, 1.0e-3]]
    perdata = [
        [0, "INFLOW", 50.0],
        [0, "EVAPORATION", 0.0],
        [1, "RAINFALL", 1.0e-3],
        [1, "EVAPORATION", 0.0],
    ]
    flopy.mf6.ModflowGwflak(
        gwf,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage.bin",
        nlakes=2,
        noutlets=1,
        packagedata=pkd,
        connectiondata=conn,
        outlets=outlets,
        perioddata=perdata,
        surfdep=0.5,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


def _stages(ws, name):
    sf = flopy.utils.HeadFile(os.path.join(ws, f"{name}.lak.stage.bin"), text="STAGE")
    return sf.get_data().flatten()


@pytest.mark.developmode
def test_two_lakes_outlet(function_tmpdir, targets):
    # two lakes joined by a lake-to-lake outlet (simoutrate routing). The
    # implicit formulation must converge, route the outlet flow, match the legacy
    # substitution solver on both lake stages and on the heads, and close its
    # budget. Also confirms the upper lake actually spills (stage above invert).
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, mname = _build_twolake(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_twolake(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    assert _run(sim_l), "legacy solver failed for the two-lake outlet model"
    assert _run(sim_i), "IMPLICIT failed to converge for the two-lake outlet model"

    _assert_budget_closes(str(ws_impl), mname)

    hl = _heads(str(ws_legacy), mname)
    hi = _heads(str(ws_impl), mname)
    maxdiff = float(np.nanmax(np.abs(hl - hi)))
    assert maxdiff < 1e-3, f"two-lake head mismatch: {maxdiff}"

    sl = _stages(str(ws_legacy), mname)
    si = _stages(str(ws_impl), mname)
    assert np.allclose(sl, si, atol=1e-3), f"two-lake stage mismatch: {sl} vs {si}"
    # the upper lake must spill for the outlet routing to be exercised
    assert sl[0] > 1.0, f"upper lake did not reach the outlet invert: stage {sl[0]}"


def _build_multiperiod(ws, exe):
    # a lake over two stress periods: a steady-state period followed by a
    # transient period in which the rainfall jumps, driving a transient
    # lake-stage rise. Exercises the implicit assembly across multiple stress
    # periods and time steps, and the per-time-step reset of the substitution
    # fallback flags (ifallback is cleared at the start of each time step), which
    # the single-period cases never reach.
    name = "lk"
    nlay, nrow, ncol = 1, 11, 11
    delr = delc = 100.0
    top, botm = 10.0, [-100.0]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (200.0, 10, 1.0)])
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
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-5,
        sy=0.2,
        steady_state={0: True},
        transient={1: True},
    )
    chd = [[(0, i, 0), 3.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 1.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    lr, lc = range(4, 7), range(4, 7)
    conn = [
        [0, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]
        for k, (r, c) in enumerate((r, c) for r in lr for c in lc)
    ]
    pkd = [[0, 5.0, len(conn)]]
    # rainfall jumps in the transient period, driving a transient stage rise
    perdata = {
        0: [[0, "RAINFALL", 1.0e-3], [0, "EVAPORATION", 0.0]],
        1: [[0, "RAINFALL", 2.0e-2]],
    }
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
def test_multiperiod_ss_to_transient(function_tmpdir, targets):
    # the implicit formulation must assemble and converge across multiple stress
    # periods (steady-state then transient), reset its fallback flags each time
    # step, and match the legacy substitution solver, with the budget closing.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, mname = _build_multiperiod(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_multiperiod(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    assert _run(sim_l), "legacy solver failed for the multi-period lake"
    assert _run(sim_i), "IMPLICIT failed to converge for the multi-period lake"

    _assert_budget_closes(str(ws_impl), mname)

    hl = _heads(str(ws_legacy), mname)
    hi = _heads(str(ws_impl), mname)
    maxdiff = float(np.nanmax(np.abs(hl - hi)))
    assert maxdiff < 1e-3, f"multi-period head mismatch: {maxdiff}"
    sl = _stage(str(ws_legacy), mname)
    si = _stage(str(ws_impl), mname)
    assert abs(sl - si) < 1e-3, f"multi-period stage mismatch: {sl} vs {si}"


@pytest.mark.developmode
def test_nonnewton_confined(function_tmpdir, targets):
    # the implicit formulation must assemble and converge under the standard
    # (non-Newton) formulation, not only under NEWTON. With confined cells that
    # stay saturated, the implicit result must match the legacy substitution
    # solver, and both budgets must close.
    exe = targets["mf6"]
    ws_legacy = function_tmpdir / "legacy"
    ws_impl = function_tmpdir / "implicit"
    ws_legacy.mkdir()
    ws_impl.mkdir()

    sim_l, mname = _build_confined(str(ws_legacy), exe)
    sim_l.write_simulation(silent=True)
    sim_i, _ = _build_confined(str(ws_impl), exe)
    sim_i.write_simulation(silent=True)
    _enable_implicit(str(ws_impl / f"{mname}.lak"))

    assert _run(sim_l), "legacy solver failed for the confined lake"
    assert _run(sim_i), "IMPLICIT failed to converge for the confined lake"

    _assert_budget_closes(str(ws_impl), mname)

    hl = _heads(str(ws_legacy), mname)
    hi = _heads(str(ws_impl), mname)
    maxdiff = float(np.nanmax(np.abs(hl - hi)))
    assert maxdiff < 1e-3, f"non-Newton head mismatch: {maxdiff}"
    sl = _stage(str(ws_legacy), mname)
    si = _stage(str(ws_impl), mname)
    assert abs(sl - si) < 1e-3, f"non-Newton stage mismatch: {sl} vs {si}"


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

    _assert_budget_closes(str(ws_impl), mname)
    _assert_budget_closes(str(ws_fb), mname)

    hl = _heads(str(ws_legacy), mname)
    for label, ws in (("implicit", ws_impl), ("fallback", ws_fb)):
        hx = _heads(str(ws), mname)
        maxdiff = float(np.nanmax(np.abs(hl - hx)))
        assert maxdiff < 1e-6, f"{label} head mismatch vs legacy: {maxdiff}"
        sl = _stage(str(ws_legacy), mname)
        sx = _stage(str(ws), mname)
        assert abs(sl - sx) < 1e-6, f"{label} stage mismatch: {sl} vs {sx}"
