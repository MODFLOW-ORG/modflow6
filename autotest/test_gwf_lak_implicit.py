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
from framework import TestFramework

cases = ["drylake", "drain", "fill"]


def _write_implicit(test, builder, force_fallback=False):
    # build and write the implicit simulation in the test workspace, enabling
    # the IMPLICIT option (and optionally the forced substitution fallback)
    sim, name = builder(str(test.workspace), test.targets["mf6"])
    _set_implicit(sim, name)
    sim.write_simulation(silent=True)
    lak_file = str(test.workspace / f"{name}.lak")
    _assert_implicit_written(lak_file)
    if force_fallback:
        _enable_dev_fallback(lak_file)
    return sim, name


def _write_legacy(test, builder):
    # build and write the legacy comparison simulation in the "mf6" subdirectory
    # so the framework can compare it (compare="mf6") against the implicit run
    sim, _ = builder(str(test.workspace / "mf6"), test.targets["mf6"])
    sim.write_simulation(silent=True)
    return sim


def _framework(function_tmpdir, targets, build, check, compare="mf6", xfail=False):
    # all implicit LAK tests share the same harness: the input files are written
    # by hand (to add the IMPLICIT keyword and other dev options flopy cannot
    # write), so overwrite=False keeps the framework from rewriting them
    return TestFramework(
        name="lak_implicit",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        compare=compare,
        xfail=xfail,
        overwrite=False,
    )


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


def _set_implicit(sim, modelname):
    # enable the IMPLICIT option on the model's LAK package through the flopy
    # keyword (the same way the other LAK implicit tests do)
    gwf = sim.get_model(modelname)
    lak = next(p for p in gwf.packagelist if p.package_type.lower() == "lak")
    lak.implicit = True


def _assert_implicit_written(lak_file):
    # guard against a silent no-op: if flopy predates the IMPLICIT keyword it is
    # dropped and the implicit run would quietly match the default formulation.
    # CI runs update-flopy before the tests, so this only trips a stale flopy.
    assert "IMPLICIT" in open(lak_file).read().upper(), (
        f"IMPLICIT was not written to {os.path.basename(lak_file)}; run update-flopy"
    )


def _enable_dev_fallback(lak_file):
    # add the hidden DEV_FORCE_FALLBACK option to the LAK OPTIONS block, which
    # routes every active lake through the substitution fallback path. It is a
    # development-only option that flopy does not expose, so unlike IMPLICIT it
    # is written directly to the input file.
    lines = open(lak_file).read().splitlines()
    out = []
    done = False
    for ln in lines:
        out.append(ln)
        if ln.strip().lower() == "begin options" and not done:
            out.append("  DEV_FORCE_FALLBACK")
            done = True
    assert done, f"could not find OPTIONS block in {lak_file}"
    open(lak_file, "w").write("\n".join(out) + "\n")


def _make_connectionless(lak_file):
    # rewrite the LAK input so its single lake has no groundwater connections:
    # set NLAKECONN to 0 in the PACKAGEDATA row and empty the CONNECTIONDATA
    # block. flopy will not build such a lake, so it is edited into the file.
    text = open(lak_file).read()
    # zero the NLAKECONN field (last column of the single packagedata row)
    text = re.sub(
        r"(BEGIN packagedata\s*\n\s+\S+\s+\S+\s+)\d+",
        r"\g<1>0",
        text,
    )
    # drop every connection from the CONNECTIONDATA block
    text = re.sub(
        r"BEGIN connectiondata.*?END connectiondata",
        "BEGIN connectiondata\nEND connectiondata",
        text,
        flags=re.DOTALL,
    )
    open(lak_file, "w").write(text)


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
        stage_filerecord=f"{name}.lak.stage",
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
        stage_filerecord=f"{name}.lak.stage",
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
        stage_filerecord=f"{name}.lak.stage",
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
        stage_filerecord=f"{name}.lak.stage",
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
    sf = flopy.utils.HeadFile(os.path.join(ws, f"{name}.lak.stage"), text="STAGE")
    return float(sf.get_data().flatten()[0])


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    # the implicit formulation must converge for these dry-lake cases and, where
    # the legacy solver also converges, give the same heads (compared by the
    # framework) and a closed budget. The legacy substitution solver does not
    # converge the transient "drain" case, so it has no comparison run.
    builder = _BUILDERS[name]
    legacy_converges = name != "drain"

    def build(test):
        sim_i, _ = _write_implicit(test, builder)
        if legacy_converges:
            return sim_i, _write_legacy(test, builder)
        return sim_i, None

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")

    _framework(
        function_tmpdir,
        targets,
        build,
        check,
        compare="mf6" if legacy_converges else None,
    ).run()


def test_perched_disconnected_implicit(function_tmpdir, targets):
    # a steady-state lake completely disconnected from the aquifer (every
    # connection perched). Holding the connected-cell head at the lake bottom
    # keeps the perched leakage on the lake-stage diagonal, so the IMPLICIT
    # formulation converges and matches the legacy substitution solver.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_perched)
        return sim_i, _write_legacy(test, _build_perched)

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")

    _framework(function_tmpdir, targets, build, check).run()


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
        stage_filerecord=f"{name}.lak.stage",
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
        stage_filerecord=f"{name}.lak.stage",
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
        stage_filerecord=f"{gwfname}.lak.stage",
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


def test_transport_lkt(function_tmpdir, targets):
    # an active lake coupled to lake transport (LKT). The implicit LAK must feed
    # LKT the same lake-aquifer flows as the legacy solver, so the lake and
    # aquifer concentrations must match between the two formulations. The LKT
    # concentration output (.lkt.bin) would trip the framework's built-in head
    # comparison, so both simulations are run by the framework (compare=None)
    # and the concentrations are compared here.
    gwtname = "gwt_lk"

    def build(test):
        sim_i, gwfname, _ = _build_lkt(str(test.workspace), test.targets["mf6"])
        _set_implicit(sim_i, gwfname)
        sim_i.write_simulation(silent=True)
        _assert_implicit_written(str(test.workspace / f"{gwfname}.lak"))
        sim_l, _, _ = _build_lkt(str(test.workspace / "mf6"), test.targets["mf6"])
        sim_l.write_simulation(silent=True)
        return sim_i, sim_l

    def check(test):
        ws_i = str(test.workspace)
        ws_l = str(test.workspace / "mf6")
        _assert_budget_closes(ws_i, "gwf_lk")

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

        cl, ci = lake_conc(ws_l), lake_conc(ws_i)
        assert np.allclose(cl, ci, atol=1e-4), (
            f"lake concentration mismatch: {cl} vs {ci}"
        )
        # the lake must actually exchange solute (concentration evolves from 100)
        assert cl[-1] < 99.0, f"lake concentration did not evolve: {cl[-1]}"

        al, ai = aquifer_conc(ws_l), aquifer_conc(ws_i)
        maxdiff = float(np.nanmax(np.abs(al - ai)))
        assert maxdiff < 1e-4, f"aquifer concentration mismatch: {maxdiff}"
        assert float(np.nanmax(ai)) > 1e-3, "no solute reached the aquifer"

    _framework(function_tmpdir, targets, build, check, compare=None).run()


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
        stage_filerecord=f"{name}.lak.stage",
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


def test_constant_stage_lake(function_tmpdir, targets):
    # a constant-stage lake exercises the implicit constant-stage assembly branch
    # (hold the lake row at the specified stage). The aquifer heads and the held
    # stage must match the legacy solver, and the budget must close.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_constant_stage)
        return sim_i, _write_legacy(test, _build_constant_stage)

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")
        # the held stage must be reproduced exactly
        assert abs(_stage(str(test.workspace), "lk") - 5.0) < 1e-6, (
            "stage not held constant"
        )

    _framework(function_tmpdir, targets, build, check).run()


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
        stage_filerecord=f"{name}.lak.stage",
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
    sf = flopy.utils.HeadFile(os.path.join(ws, f"{name}.lak.stage"), text="STAGE")
    return sf.get_data().flatten()


def test_two_lakes_outlet(function_tmpdir, targets):
    # two lakes joined by a lake-to-lake outlet (simoutrate routing). The
    # implicit formulation must converge, route the outlet flow, match the legacy
    # substitution solver on both lake stages and on the heads, and close its
    # budget. Also confirms the upper lake actually spills (stage above invert).
    def build(test):
        sim_i, _ = _write_implicit(test, _build_twolake)
        return sim_i, _write_legacy(test, _build_twolake)

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")
        # the upper lake must spill for the outlet routing to be exercised
        si = _stages(str(test.workspace), "lk")
        assert si[0] > 1.0, f"upper lake did not reach the outlet invert: stage {si[0]}"

    _framework(function_tmpdir, targets, build, check).run()


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
        stage_filerecord=f"{name}.lak.stage",
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


def test_multiperiod_ss_to_transient(function_tmpdir, targets):
    # the implicit formulation must assemble and converge across multiple stress
    # periods (steady-state then transient), reset its fallback flags each time
    # step, and match the legacy substitution solver, with the budget closing.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_multiperiod)
        return sim_i, _write_legacy(test, _build_multiperiod)

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")

    _framework(function_tmpdir, targets, build, check).run()


def test_nonnewton_confined(function_tmpdir, targets):
    # the implicit formulation must assemble and converge under the standard
    # (non-Newton) formulation, not only under NEWTON. With confined cells that
    # stay saturated, the implicit result must match the legacy substitution
    # solver, and both budgets must close.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_confined)
        return sim_i, _write_legacy(test, _build_confined)

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")

    _framework(function_tmpdir, targets, build, check).run()


def _build_lak_mvr(ws, exe):
    # two lakes where the upper lake's outlet discharges to the water mover
    # (lakeout = -1) instead of directly to the lower lake, and a MVR package
    # routes that outlet flow into the lower lake. Exercises the implicit
    # formulation's mover-provider accumulation (qformvr) for lake outlets, which
    # the direct lake-to-lake outlet path does not reach.
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
    up = [(r, c) for r in range(2, 4) for c in range(4, 7)]
    lo = [(r, c) for r in range(7, 9) for c in range(4, 7)]
    conn = []
    for k, (r, c) in enumerate(up):
        conn.append([0, k, (0, r, c), "VERTICAL", 1.0e-2, 0.0, 0.0, 0.0, 0.0])
    for k, (r, c) in enumerate(lo):
        conn.append([1, k, (0, r, c), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0])
    pkd = [[0, 3.0, len(up)], [1, 2.0, len(lo)]]
    # outlet from lake 0 to the mover (lakeout = -1)
    outlets = [[0, 0, -1, "MANNING", 1.0, 10.0, 0.03, 1.0e-3]]
    perdata = [
        [0, "INFLOW", 50.0],
        [0, "EVAPORATION", 0.0],
        [1, "RAINFALL", 1.0e-3],
        [1, "EVAPORATION", 0.0],
    ]
    flopy.mf6.ModflowGwflak(
        gwf,
        mover=True,
        print_stage=True,
        stage_filerecord=f"{name}.lak.stage",
        nlakes=2,
        noutlets=1,
        packagedata=pkd,
        connectiondata=conn,
        outlets=outlets,
        perioddata=perdata,
        pname="LAK-1",
        surfdep=0.5,
    )
    # route the lake-0 outlet (outlet index 0) into lake 1 via the mover
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=1,
        maxpackages=1,
        packages=[("LAK-1",)],
        perioddata=[("LAK-1", 0, "LAK-1", 1, "FACTOR", 1.0)],
        pname="MVR-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, name


def test_lake_outlet_mover(function_tmpdir, targets):
    # a lake outlet routed through the water mover (MVR). The implicit
    # formulation's mover-provider accumulation must match the legacy solver on
    # the lake stages and heads, with the budget closing. The legacy comparison
    # is made only when the legacy solver also converges, so it is done here
    # rather than through the framework's automatic comparison.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_lak_mvr)
        return sim_i, None

    def check(test):
        ws_i = str(test.workspace)
        _assert_budget_closes(ws_i, "lk")
        # the upper lake must spill into the mover (stage above the outlet invert)
        si = _stages(ws_i, "lk")
        assert si[0] > 1.0, f"upper lake did not spill to the mover: {si[0]}"

        # compare against the legacy substitution solver only if it converges
        ws_l = test.workspace / "mf6"
        ws_l.mkdir(exist_ok=True)
        sim_l, _ = _build_lak_mvr(str(ws_l), test.targets["mf6"])
        sim_l.write_simulation(silent=True)
        if _run(sim_l):
            sl = _stages(str(ws_l), "lk")
            assert np.allclose(sl, si, atol=1e-3), f"mover stage mismatch: {sl} vs {si}"
            hl = _heads(str(ws_l), "lk")
            hi = _heads(ws_i, "lk")
            assert float(np.nanmax(np.abs(hl - hi))) < 1e-3, "mover head mismatch"

    _framework(function_tmpdir, targets, build, check, compare=None).run()


def _build_withdraw(ws, exe):
    # a transient table lake emptied by a specified WITHDRAWAL (area-independent)
    # rather than by evaporation, with the aquifer held below the lake bottom so
    # the lake drains. As the lake approaches its bottom the implicit formulation
    # ramps the withdrawal toward zero (the surfdep factor), so this exercises the
    # withdrawal branch of the implicit budget reporting and its closure.
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
    chd = [[(0, i, 0), 16.0] for i in range(nrow)] + [
        [(0, i, ncol - 1), 14.0] for i in range(nrow)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    conn = [[0, 0, (0, 7, 7), "embeddedv", 10.0, 0.0, 0.0, 1.0, 0.0]]
    pkd = [[0, 30.0, 1]]
    # small rainfall in, a steady withdrawal that exceeds it -> the lake drains
    perdata = [
        [0, "RAINFALL", 1.0e-4],
        [0, "EVAPORATION", 0.0],
        [0, "WITHDRAWAL", 10.0],
    ]
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
        stage_filerecord=f"{name}.lak.stage",
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


def test_withdrawal_dries(function_tmpdir, targets):
    # a lake emptied by a WITHDRAWAL must converge under the implicit formulation
    # and close its budget -- the withdrawal term is ramped near the lake bottom
    # the same way the matrix assembles it. Only the implicit run and its budget
    # are checked (the legacy solver is not required to converge this aggressively
    # drained case).
    def build(test):
        sim_i, _ = _write_implicit(test, _build_withdraw)
        return sim_i, None

    def check(test):
        _assert_budget_closes(str(test.workspace), "lk")
        # the lake should be drawn down well below its initial stage (30) toward
        # the lake bottom (20), exercising the near-bottom withdrawal ramp
        assert _stage(str(test.workspace), "lk") < 24.0, "lake was not drawn down"

    _framework(function_tmpdir, targets, build, check, compare=None).run()


def test_requires_bicgstab(function_tmpdir, targets):
    # the IMPLICIT formulation makes the coefficient matrix asymmetric, so the
    # symmetric CG linear acceleration must be rejected with a clear error
    # rather than silently producing a wrong or non-converged result.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_confined)
        # downgrade the linear acceleration to CG (incompatible with IMPLICIT)
        # and rewrite; IMPLICIT persists because it is now a flopy keyword
        for p in sim_i.sim_package_list:
            if p.package_type.lower() == "ims":
                p.linear_acceleration = "CG"
        sim_i.write_simulation(silent=True)
        return sim_i

    def check(test):
        msg = (test.workspace / "mfsim.lst").read_text().upper()
        assert "ASYMMETRIC" in msg and "BICGSTAB" in msg, (
            "expected an asymmetric-matrix / use-BICGSTAB error"
        )

    _framework(function_tmpdir, targets, build, check, compare=None, xfail=True).run()


def test_connectionless_lake_errors(function_tmpdir, targets):
    # a lake with no groundwater connections cannot be solved as a matrix
    # unknown (its lake row would have no off-diagonal coupling and a singular
    # diagonal), so the IMPLICIT option must terminate with an informative error
    # rather than producing a meaningless stage.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_drylake)
        _make_connectionless(str(test.workspace / "lk.lak"))
        return sim_i

    def check(test):
        msg = (test.workspace / "mfsim.lst").read_text().lower()
        assert "has no connections" in msg, msg
        assert "implicit" in msg, msg

    _framework(function_tmpdir, targets, build, check, compare=None, xfail=True).run()


@pytest.mark.developmode
def test_two_lakes_outlet_fallback(function_tmpdir, targets):
    # the two-lake outlet model with every lake forced onto the substitution
    # fallback (DEV_FORCE_FALLBACK). The fallback path handles the lake-to-lake
    # outlet routing (it zeroes and recomputes simoutrate before assembling the
    # fallback lakes), so it must reproduce the legacy result exactly. The strict
    # 1e-6 agreement is checked here rather than through the framework's looser
    # head comparison.
    def build(test):
        sim_f, _ = _write_implicit(test, _build_twolake, force_fallback=True)
        return sim_f, None

    def check(test):
        ws_f = str(test.workspace)
        _assert_budget_closes(ws_f, "lk")

        ws_l = test.workspace / "mf6"
        ws_l.mkdir(exist_ok=True)
        sim_l, _ = _build_twolake(str(ws_l), test.targets["mf6"])
        sim_l.write_simulation(silent=True)
        assert _run(sim_l), "legacy solver failed for the two-lake outlet model"

        sl = _stages(str(ws_l), "lk")
        sf = _stages(ws_f, "lk")
        assert np.allclose(sl, sf, atol=1e-6), f"fallback stage mismatch: {sl} vs {sf}"
        hl = _heads(str(ws_l), "lk")
        hf = _heads(ws_f, "lk")
        assert float(np.nanmax(np.abs(hl - hf))) < 1e-6, "fallback head mismatch"

    _framework(function_tmpdir, targets, build, check, compare=None).run()


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
    # here to give a deterministic regression test of that assembly. The strict
    # 1e-6 three-way agreement is checked here rather than through the framework's
    # looser head comparison.
    def build(test):
        sim_i, _ = _write_implicit(test, _build_weak)
        return sim_i, None

    def check(test):
        ws_i = str(test.workspace)
        _assert_budget_closes(ws_i, "lk")

        # legacy substitution reference
        ws_l = test.workspace / "mf6"
        ws_l.mkdir(exist_ok=True)
        sim_l, _ = _build_weak(str(ws_l), test.targets["mf6"])
        sim_l.write_simulation(silent=True)
        assert _run(sim_l), "legacy solver failed for the weak lake"

        # IMPLICIT with every lake forced onto the substitution fallback
        ws_fb = test.workspace / "fallback"
        ws_fb.mkdir(exist_ok=True)
        sim_f, _ = _build_weak(str(ws_fb), test.targets["mf6"])
        _set_implicit(sim_f, "lk")
        sim_f.write_simulation(silent=True)
        _assert_implicit_written(str(ws_fb / "lk.lak"))
        _enable_dev_fallback(str(ws_fb / "lk.lak"))
        assert _run(sim_f), "IMPLICIT with forced fallback failed for the weak lake"
        _assert_budget_closes(str(ws_fb), "lk")

        hl = _heads(str(ws_l), "lk")
        sl = _stage(str(ws_l), "lk")
        for label, ws in (("implicit", ws_i), ("fallback", str(ws_fb))):
            hx = _heads(ws, "lk")
            maxdiff = float(np.nanmax(np.abs(hl - hx)))
            assert maxdiff < 1e-6, f"{label} head mismatch vs legacy: {maxdiff}"
            sx = _stage(ws, "lk")
            assert abs(sl - sx) < 1e-6, f"{label} stage mismatch: {sl} vs {sx}"

    _framework(function_tmpdir, targets, build, check, compare=None).run()
