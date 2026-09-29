"""
Tests for the LAK IMPLICIT option that use its development options, which
force lakes onto the substitution fallback and are only available in a
development-mode build.

Cases:
  - two_lakes_outlet_fallback   : two lakes joined by an outlet, both forced onto
                                  the fallback.
  - outlet_mover_fallback       : an outlet to the mover with both lakes forced
                                  onto the fallback; the mover must receive the
                                  legacy solver's outlet flow once.
  - outlet_mover_mixed_fallback : the same model with only the downstream lake
                                  forced, so the outlet's implicit lake keeps its
                                  rate through the fallback-only solve.
  - force_fallback_lake_errors  : DEV_FORCE_FALLBACK_LAKE with DEV_FORCE_FALLBACK,
                                  or naming a lake that does not exist, is an
                                  input error (xfail).
  - fallback_matches_legacy     : a weakly connected lake solved by the legacy
                                  solver, IMPLICIT, and IMPLICIT with the forced
                                  fallback must agree.
"""

import numpy as np
import pytest
from test_gwf_lak_implicit import (
    _assert_budget_closes,
    _assert_implicit_written,
    _build_twolake,
    _build_twolake_mvr,
    _build_weak,
    _framework,
    _heads,
    _lak_flows,
    _run,
    _set_implicit,
    _stage,
    _stages,
    _write_implicit,
)


def _enable_dev_fallback(lak_file, lake=None):
    # add the hidden DEV_FORCE_FALLBACK option to the LAK OPTIONS block, which
    # routes every active lake through the substitution fallback path, or
    # DEV_FORCE_FALLBACK_LAKE to route only the given (one-based) lake. They are
    # development-only options that flopy does not expose, so unlike IMPLICIT
    # they are written directly to the input file.
    if lake is None:
        option = "DEV_FORCE_FALLBACK"
    else:
        option = f"DEV_FORCE_FALLBACK_LAKE {lake}"
    lines = open(lak_file).read().splitlines()
    out = []
    done = False
    for ln in lines:
        out.append(ln)
        if ln.strip().lower() == "begin options" and not done:
            out.append(f"  {option}")
            done = True
    assert done, f"could not find OPTIONS block in {lak_file}"
    open(lak_file, "w").write("\n".join(out) + "\n")


def _write_dev(test, builder, lake=None):
    # write the implicit simulation with every lake, or only the given
    # (one-based) lake, forced onto the substitution fallback
    sim, name = _write_implicit(test, builder)
    _enable_dev_fallback(str(test.workspace / f"{name}.lak"), lake=lake)
    return sim, name


@pytest.mark.developmode
def test_two_lakes_outlet_fallback(function_tmpdir, targets):
    # the two-lake outlet model with every lake forced onto the substitution
    # fallback (DEV_FORCE_FALLBACK). The fallback path handles the lake-to-lake
    # outlet routing (it zeroes and recomputes simoutrate before assembling the
    # fallback lakes), so it must reproduce the legacy result exactly. The strict
    # 1e-6 agreement is checked here rather than through the framework's looser
    # head comparison.
    def build(test):
        sim_f, _ = _write_dev(test, _build_twolake)
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
def test_outlet_mover_fallback(function_tmpdir, targets):
    # regression test for #2968. test_lake_outlet_mover already covers an outlet
    # routed through the mover, but not with a lake on the substitution fallback,
    # which is what triggered the defect: the mover provider term was accumulated
    # once in lak_fc_implicit and again at the end of the fallback-only
    # lak_solve, so any model with a fallback lake advertised twice the outlet
    # flow to the mover. TO-MVR came out roughly doubled and the surplus appeared
    # as a positive EXT-OUTFLOW, an apparent inflow that let the lake budget
    # close on the wrong flows.
    def build(test):
        sim_f, _ = _write_dev(test, _build_twolake_mvr)
        return sim_f, None

    _framework(function_tmpdir, targets, build, _check_outlet_mover, compare=None).run()


@pytest.mark.developmode
def test_outlet_mover_mixed_fallback(function_tmpdir, targets):
    # the outlet mover model with only the downstream lake on the substitution
    # fallback. The fallback-only lak_solve cleared the outlet rates of every
    # lake, so the implicit upstream lake then provided no flow to the mover.
    def build(test):
        sim_f, _ = _write_dev(test, _build_twolake_mvr, lake=2)
        return sim_f, None

    _framework(function_tmpdir, targets, build, _check_outlet_mover, compare=None).run()


@pytest.mark.developmode
@pytest.mark.parametrize(
    "lake, message",
    [
        (2, "cannot both be specified"),
        (3, "must be between 1 and nlakes"),
    ],
)
def test_force_fallback_lake_errors(function_tmpdir, targets, lake, message):
    # DEV_FORCE_FALLBACK_LAKE cannot be combined with DEV_FORCE_FALLBACK and must
    # name a lake in the package; either is an input error
    def build(test):
        sim, name = _write_dev(test, _build_twolake, lake=lake)
        if message == "cannot both be specified":
            _enable_dev_fallback(str(test.workspace / f"{name}.lak"))
        return sim, None

    def check(test):
        text = open(test.workspace / "mfsim.lst").read().lower()
        assert message in text, f"expected error not reported: {message}"

    _framework(function_tmpdir, targets, build, check, compare=None, xfail=True).run()


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


def _check_outlet_mover(test):
    # the implicit run must close its lake budget and send the legacy solver's
    # outlet flow to the mover, with no spurious external outflow
    ws_f = str(test.workspace)
    _assert_budget_closes(ws_f, "lk")

    ws_l = test.workspace / "mf6"
    ws_l.mkdir(exist_ok=True)
    sim_l, _ = _build_twolake_mvr(str(ws_l), test.targets["mf6"])
    sim_l.write_simulation(silent=True)
    assert _run(sim_l), "legacy solver failed for the two-lake mover model"

    ql = _lak_flows(str(ws_l), "lk", "TO-MVR")
    qf = _lak_flows(ws_f, "lk", "TO-MVR")
    assert np.abs(ql).max() > 1.0e-6, "the outlet did not spill to the mover"
    assert np.allclose(ql, qf, rtol=1e-6, atol=1e-9), (
        f"mover flow does not match the legacy solver: {ql} vs {qf}"
    )

    # a doubled provider term shows up as a positive (inflow) EXT-OUTFLOW
    eo = _lak_flows(ws_f, "lk", "EXT-OUTFLOW")
    assert np.abs(eo).max() < 1e-6 * np.abs(ql).max(), (
        f"unexpected external outflow with a mover on the outlet: {eo}"
    )
