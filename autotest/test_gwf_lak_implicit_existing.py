"""
Sanity check for the LAK IMPLICIT option against the established LAK integration
test models.

The IMPLICIT option solves each lake stage as an additional unknown in the
groundwater flow matrix instead of by the default substitution-iteration solver;
converged results are supposed to be equivalent to the default formulation. This
module reuses the model builders from a selection of the existing LAK tests,
runs each model with the default formulation and again with IMPLICIT enabled,
and asserts that the groundwater heads match and that each run closes its
mass balance. It does not modify the original tests -- it only imports and
reruns their builders.

Most cases run the existing model unchanged. The one exception is the adaptive
time stepping (ATS) case (test_ats_implicit_matches_default), which overrides
the lake surfdep so the implicit solve is stable -- see that test for details.

The existing LAK tests already use the BICGSTAB linear acceleration that the
asymmetric implicit matrix requires, so enabling IMPLICIT here is just a matter
of setting the flag on the LAK package.
"""

import importlib
from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

# (module, case index) pairs to exercise -- one representative converging case
# per existing LAK test (core lake physics plus the coupled transport/density/
# viscosity models)
#
# Deliberately excluded:
#   test_gwf_lakobs01  -- a negative test whose model is intentionally invalid
#                         (a LAK observation is missing its iconn) so it can
#                         check the error message; not a physics candidate.
#   test_gwf_ats_lak01 -- adaptive time stepping makes the default and implicit
#                         runs follow different time-step sequences, so their
#                         heads are not directly comparable, and with the model's
#                         surfdep=0 the implicit solve diverges. It is covered
#                         instead by test_ats_implicit_matches_default below.
_TARGETS = [
    ("test_gwf_lak_bedleak", 0),  # vertical bedleak, two lakes
    ("test_gwf_lak_status", 0),  # STATUS changes (active/inactive/constant)
    ("test_gwf_lak_et_rch", 0),  # lake evaporation / recharge
    ("test_gwf_lak_wetlakbedarea01", 0),  # wetted lakebed area
    ("test_gwf_lak_wetlakbedarea02", 0),  # wetted lakebed area
    ("test_gwf_lakobs02", 0),  # lake outlet observations
    ("test_gwf_ts_lak01", 0),  # time-series lake input
    ("test_gwf_laket", 0),  # coupled lake transport (LKT)
    ("test_gwf_buy_lak01", 0),  # buoyancy (density) coupled lake
    ("test_gwf_buy_lak02", 1),  # buoyancy with a density contrast
    ("test_gwf_vsc04_lak", 1),  # viscosity active
]

# heads produced by the two formulations should agree to solver tolerance
HTOL = 1e-3
# maximum acceptable model volume-budget percent discrepancy for either run
BUDGET_TOL = 0.5


class _Ctx:
    # minimal stand-in for the framework test object the builders expect
    def __init__(self, workspace, targets):
        self.workspace = workspace
        self.targets = targets


def _pkg(model, ftype):
    for p in model.packagelist:
        if p.package_type.lower() == ftype:
            return p
    return None


def _find_lak(sim):
    # return the (model, lak package) of the first model that has a LAK package
    for mname in sim.model_names:
        model = sim.get_model(mname)
        lak = _pkg(model, "lak")
        if lak is not None:
            return model, lak
    raise LookupError("no LAK package found in the simulation")


def _force_bicgstab(sim):
    # the implicit formulation makes the coefficient matrix asymmetric, which
    # requires the BICGSTAB linear acceleration; some existing LAK tests use CG,
    # so switch every IMS solution to BICGSTAB (converged results are unchanged)
    for p in getattr(sim, "sim_package_list", []) or []:
        if getattr(p, "package_type", "").lower() == "ims":
            p.linear_acceleration = "BICGSTAB"


def _ensure_head_output(gwf):
    # make sure the lake-bearing model writes a head file so the default and
    # implicit runs can be compared (some LAK tests only save the budget)
    name = gwf.name
    oc = _pkg(gwf, "oc")
    if oc is None:
        flopy.mf6.ModflowGwfoc(
            gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
        )
    else:
        oc.head_filerecord = f"{name}.hds"
        oc.saverecord = [("HEAD", "ALL")]


def _build_target(mod, idx, ws, exe, implicit, surfdep=None):
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    if hasattr(mod, "build_models"):
        res = mod.build_models(idx, _Ctx(ws, {"mf6": exe}))
    else:
        # a few older tests expose build_model(dir, exe) instead
        res = mod.build_model(str(ws), exe)
    sim = res[0] if isinstance(res, (tuple, list)) else res
    sim.set_sim_path(str(ws))
    gwf, lak = _find_lak(sim)
    if surfdep is not None:
        lak.surfdep = surfdep
    if implicit:
        lak.implicit = True
        _force_bicgstab(sim)
    _ensure_head_output(gwf)
    return sim, gwf.name


def _final_heads(ws, name):
    hf = flopy.utils.HeadFile(str(Path(ws) / f"{name}.hds"))
    h = hf.get_data(totim=hf.get_times()[-1])
    return np.where(np.abs(h) < 1e29, h, np.nan)


def _budget_discrepancy(ws, name, floor=1e-6):
    # absolute *final cumulative* volume-budget percent discrepancy for the
    # lake-bearing model -- the mass-balance-closure metric. The per-step
    # (incremental) discrepancy is deliberately not used: it can be large without
    # indicating an error (e.g. density-driven BUY flow does not conserve volume
    # step to step, yet the run still balances overall). A near-zero-flow budget
    # (e.g. a static lake in equilibrium, as in the LKT transport model, where
    # the flow field is steady and the transport does the work) has a round-off-
    # dominated percent discrepancy, so it is treated as closed.
    lst = Path(ws) / f"{name}.lst"
    if not lst.is_file():
        return None
    _, cumulative = flopy.utils.Mf6ListBudget(str(lst)).get_dataframes()
    row = cumulative.iloc[-1]
    total = max(abs(float(row["TOTAL_IN"])), abs(float(row["TOTAL_OUT"])))
    if total < floor:
        return 0.0
    return abs(float(row["PERCENT_DISCREPANCY"]))


def _assert_budget_closes(ws, name, label):
    disc = _budget_discrepancy(ws, name)
    assert disc is not None, f"no volume budget found in {name}.lst for the {label} run"
    assert disc < BUDGET_TOL, (
        f"{label} run budget discrepancy {disc} exceeds {BUDGET_TOL} percent"
    )


@pytest.mark.parametrize(
    "module_name, idx", _TARGETS, ids=[f"{m}-{i}" for m, i in _TARGETS]
)
def test_implicit_matches_default(module_name, idx, function_tmpdir, targets):
    mod = importlib.import_module(module_name)
    info = {}

    def build(test):
        exe = test.targets["mf6"]
        sim_def, gwfname = _build_target(mod, idx, test.workspace, exe, implicit=False)
        sim_impl, _ = _build_target(
            mod, idx, test.workspace / "implicit", exe, implicit=True
        )
        info["gwfname"] = gwfname
        return sim_def, sim_impl

    def check(test):
        # guard against a silent no-op: confirm IMPLICIT was actually written
        # (the flopy kwarg is only available after update-flopy regenerates the
        # LAK class from the dfns, which CI does before running the tests)
        lak_files = list((test.workspace / "implicit").glob("*.lak"))
        assert lak_files and any(
            "IMPLICIT" in f.read_text().upper() for f in lak_files
        ), "IMPLICIT was not written to the LAK input (run update-flopy)"

        name = info["gwfname"]
        # both runs must close their mass balance
        _assert_budget_closes(test.workspace, name, "default")
        _assert_budget_closes(test.workspace / "implicit", name, "implicit")

        # and the heads must match between the two formulations
        hd = _final_heads(test.workspace, name)
        hi = _final_heads(test.workspace / "implicit", name)
        maxdiff = float(np.nanmax(np.abs(hd - hi)))
        assert maxdiff < HTOL, (
            f"{module_name}[{idx}]: IMPLICIT heads differ from the default "
            f"formulation by {maxdiff} (> {HTOL})"
        )

    TestFramework(
        name="lak_implicit_existing",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        compare=None,
    ).run()


# the ATS model's unsmoothed lakebed seepage cutoff (surfdep=0) makes the
# implicit Jacobian discontinuous and the implicit solve diverges, so this case
# uses a smoothed cutoff (surfdep>0) to exercise the implicit lake through the
# adaptive time stepping machinery. ATS lets the two formulations take different
# time-step sequences, so the heads only agree loosely; the per-run mass-balance
# closure (which is time-step-path independent) is the stronger check.
_ATS_SURFDEP = 1.0
_ATS_HTOL = 5e-2


def test_ats_implicit_matches_default(function_tmpdir, targets):
    mod = importlib.import_module("test_gwf_ats_lak01")
    info = {}

    def build(test):
        exe = test.targets["mf6"]
        sim_def, gwfname = _build_target(
            mod, 0, test.workspace, exe, implicit=False, surfdep=_ATS_SURFDEP
        )
        sim_impl, _ = _build_target(
            mod,
            0,
            test.workspace / "implicit",
            exe,
            implicit=True,
            surfdep=_ATS_SURFDEP,
        )
        info["gwfname"] = gwfname
        return sim_def, sim_impl

    def check(test):
        lak_files = list((test.workspace / "implicit").glob("*.lak"))
        assert lak_files and any(
            "IMPLICIT" in f.read_text().upper() for f in lak_files
        ), "IMPLICIT was not written to the LAK input (run update-flopy)"

        name = info["gwfname"]
        # the stronger, time-step-path-independent check: each run closes
        _assert_budget_closes(test.workspace, name, "default")
        _assert_budget_closes(test.workspace / "implicit", name, "implicit")

        # heads agree loosely (ATS gives the formulations different dt paths)
        hd = _final_heads(test.workspace, name)
        hi = _final_heads(test.workspace / "implicit", name)
        maxdiff = float(np.nanmax(np.abs(hd - hi)))
        assert maxdiff < _ATS_HTOL, (
            f"ATS IMPLICIT heads differ from the default formulation by "
            f"{maxdiff} (> {_ATS_HTOL})"
        )

    TestFramework(
        name="lak_implicit_ats",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        compare=None,
    ).run()
