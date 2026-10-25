"""
Test that a package listed more than once in the PACKAGES block of a
model mover or an exchange mover is reported as an input error, and
that packages with the same name in different models are accepted.
"""

import flopy
import pytest
from framework import TestFramework
from test_gwf_exgmvr01 import build_simulation
from test_gwf_mvr01 import get_model

errmsg = "listed more than once in PACKAGES block."


def build_model_mvr(ws):
    sim = get_model(0, ws)
    return sim, sim.get_model().get_package("mvr")


def build_exchange_mvr(ws):
    sim = build_simulation(0, ws, sim_type="split")
    exg = next(iter(sim.exchange_files))
    # add a package with the same name to both models: these are
    # different packages, so the exchange mover must accept both
    packages = [tuple(row) for row in exg.mvr.packages.get_data()]
    for name in ("left", "right"):
        gwf = sim.get_model(name)
        flopy.mf6.ModflowGwfdrn(
            gwf,
            mover=True,
            pname="drn",
            stress_period_data=[[(0, 0, 0), 1.0, 1.0]],
        )
        packages.append((name, "drn"))
    exg.mvr.packages.set_data(packages)
    exg.mvr.maxpackages.set_data(len(packages))
    return sim, exg.mvr


def build_models(test, mover_type, duplicate):
    build = build_model_mvr if mover_type == "model" else build_exchange_mvr
    sim, mvr = build(test.workspace)
    if duplicate:
        # repeat the first package, differing only in case
        packages = [tuple(row) for row in mvr.packages.get_data()]
        packages.append(tuple(v.upper() for v in packages[0]))
        mvr.packages.set_data(packages)
        mvr.maxpackages.set_data(len(packages))
    return sim, None


def check_output(test, duplicate):
    buff = "\n".join(test.buffs[0])
    assert (errmsg in buff) == duplicate


@pytest.mark.parametrize("duplicate", [False, True], ids=["unique", "duplicate"])
@pytest.mark.parametrize("mover_type", ["model", "exchange"], ids=lambda t: f"{t}_mvr")
def test_mvr_duplicate_packages(function_tmpdir, targets, mover_type, duplicate):
    test = TestFramework(
        name=f"{mover_type}-mvr",
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(t, mover_type, duplicate),
        check=lambda t: check_output(t, duplicate),
        compare=None,
        xfail=duplicate,
    )
    test.run()
