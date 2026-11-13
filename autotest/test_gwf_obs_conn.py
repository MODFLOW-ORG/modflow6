"""
Test LAK and MAW observations that give a connection number (ID2).

ID2 must be between 1 and the feature's number of connections, and ID must be
a valid lake or well number. Otherwise the simulation stops with an error.

Lake 1 and well 1 each have two connections, lake 2 and well 2 have one.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

FEATURES = {"lak": "lake", "maw": "multi-aquifer well"}

# name: (package, observations, expected error or None)
cases = {}
for pkg, obstype in [
    ("lak", "lak"),
    ("lak", "conductance"),
    ("lak", "wetted-area"),
    ("maw", "maw"),
    ("maw", "conductance"),
]:
    for iconn in [0, 3]:
        cases[f"{pkg}-{obstype}-conn{iconn}"] = (
            pkg,
            [("x", obstype, 1, iconn)],
            f"was assigned connection {iconn} of {FEATURES[pkg]} 1. "
            "Connection must be >= 1 and <= 2.",
        )
for pkg, feature in FEATURES.items():
    cases[f"{pkg}-feature3"] = (
        pkg,
        [("x", pkg, 3, 1)],
        f"was assigned {feature} 3.",
    )
    cases[f"{pkg}-valid"] = (
        pkg,
        [
            ("c11", pkg, 1, 1),
            ("c12", pkg, 1, 2),
            ("c21", pkg, 2, 1),
            ("cond12", "conductance", 1, 2),
        ],
        None,
    )


def add_lak(gwf, obs):
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=2,
        packagedata=[[0, 9.0, 2], [1, 9.0, 1]],
        connectiondata=[
            [0, 0, (0, 0, 0), "horizontal", 0.1, 5.0, 10.0, 5.0, 10.0],
            [0, 1, (0, 0, 1), "horizontal", 0.1, 5.0, 10.0, 5.0, 10.0],
            [1, 0, (0, 0, 2), "horizontal", 0.1, 5.0, 10.0, 5.0, 10.0],
        ],
        perioddata=[[0, "inflow", 10.0], [1, "inflow", 5.0]],
        observations={"obs.csv": obs},
    )


def add_maw(gwf, obs):
    flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=2,
        packagedata=[[0, 0.1, 0.0, 8.0, "thiem", 2], [1, 0.1, 0.0, 8.0, "thiem", 1]],
        connectiondata=[
            [0, 0, (0, 0, 0), 10.0, 5.0, 0.0, 0.0],
            [0, 1, (1, 0, 0), 5.0, 0.0, 0.0, 0.0],
            [1, 0, (1, 0, 1), 5.0, 0.0, 0.0, 0.0],
        ],
        perioddata=[[0, "rate", -10.0], [1, "rate", -5.0]],
        observations={"obs.csv": obs},
    )


def build_models(name, test):
    pkg, obs, _ = cases[name]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, complexity="complex")
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", newtonoptions="newton", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=2, ncol=4, delr=10.0, delc=10.0, top=10.0, botm=[5.0, 0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 3), 8.0]])
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord="gwf.cbc", saverecord=[("budget", "all")]
    )
    {"lak": add_lak, "maw": add_maw}[pkg](gwf, obs)
    return sim


def check_output(name, test):
    pkg, _, error = cases[name]
    if error is not None:
        output = " ".join(" ".join(test.buffs[0]).split())
        assert error in output
        assert f"gwf.{pkg}" in output
        return

    # each connection reports its own flow
    values = np.genfromtxt(test.workspace / "obs.csv", delimiter=",", names=True)
    c11, c12, c21 = (float(values[n]) for n in ("C11", "C12", "C21"))
    assert float(values["COND12"]) > 0.0
    assert c12 != c21

    # the observations match the connection flows in the budget file
    # (MAW observations are signed from the well's side, the budget's q
    # from the aquifer's)
    bud = test.sims[0].get_model("gwf").output.budget()
    flows = bud.get_data(text=pkg.upper())[-1]["q"]
    if pkg == "maw":
        flows = -flows
    assert np.allclose([c11, c12, c21], flows)


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(name, t),
        check=lambda t: check_output(name, t),
        xfail=cases[name][2] is not None,
    )
    test.run()
