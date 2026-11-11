"""
Test LAK and MAW observations that give a connection number (ID2).

ID2 must be between 1 and the feature's number of connections, and ID must be
a valid lake or well number. Otherwise the simulation stops with an error.

Lake 1 and well 1 each have two connections, lake 2 and well 2 have one.
"""

import flopy
import numpy as np
import pytest


def gwf_model(ws, exe):
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name=exe)
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
    return sim, gwf


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


PACKAGES = {
    "lak": (add_lak, "lake"),
    "maw": (add_maw, "multi-aquifer well"),
}


@pytest.mark.parametrize(
    "pkg, obstype",
    [
        ("lak", "lak"),
        ("lak", "conductance"),
        ("lak", "wetted-area"),
        ("maw", "maw"),
        ("maw", "conductance"),
    ],
)
@pytest.mark.parametrize("iconn", [0, 3])
def test_obs_invalid_conn(function_tmpdir, targets, pkg, obstype, iconn):
    add, feature = PACKAGES[pkg]
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    add(gwf, [("x", obstype, 1, iconn)])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert not success
    output = " ".join(" ".join(buff).split())
    assert (
        f"was assigned connection {iconn} of {feature} 1. "
        "Connection must be >= 1 and <= 2." in output
    )
    assert f"gwf.{pkg}" in output


@pytest.mark.parametrize("pkg", ["lak", "maw"])
def test_obs_invalid_feature(function_tmpdir, targets, pkg):
    add, feature = PACKAGES[pkg]
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    add(gwf, [("x", pkg, 3, 1)])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert not success
    output = " ".join(" ".join(buff).split())
    assert f"was assigned {feature} 3." in output
    assert f"gwf.{pkg}" in output


@pytest.mark.parametrize("pkg", ["lak", "maw"])
def test_obs_valid_conn(function_tmpdir, targets, pkg):
    add, _ = PACKAGES[pkg]
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    obs = [
        ("c11", pkg, 1, 1),
        ("c12", pkg, 1, 2),
        ("c21", pkg, 2, 1),
        ("cond12", "conductance", 1, 2),
    ]
    add(gwf, obs)
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff)

    # each connection reports its own flow
    values = np.genfromtxt(function_tmpdir / "obs.csv", delimiter=",", names=True)
    c11, c12, c21 = (float(values[n]) for n in ("C11", "C12", "C21"))
    assert float(values["COND12"]) > 0.0
    assert c12 != c21

    # the observations match the connection flows in the budget file
    # (MAW observations are signed from the well's side, the budget's q
    # from the aquifer's)
    bud = gwf.output.budget()
    flows = bud.get_data(text=pkg.upper())[-1]["q"]
    if pkg == "maw":
        flows = -flows
    assert np.allclose([c11, c12, c21], flows)
