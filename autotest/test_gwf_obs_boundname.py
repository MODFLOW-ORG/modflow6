"""
Test observations that name a boundary which doesn't exist.

In a list package, a name is an error only without BOUNDNAMES, since the list can
change between periods. The name array has no elements then, so the search read out
of bounds, and the observation wrote the no-data value with no error.

CSUB and UZF boundaries are fixed for the simulation, so any unmatched name is an
error, with or without BOUNDNAMES. CSUB wrote no-data, and UZF wrote zero.
"""

import flopy
import pytest

errmsg = "but BOUNDNAMES is not specified"
invalid = 'Boundary "NOPE" for observation "X" is invalid'


@pytest.mark.parametrize("boundnames", [True, False])
def test_obs_boundname(function_tmpdir, targets, boundnames):
    sim = flopy.mf6.MFSimulation(sim_ws=str(function_tmpdir), exe_name=targets["mf6"])
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, ncol=3)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 2), 1.0]])
    row = [(0, 0, 0), -1.0] + (["mywell"] if boundnames else [])
    flopy.mf6.ModflowGwfwel(
        gwf,
        boundnames=boundnames,
        stress_period_data=[row],
        observations={"wel.csv": [("x", "wel", "mywell")]},
    )
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    if boundnames:
        assert success, "\n".join(buff)
    else:
        assert not success
        output = " ".join(" ".join(buff).split())
        assert errmsg in output
        assert "gwf.wel" in output


def gwf_model(ws, exe):
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, complexity="complex")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", newtonoptions="newton")
    flopy.mf6.ModflowGwfdis(gwf, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 2), 5.0]])
    return sim, gwf


def add_csub(gwf, boundnames, obs):
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=0.1, transient={0: True})
    flopy.mf6.ModflowGwfcsub(
        gwf,
        boundnames=boundnames,
        ninterbeds=2,
        cg_ske_cr=1e-5,
        cg_theta=0.3,
        packagedata=[
            [i, (0, 0, i), "nodelay", 0.0, 0.1, 1.0, 1e-3, 1e-4, 0.3, 0.0, 5.0]
            + ([f"bnd{i}"] if boundnames else [])
            for i in range(2)
        ],
        observations={"obs.csv": obs},
    )


def add_uzf(gwf, boundnames, obs):
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.1, transient={0: True})
    flopy.mf6.ModflowGwfuzf(
        gwf,
        boundnames=boundnames,
        ntrailwaves=7,
        nwavesets=40,
        nuzfcells=2,
        packagedata=[
            [i, (0, 0, i), 1, -1, 0.1, 1.0, 0.05, 0.3, 0.1, 4.0]
            + ([f"bnd{i}"] if boundnames else [])
            for i in range(2)
        ],
        perioddata={0: [[i, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in range(2)]},
        observations={"obs.csv": obs},
    )


PACKAGES = {"csub": add_csub, "uzf": add_uzf}


@pytest.mark.parametrize(
    "pkg, obstype, depth",
    [
        ("csub", "csub", ()),
        ("csub", "delay-flowtop", ()),
        ("uzf", "uzf-gwrch", ()),
        ("uzf", "water-content", (0.5,)),
    ],
)
@pytest.mark.parametrize("boundnames", [True, False])
def test_obs_invalid_boundname(
    function_tmpdir, targets, pkg, obstype, depth, boundnames
):
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    PACKAGES[pkg](gwf, boundnames, [("x", obstype, "nope", *depth)])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert not success
    output = " ".join(" ".join(buff).split())
    assert invalid in output
    assert f"gwf.{pkg}" in output


@pytest.mark.parametrize(
    "pkg, obstype, depth",
    [
        ("csub", "csub", ()),
        ("uzf", "uzf-gwrch", ()),
        ("uzf", "water-content", (0.5,)),
    ],
)
def test_obs_valid_boundname(function_tmpdir, targets, pkg, obstype, depth):
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    PACKAGES[pkg](gwf, True, [("x", obstype, "bnd0", *depth)])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff)


def test_uzf_obs_invalid_node(function_tmpdir, targets):
    sim, gwf = gwf_model(function_tmpdir, targets["mf6"])
    add_uzf(gwf, False, [("x", "uzf-gwrch", 99)])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert not success
    output = " ".join(" ".join(buff).split())
    assert "Invalid node number in OBS input: 99" in output
    assert "gwf.uzf" in output
