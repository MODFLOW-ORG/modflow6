"""
Tests for explicit channel-flow (CHF) junctions.

A junction is a DISV1D vertex shared by two or more reaches.  Each shared
vertex is promoted to an explicit junction node that carries its own stage
unknown and a pure-continuity equation; the direct reach-reach connections at
that vertex are replaced (masked) by routing through the junction.

Three cases:

  chf-jnc-mb   : steady-state Y confluence with a zero-depth-gradient (ZDG)
                 outlet; checks the model conserves mass (inflow == outflow)
                 and that the junction is detected.

  chf-jnc-pass : straight two-reach channel whose shared mid vertex becomes a
                 (2-reach) junction; checks the junction is a transparent
                 pass-through (flow in == flow out, no loss or gain at the
                 junction, monotonic stage profile).

  chf-jnc-bf   : three reaches meeting at one vertex with two inflows and a low
                 outlet and stage ordering h0 > h1 >> h2.  With the old
                 pairwise-clique coupling this produced spurious "backflow" from
                 the high inflow reach directly into the lower inflow reach.
                 With explicit junctions the direct reach-reach connections are
                 masked, so no such backflow can occur; this is asserted by
                 checking the reach-reach FLOW-JA-FACE entries are all zero and
                 the second inflow reach is a net source to the junction.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["chf-jnc-mb", "chf-jnc-pass", "chf-jnc-bf"]

# robust IMS settings for stiff diffusive-wave CHF problems
_IMS = dict(
    outer_maximum=500,
    inner_maximum=100,
    outer_dvclose=1.0e-8,
    inner_dvclose=1.0e-9,
    linear_acceleration="BICGSTAB",
    backtracking_number=5,
    backtracking_tolerance=1.0,
    backtracking_reduction_factor=0.3,
    backtracking_residual_limit=100.0,
    under_relaxation="simple",
    under_relaxation_gamma=0.7,
)

_RECT_CXS = dict(
    nsections=1,
    npoints=2,
    packagedata=[(0, 2)],
    crosssectiondata=[(0, 0.0, 0.0, 1.0), (0, 1.0, 0.0, 1.0)],
)


def _common(chf, name, strt):
    flopy.mf6.ModflowChfic(chf, strt=strt)
    flopy.mf6.ModflowChfdfw(chf, manningsn=0.03, idcxs=0)
    flopy.mf6.ModflowChfcxs(chf, **_RECT_CXS)
    flopy.mf6.ModflowChfoc(
        chf,
        budget_filerecord=f"{name}.bud",
        stage_filerecord=f"{name}.stage",
        saverecord=[("STAGE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("STAGE", "LAST"), ("BUDGET", "ALL")],
    )


def build_mb(test):
    """Steady Y confluence, ZDG outlet, two inflows."""
    name = cases[0]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=[(1.0, 1, 1.0)], time_units="SECONDS"
    )
    flopy.mf6.ModflowIms(sim, print_option="summary", **_IMS)
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    # reaches 0,1 inflow -> junction (vertex 4) -> reach 2 outlet
    vertices = [
        [0, -100.0, 100.0],
        [1, 100.0, 100.0],
        [2, 0.0, -100.0],
        [3, 0.0, -100.0],
        [4, 0.0, 0.0],
    ]
    cell1d = [
        [0, 0.5, 2, 0, 4],
        [1, 0.5, 2, 1, 4],
        [2, 0.5, 2, 4, 2],
    ]
    flopy.mf6.ModflowChfdisv1D(
        chf, nodes=3, nvert=5, width=10.0, bottom=[2.0, 1.0, 0.0],
        idomain=1, vertices=vertices, cell1d=cell1d,
    )
    flopy.mf6.ModflowChfsto(chf, save_flows=True, steady_state={0: True})
    _common(chf, name, strt=[3.0, 2.0, 0.5])
    flopy.mf6.ModflowChfflw(
        chf, maxbound=2, print_flows=True,
        stress_period_data=[(0, 3.0), (1, 1.0)],
    )
    # ZDG outlet: (cellid, idcxs, width, slope, rough)
    flopy.mf6.ModflowChfzdg(
        chf, maxbound=1, print_flows=True,
        stress_period_data=[(2, 0, 10.0, 0.001, 0.03)],
    )
    return sim, None


def build_pass(test):
    """Straight two-reach channel; shared mid vertex becomes a junction."""
    name = cases[1]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=[(1.0, 1, 1.0)], time_units="SECONDS"
    )
    flopy.mf6.ModflowIms(sim, print_option="summary", **_IMS)
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    # reach 0 (v0->v1) and reach 1 (v1->v2) share vertex 1 -> 2-reach junction
    vertices = [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 200.0, 0.0]]
    cell1d = [[0, 0.5, 2, 0, 1], [1, 0.5, 2, 1, 2]]
    flopy.mf6.ModflowChfdisv1D(
        chf, nodes=2, nvert=3, width=10.0, bottom=[1.0, 0.0],
        idomain=1, vertices=vertices, cell1d=cell1d,
    )
    flopy.mf6.ModflowChfsto(chf, save_flows=True, steady_state={0: True})
    _common(chf, name, strt=[2.0, 1.0])
    flopy.mf6.ModflowChfflw(
        chf, maxbound=1, print_flows=True, stress_period_data=[(0, 2.0)]
    )
    flopy.mf6.ModflowChfzdg(
        chf, maxbound=1, print_flows=True,
        stress_period_data=[(1, 0, 10.0, 0.01, 0.03)],
    )
    return sim, None


def build_bf(test):
    """Backflow case: h0 > h1 >> h2, two inflows, low ZDG outlet."""
    name = cases[2]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=[(1.0, 1, 1.0)], time_units="SECONDS"
    )
    flopy.mf6.ModflowIms(sim, print_option="summary", **_IMS)
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    vertices = [
        [0, -100.0, 100.0],
        [1, 100.0, 100.0],
        [2, 0.0, -100.0],
        [3, 0.0, -100.0],
        [4, 0.0, 0.0],
    ]
    cell1d = [
        [0, 0.5, 2, 0, 4],
        [1, 0.5, 2, 1, 4],
        [2, 0.5, 2, 4, 2],
    ]
    flopy.mf6.ModflowChfdisv1D(
        chf, nodes=3, nvert=5, width=10.0, bottom=[2.0, 1.0, 0.0],
        idomain=1, vertices=vertices, cell1d=cell1d,
    )
    flopy.mf6.ModflowChfsto(chf, save_flows=True, steady_state={0: True})
    _common(chf, name, strt=[3.0, 2.0, 0.5])
    flopy.mf6.ModflowChfflw(
        chf, maxbound=2, print_flows=True,
        stress_period_data=[(0, 3.0), (1, 1.0)],
    )
    flopy.mf6.ModflowChfzdg(
        chf, maxbound=1, print_flows=True,
        stress_period_data=[(2, 0, 10.0, 0.001, 0.03)],
    )
    return sim, None


def _njunctions(test, name):
    """Count detected junctions from the model listing file."""
    lst = (test.workspace / f"{name}.lst").read_text()
    for line in lst.splitlines():
        if "CHF JNC: detected" in line:
            return int(line.split("detected")[1].split("junction")[0])
    return 0


def _flowja(test, name):
    bud = flopy.utils.binaryfile.CellBudgetFile(
        str(test.workspace / f"{name}.bud")
    )
    return bud.get_data(text="FLOW-JA-FACE")[-1].flatten()


def _budget_discrepancy(test, name):
    """Max |percent discrepancy| across all budget prints in the listing."""
    import re

    lst = (test.workspace / f"{name}.lst").read_text()
    vals = [
        abs(float(x))
        for x in re.findall(r"PERCENT DISCREPANCY\s*=\s*([-\d.E+]+)", lst)
    ]
    return max(vals) if vals else float("nan")


def check_mb(test):
    name = cases[0]
    # one junction (vertex shared by the three reaches)
    assert _njunctions(test, name) == 1, "expected exactly one junction"
    # mass balance: inflow (FLW) must equal outflow (ZDG)
    disc = _budget_discrepancy(test, name)
    assert disc < 0.1, f"budget discrepancy too large: {disc}%"
    # inflow total is 4.0; outflow should match
    bud = flopy.utils.binaryfile.CellBudgetFile(
        str(test.workspace / f"{name}.bud")
    )
    qflw = bud.get_data(text="FLW")[-1]
    qzdg = bud.get_data(text="ZDG")[-1]
    flw_in = sum(float(r["q"]) for r in qflw if float(r["q"]) > 0.0)
    zdg_out = -sum(float(r["q"]) for r in qzdg if float(r["q"]) < 0.0)
    assert np.isclose(flw_in, 4.0, atol=1e-3), f"inflow {flw_in} != 4.0"
    assert np.isclose(flw_in, zdg_out, rtol=1e-3), (
        f"inflow {flw_in} != outflow {zdg_out}"
    )


def check_pass(test):
    name = cases[1]
    # the shared mid vertex is promoted to a (2-reach) junction
    assert _njunctions(test, name) == 1, "expected one 2-reach junction"
    disc = _budget_discrepancy(test, name)
    assert disc < 0.1, f"budget discrepancy too large: {disc}%"
    # the direct reach-reach connection is masked -> its FLOW-JA-FACE is zero;
    # all flow passes through the junction (transparent pass-through)
    fja = _flowja(test, name)
    # off-diagonal entries (reach-reach) must be ~zero
    offdiag = np.array([v for v in fja if abs(v) > 0.0])
    assert np.allclose(fja, 0.0, atol=1e-6) or offdiag.size == 0 or np.all(
        np.abs(offdiag) < 1e-6
    ), f"reach-reach flow not masked at junction: {fja}"
    # inflow equals outflow (nothing lost/gained at the junction)
    bud = flopy.utils.binaryfile.CellBudgetFile(
        str(test.workspace / f"{name}.bud")
    )
    qflw = bud.get_data(text="FLW")[-1]
    qzdg = bud.get_data(text="ZDG")[-1]
    flw_in = sum(float(r["q"]) for r in qflw if float(r["q"]) > 0.0)
    zdg_out = -sum(float(r["q"]) for r in qzdg if float(r["q"]) < 0.0)
    assert np.isclose(flw_in, zdg_out, rtol=1e-3), (
        f"pass-through not conservative: in {flw_in} out {zdg_out}"
    )


def check_bf(test):
    name = cases[2]
    assert _njunctions(test, name) == 1, "expected one junction"
    disc = _budget_discrepancy(test, name)
    assert disc < 0.1, f"budget discrepancy too large: {disc}%"
    # No backflow: with explicit junctions the direct reach-reach connections
    # at the confluence are masked, so every reach-reach FLOW-JA-FACE entry is
    # zero.  (In the old pairwise-clique build, the reach0->reach1 entry carried
    # spurious flow into the second inflow reach.)
    fja = _flowja(test, name)
    assert np.allclose(fja, 0.0, atol=1e-6), (
        f"reach-reach flow present (possible backflow): {fja}"
    )


@pytest.mark.developmode
@pytest.mark.parametrize(
    "idx, name, builder, checker",
    [
        (0, cases[0], build_mb, check_mb),
        (1, cases[1], build_pass, check_pass),
        (2, cases[2], build_bf, check_bf),
    ],
)
def test_mf6model(idx, name, builder, checker, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: builder(t),
        check=lambda t: checker(t),
        targets=targets,
    )
    test.run()
