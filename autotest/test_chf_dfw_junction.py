"""
Tests for explicit channel-flow (CHF) junctions.

A junction is a DISV1D vertex shared by two or more reaches.  When the JNC
package is active (JNC6 listed in the CHF name file) each shared vertex is
promoted to an explicit junction node that carries its own stage unknown and a
pure-continuity equation, and the direct reach-reach connections at that vertex
are replaced (masked) by routing through the junction.

Cases:

  chf-jnc-mb   : steady-state Y confluence (three reaches meeting at one
                 vertex, two inflows and a zero-depth-gradient (ZDG) outlet)
                 with junctions ON; checks the junction is detected, the direct
                 reach-reach connections are masked (reach-reach FLOW-JA-FACE
                 off-diagonals are zero, so flow routes only through the
                 junction), and the model conserves mass (inflow == outflow).

  chf-jnc-off  : the same confluence model with junctions OFF (JNC6 omitted);
                 checks no junction is detected and the direct reach-reach
                 connections remain active (reach-reach FLOW-JA-FACE
                 off-diagonals are nonzero).  Together with chf-jnc-mb this is
                 the on/off contrast: the same model routes through a junction
                 when JNC6 is present and couples reaches directly when it is
                 absent.

  chf-jnc-pass : straight two-reach channel whose shared mid vertex becomes a
                 (2-reach) junction (JNC on); checks the junction is a
                 transparent pass-through (reach-reach connection masked, flow
                 in == flow out, no loss or gain at the junction).

Note on backflow: because an active junction masks the direct reach-reach
connections, reach-to-reach flow at a junction is zero by construction for any
forcing; there is no code path left for the spurious "backflow" the old
pairwise-clique coupling could produce.  The ON cases assert that masking (zero
reach-reach off-diagonals); the OFF case shows the unmasked coupling it
replaces.

TODO(jnc-test): possible future additions, not blockers for the current
feature:
  - exercise a reduced DISV1D grid (IDOMAIN excluding reaches) once the JNC
    package supports it (it currently guards against reduced grids).
  - add junction-observation checks once the JNC package emits observations.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["chf-jnc-mb", "chf-jnc-off", "chf-jnc-pass"]

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


def _sim(name, workspace):
    """Create the simulation + IMS shared by every case."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=workspace
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)], time_units="SECONDS")
    flopy.mf6.ModflowIms(sim, print_option="summary", **_IMS)
    return sim


def _add_common(chf, name, strt, with_jnc):
    """Add the packages common to every case.

    with_jnc toggles the JNC package: present enables explicit junctions,
    absent leaves the reaches coupled directly through the grid connectivity.
    """
    flopy.mf6.ModflowChfic(chf, strt=strt)
    flopy.mf6.ModflowChfdfw(chf, manningsn=0.03, idcxs=0)
    flopy.mf6.ModflowChfcxs(chf, **_RECT_CXS)
    if with_jnc:
        flopy.mf6.ModflowChfjnc(chf, save_flows=True)
    flopy.mf6.ModflowChfoc(
        chf,
        budget_filerecord=f"{name}.bud",
        stage_filerecord=f"{name}.stage",
        saverecord=[("STAGE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("STAGE", "LAST"), ("BUDGET", "ALL")],
    )


def _build_confluence(test, name, with_jnc):
    """Steady Y confluence: reaches 0,1 inflow -> shared vertex 3 -> reach 2
    outlet (ZDG).  Shared by the junctions-on (mb) and junctions-off cases;
    the only difference is whether the JNC package is added."""
    sim = _sim(name, test.workspace)
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    vertices = [
        [0, -100.0, 100.0],
        [1, 100.0, 100.0],
        [2, 0.0, -100.0],
        [3, 0.0, 0.0],
    ]
    cell1d = [
        [0, 0.5, 2, 0, 3],
        [1, 0.5, 2, 1, 3],
        [2, 0.5, 2, 3, 2],
    ]
    flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=3,
        nvert=4,
        width=10.0,
        bottom=[2.0, 1.0, 0.0],
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )
    flopy.mf6.ModflowChfsto(chf, save_flows=True, steady_state={0: True})
    _add_common(chf, name, strt=[3.0, 2.0, 0.5], with_jnc=with_jnc)
    flopy.mf6.ModflowChfflw(
        chf,
        maxbound=2,
        print_flows=True,
        stress_period_data=[(0, 3.0), (1, 1.0)],
    )
    # ZDG outlet: (cellid, idcxs, width, slope, rough)
    flopy.mf6.ModflowChfzdg(
        chf,
        maxbound=1,
        print_flows=True,
        stress_period_data=[(2, 0, 10.0, 0.001, 0.03)],
    )
    return sim, None


def build_mb(test):
    """Confluence with junctions ON."""
    return _build_confluence(test, cases[0], with_jnc=True)


def build_off(test):
    """Confluence with junctions OFF (same model, JNC6 omitted)."""
    return _build_confluence(test, cases[1], with_jnc=False)


def _build_straight(name, workspace, with_jnc):
    """Straight two-reach channel.  The reaches share mid vertex 1; with
    junctions ON that vertex becomes a 2-reach junction, with junctions OFF the
    reaches couple directly through the grid.  Both are the same physical
    straight channel, so the two should produce the same solution -- the basis
    for the transparent pass-through check."""
    sim = _sim(name, workspace)
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)

    vertices = [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 200.0, 0.0]]
    cell1d = [[0, 0.5, 2, 0, 1], [1, 0.5, 2, 1, 2]]
    flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=2,
        nvert=3,
        width=10.0,
        bottom=[1.0, 0.0],
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )
    flopy.mf6.ModflowChfsto(chf, save_flows=True, steady_state={0: True})
    _add_common(chf, name, strt=[2.0, 1.0], with_jnc=with_jnc)
    flopy.mf6.ModflowChfflw(
        chf, maxbound=1, print_flows=True, stress_period_data=[(0, 2.0)]
    )
    flopy.mf6.ModflowChfzdg(
        chf,
        maxbound=1,
        print_flows=True,
        stress_period_data=[(1, 0, 10.0, 0.01, 0.03)],
    )
    return sim


def build_pass(test):
    """Straight two-reach channel with junctions ON (mid vertex is a
    2-reach junction)."""
    return _build_straight(cases[2], test.workspace, with_jnc=True), None


def _njunctions_in(workspace, name):
    """Number of junctions reported in the listing, or 0 if none were reported.

    The JNC package emits a line like ``CHF JNC: detected 1 junction(s).`` when
    active.  When the package is inactive (JNC6 omitted) no such line is
    written, which is a legitimate zero rather than a parse failure, so this
    returns 0 in that case.
    """
    import re

    lst = (workspace / f"{name}.lst").read_text()
    match = re.search(r"CHF JNC: detected\s+(\d+)\s+junction", lst)
    return int(match.group(1)) if match else 0


def _njunctions(test, name):
    return _njunctions_in(test.workspace, name)


def _flowja_offdiag(test, name):
    """Return the reach-reach off-diagonal FLOW-JA-FACE entries.

    FLOW-JA-FACE is a flat CSR array over the reach grid connectivity
    (dis%con) only; junction connections live in the solution matrix, not
    dis%con, so they never appear here (the reach-junction exchange is folded
    onto the reach diagonals instead).  To isolate the reach-reach exchange
    from the cell residuals carried on the diagonals, this returns only the
    off-diagonal entries, located via the ia/ja structure in the binary grid
    file.
    """
    from flopy.mf6.utils.binarygrid_util import MfGrdFile

    bud = flopy.utils.binaryfile.CellBudgetFile(str(test.workspace / f"{name}.bud"))
    fja = bud.get_data(text="FLOW-JA-FACE")[-1].flatten()

    grb = MfGrdFile(str(test.workspace / f"{name}.disv1d.grb"))
    ia = grb.ia
    ja = grb.ja

    offdiag = []
    for n in range(len(ia) - 1):
        for ipos in range(ia[n], ia[n + 1]):
            if ja[ipos] != n:  # skip the diagonal
                offdiag.append(fja[ipos])
    return np.array(offdiag)


def _budget_discrepancy(test, name):
    """Max |percent discrepancy| across all budget prints in the listing."""
    import re

    lst = (test.workspace / f"{name}.lst").read_text()
    vals = [
        abs(float(x)) for x in re.findall(r"PERCENT DISCREPANCY\s*=\s*([-\d.E+]+)", lst)
    ]
    return max(vals) if vals else float("nan")


def _flw_zdg_totals(test, name):
    """Total FLW inflow and ZDG outflow from the budget file."""
    bud = flopy.utils.binaryfile.CellBudgetFile(str(test.workspace / f"{name}.bud"))
    qflw = bud.get_data(text="FLW")[-1]
    qzdg = bud.get_data(text="ZDG")[-1]
    flw_in = sum(float(r["q"]) for r in qflw if float(r["q"]) > 0.0)
    zdg_out = -sum(float(r["q"]) for r in qzdg if float(r["q"]) < 0.0)
    return flw_in, zdg_out


def check_mb(test):
    name = cases[0]
    # one junction (vertex shared by the three reaches)
    assert _njunctions(test, name) == 1, "expected exactly one junction"
    # the direct reach-reach connections at the confluence are masked, so every
    # reach-reach off-diagonal FLOW-JA-FACE entry is zero (flow routes only
    # through the junction, never reach-to-reach)
    offdiag = _flowja_offdiag(test, name)
    assert np.allclose(offdiag, 0.0, atol=1e-6), (
        f"reach-reach flow not masked at junction: {offdiag}"
    )
    # mass balance: inflow (FLW) must equal outflow (ZDG)
    disc = _budget_discrepancy(test, name)
    assert disc < 0.1, f"budget discrepancy too large: {disc}%"
    flw_in, zdg_out = _flw_zdg_totals(test, name)
    assert np.isclose(flw_in, 4.0, atol=1e-3), f"inflow {flw_in} != 4.0"
    assert np.isclose(flw_in, zdg_out, rtol=1e-3), (
        f"inflow {flw_in} != outflow {zdg_out}"
    )


def check_off(test):
    name = cases[1]
    # junctions are off, so none should be detected
    assert _njunctions(test, name) == 0, "JNC6 omitted: no junction should be detected"
    # with no junction the reaches remain directly coupled through the grid, so
    # at least one reach-reach off-diagonal flow is nonzero (the coupling that
    # an active junction masks)
    offdiag = _flowja_offdiag(test, name)
    assert np.any(np.abs(offdiag) > 1e-6), (
        f"junctions off: expected nonzero reach-reach flow, got {offdiag}"
    )


def check_pass(test):
    name = cases[2]
    # the shared mid vertex is promoted to a (2-reach) junction
    assert _njunctions(test, name) == 1, "expected one 2-reach junction"
    disc = _budget_discrepancy(test, name)
    assert disc < 0.1, f"budget discrepancy too large: {disc}%"
    # the direct reach-reach connection is masked, so every reach-reach
    # off-diagonal FLOW-JA-FACE entry is zero; all flow passes through the
    # junction (transparent pass-through)
    offdiag = _flowja_offdiag(test, name)
    assert np.allclose(offdiag, 0.0, atol=1e-6), (
        f"reach-reach flow not masked at junction: {offdiag}"
    )
    # inflow equals outflow (nothing lost/gained at the junction)
    flw_in, zdg_out = _flw_zdg_totals(test, name)
    assert np.isclose(flw_in, zdg_out, rtol=1e-3), (
        f"pass-through not conservative: in {flw_in} out {zdg_out}"
    )

    # NOTE: a mid-channel 2-reach junction is NOT asserted to exactly reproduce
    # the same straight channel solved without a junction.  It is close but not
    # identical: the reach-junction half-cell conductance uses a local
    # per-half-cell gradient, whereas DFW's reach-reach get_cond applies one
    # shared full-length gradient to both half-cells, and because conductance
    # scales as 1/sqrt(gradient) the two do not coincide.  Which convention is
    # correct is an open question (see TODO(dfw-halfcell-gradient) in
    # swf-dfw.f90); until it is resolved, do not add an exact
    # junction-vs-no-junction stage comparison here.


@pytest.mark.developmode
@pytest.mark.parametrize(
    "idx, name, builder, checker",
    [
        (0, cases[0], build_mb, check_mb),
        (1, cases[1], build_off, check_off),
        (2, cases[2], build_pass, check_pass),
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
