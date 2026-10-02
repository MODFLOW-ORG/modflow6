"""
Toggle test for the explicit channel-flow (CHF) junction (JNC) package.

The JNC package is enabled by listing JNC6 in the CHF name file (presence is
the switch): present means junctions are detected and activated; absent means
no junction logic runs.  This test builds one confluence model and runs it both
ways against the same executable, asserting:

  junctions ON  (JNC6 present) : the shared confluence vertex is detected as a
                                 junction and the direct reach-reach FLOW-JA-FACE
                                 entries are masked to zero (flow routes through
                                 the junction node).

  junctions OFF (JNC6 absent)  : no junction is detected and the direct
                                 reach-reach connections remain active, so the
                                 reach-reach FLOW-JA-FACE entries are nonzero
                                 (the pre-junction pairwise coupling).

The off/on contrast is exactly what the old file-sentinel hack demonstrated,
now driven through the real MODFLOW 6 input path.
"""

from pathlib import Path

import flopy
import numpy as np

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


def _build(ws, name, exe, with_jnc):
    """Y confluence: reaches 0,1 inflow -> shared vertex 4 -> reach 2 outlet."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=exe, sim_ws=str(ws)
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
    flopy.mf6.ModflowChfic(chf, strt=[3.0, 2.0, 0.5])
    flopy.mf6.ModflowChfdfw(chf, manningsn=0.03, idcxs=0)
    flopy.mf6.ModflowChfcxs(chf, **_RECT_CXS)
    flopy.mf6.ModflowChfoc(
        chf,
        budget_filerecord=f"{name}.bud",
        stage_filerecord=f"{name}.stage",
        saverecord=[("STAGE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("STAGE", "LAST"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowChfflw(
        chf, maxbound=2, print_flows=True,
        stress_period_data=[(0, 3.0), (1, 1.0)],
    )
    flopy.mf6.ModflowChfzdg(
        chf, maxbound=1, print_flows=True,
        stress_period_data=[(2, 0, 10.0, 0.001, 0.03)],
    )

    # the toggle: JNC6 present enables junctions, absent disables them
    if with_jnc:
        flopy.mf6.ModflowChfjnc(chf, save_flows=True)

    return sim


def _njunctions(ws, name):
    lst = (ws / f"{name}.lst").read_text()
    for line in lst.splitlines():
        if "CHF JNC: detected" in line:
            return int(line.split("detected")[1].split("junction")[0])
    return 0


def _flowja(ws, name):
    bud = flopy.utils.binaryfile.CellBudgetFile(str(ws / f"{name}.bud"))
    return bud.get_data(text="FLOW-JA-FACE")[-1].flatten()


def _run(sim):
    success, _ = sim.run_simulation(silent=True)
    assert success, "mf6 run did not converge"


def test_jnc_toggle(function_tmpdir):
    exe = str(
        (Path(__file__).parents[1] / "bin" / "mf6d.exe").resolve()
    )
    if not Path(exe).exists():
        exe = str((Path(__file__).parents[1] / "bin" / "mf6.exe").resolve())

    # junctions ON
    ws_on = function_tmpdir / "on"
    ws_on.mkdir()
    sim = _build(ws_on, "jncon", exe, with_jnc=True)
    sim.write_simulation(silent=True)
    _run(sim)
    assert _njunctions(ws_on, "jncon") == 1, "JNC6 present: expected 1 junction"
    fja_on = _flowja(ws_on, "jncon")
    assert np.allclose(fja_on, 0.0, atol=1e-6), (
        f"JNC6 present: reach-reach flow should be masked, got {fja_on}"
    )

    # junctions OFF
    ws_off = function_tmpdir / "off"
    ws_off.mkdir()
    sim = _build(ws_off, "jncoff", exe, with_jnc=False)
    sim.write_simulation(silent=True)
    _run(sim)
    assert _njunctions(ws_off, "jncoff") == 0, (
        "JNC6 absent: no junction should be detected"
    )
    fja_off = _flowja(ws_off, "jncoff")
    assert np.any(np.abs(fja_off) > 1e-6), (
        "JNC6 absent: reach-reach connections should carry flow (unmasked), "
        f"got all-zero {fja_off}"
    )


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        test_jnc_toggle(Path(d))
        print("toggle test passed")
