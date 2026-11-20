"""
Iteration-count performance tests for the IMS AMG preconditioner.

This test verifies that the AMG preconditioner actually reduces linear solver
work relative to a baseline, not merely that it produces correct heads.  The
comparison is made on a 3-layer confined GWF model with a 100:1 heterogeneous
K field:

  AMG vs ILU(0):
     AMG (5 levels, adaptive damping via RELAXATION_FACTOR 0) should need
     fewer CG inner iterations than zero-fill ILU(0) (RELAXATION_FACTOR 0,
     PRECONDITIONER_LEVELS 0) on this heterogeneous system.  RELAXATION_FACTOR
     0 is used on both sides: for AMG it turns on adaptive damping (recommended
     for CG/SPD); for ILU it selects plain ILU(0) instead of the stronger
     MILU(0) (RELAXATION_FACTOR 1.0), which competes too well with AMG and can
     hide the benefit on small problems.

The test writes the IMS CSV inner output (csv_inner_output) and reads the
``total_inner_iterations`` column to count solver work.
"""

from pathlib import Path

import flopy
from amg_test_utils import base_ims_kwargs, build_checker_gwf, count_inner_iters
from framework import TestFramework


def _build_amg_sim(name, ws, accel, relax, inner_maximum=500, preconditioner_levels=5):
    """Build a GWF simulation with the AMG preconditioner and CSV inner output."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim,
        csv_inner_output_filerecord=[[f"{name}.inner.csv"]],
        **base_ims_kwargs(accel, relax, inner_maximum, preconditioner_levels),
        preconditioner_type="AMG",
    )
    gwf = build_checker_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()
    return sim


def _build_ilu0_sim(name, ws, accel, relax):
    """Build a GWF simulation with plain zero-fill ILU(0) (preconditioner_levels=0,
    relaxation_factor=0) and CSV inner output, as a baseline to compare with AMG.

    RELAXATION_FACTOR 0 selects ILU(0) rather than MILU(0); the modified
    variant (RELAXATION_FACTOR 1) is much stronger and can hide AMG's benefit
    on small problems.
    """
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim,
        csv_inner_output_filerecord=[[f"{name}.inner.csv"]],
        **base_ims_kwargs(accel, relax, inner_maximum=500, preconditioner_levels=0),
    )
    gwf = build_checker_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()
    return sim


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_amg_fewer_iters_than_ilu0(function_tmpdir, targets):
    """AMG requires fewer CG inner iterations than ILU(0) on the heterogeneous K model.

    The confined 3-layer model has a 100:1 K contrast (checkerboard within
    each layer, plus 100x contrast between layer 1 and the confining unit),
    which makes the system hard to solve.  ILU(0) (zero-fill, RELAXATION_FACTOR
    0, PRECONDITIONER_LEVELS 0) cannot remove the smooth error between layers;
    AMG removes it through its coarse levels.

    AMG uses RELAXATION_FACTOR 0 (adaptive damping), the recommended setting
    for CG on SPD systems.  ILU also uses RELAXATION_FACTOR 0, giving plain
    ILU(0) rather than the stronger MILU(0) that RELAXATION_FACTOR 1 would
    produce.

    The check is strict: AMG must take fewer inner iterations, not just equal.
    """
    name_amg = "ims_amg_cg_it"  # 14 chars
    name_ilu = "ims_ilu0_cg_it"  # 14 chars

    ws_amg = Path(function_tmpdir) / "amg"
    ws_ilu = Path(function_tmpdir) / "ilu0"
    ws_amg.mkdir(exist_ok=True)
    ws_ilu.mkdir(exist_ok=True)

    def build_amg(test):
        return _build_amg_sim(
            name_amg, test.workspace, "cg", relax=0.0, preconditioner_levels=5
        )

    def build_ilu(test):
        return _build_ilu0_sim(name_ilu, test.workspace, "cg", relax=0.0)

    TestFramework(
        name=name_amg,
        workspace=ws_amg,
        targets=targets,
        build=build_amg,
        check=None,
        compare=None,
        overwrite=False,
    ).run()

    TestFramework(
        name=name_ilu,
        workspace=ws_ilu,
        targets=targets,
        build=build_ilu,
        check=None,
        compare=None,
        overwrite=False,
    ).run()

    iters_amg = count_inner_iters(ws_amg / f"{name_amg}.inner.csv")
    iters_ilu = count_inner_iters(ws_ilu / f"{name_ilu}.inner.csv")

    assert iters_amg > 0, "Could not read AMG inner iteration count from CSV"
    assert iters_ilu > 0, "Could not read ILU(0) inner iteration count from CSV"
    assert iters_amg < iters_ilu, (
        f"AMG ({iters_amg} inner iterations) should require fewer inner "
        f"iterations than ILU(0) ({iters_ilu} inner iterations) on the "
        f"heterogeneous K model"
    )
