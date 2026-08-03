"""
Tests the AMG preconditioner (PRECONDITIONER_TYPE AMG in the LINEAR block)
by running identical GWF problems with AMG and ILU0 and verifying that the
computed heads are equivalent.  Cases exercised:

  - CG with fixed omega (RELAXATION_FACTOR 1.0, undamped V-cycle)
  - BiCGSTAB with fixed omega
  - CG with adaptive omega (RELAXATION_FACTOR 0.0, adaptive optimal damping)
    [ims_amg_cg_adp]
  - CG with NUMBER_SMOOTHING_ITERATIONS 1  [ims_amg_nsmooth1]
  - CG with NUMBER_SMOOTHING_ITERATIONS 4  [ims_amg_nsmooth4]
  - Newton-Raphson GWF + BiCGSTAB + AMG fixed omega + nsmooth=1  [ims_amg_nwt_ns1]
  - Newton-Raphson GWF + BiCGSTAB + AMG fixed omega  [ims_amg_nwt_fix]
  - CG + AMG + adaptive omega + PRECONDITIONER_STRENGTH_THRESHOLD 0.25  [ims_amg_s25]
  - CG + AMG + adaptive omega + PRECONDITIONER_STRENGTH_THRESHOLD 0.50  [ims_amg_s50]
  - Newton-Raphson GWF + BiCGSTAB + AMG fixed omega
                       + PRECONDITIONER_STRENGTH_THRESHOLD 0.25  [ims_amg_nwt_s25]
  - Newton-Raphson GWF + BiCGSTAB + AMG fixed omega
                       + PRECONDITIONER_STRENGTH_THRESHOLD 0.50  [ims_amg_nwt_s50]
  - CG + SMOOTHER_TYPE ILU0_ALL  [ims_amg_ilu0all]

The confined GWF model uses a heterogeneous K field: layer base K values span
three orders of magnitude and a checkerboard pattern within each layer places
small-K cells directly next to large-K cells (100:1 contrast).

The Newton-Raphson cases use a single unconfined layer (icelltype=1).  Their
matrix is nonsymmetric, so BiCGSTAB is required.  Adaptive omega
(RELAXATION_FACTOR 0) is not used here: the adaptive value changes at every
inner iteration, which makes the preconditioner change too.  BiCGSTAB assumes
a fixed preconditioner, so a changing one prevents convergence.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from amg_test_utils import base_ims_kwargs, build_checker_gwf
from framework import TestFramework

paktest = "ims"
cmp_prefix = "mf6"  # subdirectory used for the ILU0 reference run

cases = [
    "ims_amg_cg",
    "ims_amg_bicgstab",
    "ims_amg_cg_adp",  # adaptive omega (RELAXATION_FACTOR 0); max 16 chars
    "ims_amg_nsmooth1",  # NUMBER_SMOOTHING_ITERATIONS 1
    "ims_amg_nsmooth4",  # NUMBER_SMOOTHING_ITERATIONS 4, CG
    "ims_amg_nwt_ns1",  # Newton-Raphson + BiCGSTAB + AMG fixed omega + nsmooth=1
    "ims_amg_nwt_fix",  # Newton-Raphson + BiCGSTAB + AMG fixed omega
    "ims_amg_s25",  # CG + AMG + adaptive omega + sthresh=0.25
    "ims_amg_s50",  # CG + AMG + adaptive omega + sthresh=0.50
    "ims_amg_nwt_s25",  # Newton-Raphson + BiCGSTAB + AMG fixed omega + sthresh=0.25
    "ims_amg_nwt_s50",  # Newton-Raphson + BiCGSTAB + AMG fixed omega + sthresh=0.50
    "ims_amg_ilu0all",  # CG + SMOOTHER_TYPE ILU0_ALL (ILU0 smoother at all AMG levels)
]
linear_accel = [
    "cg",
    "bicgstab",
    "cg",
    "cg",
    "cg",
    "bicgstab",
    "bicgstab",
    "cg",
    "cg",
    "bicgstab",
    "bicgstab",
    "cg",
]
relax_factors = [
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    1.0,
    1.0,
    0.0,
    0.0,
    1.0,
    1.0,
    0.0,
]
nsmooth_values = [
    None,
    None,
    None,
    1,
    4,
    1,
    None,
    None,
    None,
    None,
    None,
    None,
]
# nsmooth=4 with CG may need more inner iterations on the
# heterogeneous K field; all other cases converge well within 500
inner_max_values = [
    500,
    500,
    500,
    500,
    1000,
    500,
    500,
    500,
    500,
    500,
    500,
    500,
]
# Newton-Raphson cases produce a nonsymmetric Jacobian (requires BiCGSTAB)
is_newton = [
    False,
    False,
    False,
    False,
    False,
    True,
    True,
    False,
    False,
    True,
    True,
    False,
]
# PRECONDITIONER_STRENGTH_THRESHOLD for the sthresh cases; None = default (0)
sthresh_values = [
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    0.25,
    0.50,
    0.25,
    0.50,
    None,
]
# AMG smoother type: None = default (ILU0), "ILU0_ALL" = ILU0 at all levels
smoother_values = [
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    "ILU0_ALL",
]

# spatial discretization ── 3-layer confined model, anisotropic K
nrow, ncol = 10, 20
delr, delc = 100.0, 100.0

# spatial discretization ── single unconfined layer for Newton-Raphson cases
nlay_nwt = 1
top_nwt = 10.0
botm_nwt = [0.0]
strt_nwt = 7.0
chd_left_nwt = 9.0
chd_right_nwt = 3.0

# tolerance: AMG and ILU0 must produce identical discrete solutions
htol = 1.0e-6


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _build_amg_sim(
    name,
    ws,
    accel,
    relax,
    nsmooth=None,
    inner_maximum=500,
    sthresh=None,
    smoother_type=None,
):
    """Build and write a GWF simulation that uses the AMG preconditioner.

    The AMG keywords are set directly through the flopy ModflowIms interface.
    nsmooth, if given, sets NUMBER_SMOOTHING_ITERATIONS.
    sthresh, if given, sets PRECONDITIONER_STRENGTH_THRESHOLD.
    """
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim,
        **base_ims_kwargs(accel, relax, inner_maximum),
        preconditioner_type="AMG",
        number_smoothing_iterations=nsmooth,
        preconditioner_strength_threshold=sthresh,
        smoother_type=smoother_type,
    )
    gwf = build_checker_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])

    sim.write_simulation()
    return sim


def _build_ilu0_sim(name, ws, accel, relax):
    """Build and write a GWF simulation that uses the ILU preconditioner
    for comparison against the AMG run."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(sim, **base_ims_kwargs(accel, relax))
    gwf = build_checker_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])

    sim.write_simulation()
    return sim


def _build_newton_gwf(sim, name):
    """Add an unconfined (Newton-Raphson) GWF model to *sim*.

    icelltype=1 makes cells convertible (unconfined).  Newton produces a
    nonsymmetric matrix that requires BiCGSTAB.
    """
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay_nwt,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top_nwt,
        botm=botm_nwt,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt_nwt)

    rows = np.arange(nrow)
    cols = np.arange(ncol)
    checker = np.where((rows[:, None] + cols[None, :]) % 2 == 0, 10.0, 0.1)
    k = checker[None, :, :]  # (1, nrow, ncol)

    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k)

    spd = [[(0, irow, 0), chd_left_nwt] for irow in range(nrow)]
    spd += [[(0, irow, ncol - 1), chd_right_nwt] for irow in range(nrow)]
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return gwf


def _build_newton_amg_sim(
    name, ws, relax, nsmooth=None, sthresh=None, smoother_type=None
):
    """Build and write a Newton GWF simulation that uses the AMG preconditioner."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim,
        **base_ims_kwargs("bicgstab", relax),
        preconditioner_type="AMG",
        number_smoothing_iterations=nsmooth,
        preconditioner_strength_threshold=sthresh,
        smoother_type=smoother_type,
    )
    gwf = _build_newton_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])

    sim.write_simulation()
    return sim


def _build_newton_ilu_sim(name, ws, relax):
    """Build and write a Newton GWF simulation that uses the ILU preconditioner
    for comparison against the AMG run."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(sim, **base_ims_kwargs("bicgstab", relax))
    gwf = _build_newton_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])

    sim.write_simulation()
    return sim


# ---------------------------------------------------------------------------
# TestFramework interface
# ---------------------------------------------------------------------------


def build_models(idx, test):
    name = cases[idx]
    accel = linear_accel[idx]
    relax = relax_factors[idx]
    nsmooth = nsmooth_values[idx]
    inner_maximum = inner_max_values[idx]
    sthresh = sthresh_values[idx]
    smoother = smoother_values[idx]
    cmp_ws = Path(test.workspace) / cmp_prefix

    if is_newton[idx]:
        amg_sim = _build_newton_amg_sim(
            name,
            test.workspace,
            relax,
            nsmooth=nsmooth,
            sthresh=sthresh,
            smoother_type=smoother,
        )
        ilu_sim = _build_newton_ilu_sim(name, cmp_ws, relax)
    else:
        amg_sim = _build_amg_sim(
            name,
            test.workspace,
            accel,
            relax,
            nsmooth=nsmooth,
            inner_maximum=inner_maximum,
            sthresh=sthresh,
            smoother_type=smoother,
        )
        ilu_sim = _build_ilu0_sim(name, cmp_ws, accel, relax)
    return amg_sim, ilu_sim


def check_output(idx, test):
    name = cases[idx]

    # AMG heads
    hds_amg = flopy.utils.HeadFile(Path(test.workspace) / f"{name}.hds").get_alldata()

    # ILU0 heads (reference)
    hds_ilu = flopy.utils.HeadFile(
        Path(test.workspace) / cmp_prefix / f"{name}.hds"
    ).get_alldata()

    assert hds_amg.shape == hds_ilu.shape, (
        f"Head array shapes differ: AMG {hds_amg.shape} vs ILU {hds_ilu.shape}"
    )

    max_diff = np.max(np.abs(hds_amg - hds_ilu))
    assert max_diff < htol, (
        f"AMG and ILU heads differ by {max_diff:.3e} m (tolerance {htol:.3e} m)"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        # build_models already wrote the files; overwrite=False stops the
        # framework from re-writing and clobbering the AMG IMS file
        overwrite=False,
        # no legacy-solver comparison here; check_output compares the AMG
        # and ILU heads directly
        compare=None,
    )
    test.run()


def test_amg_bicgstab_relax0_error(function_tmpdir, targets):
    """AMG + BICGSTAB + RELAXATION_FACTOR 0 must stop with an error.

    Adaptive damping changes the preconditioner every inner iteration, but
    BiCGSTAB needs a fixed one.  The IMS solver detects this combination and
    MODFLOW 6 stops before running the time loop.
    """
    name = "ims_amg_bgs_r0_err"

    def _build(test):
        sim = flopy.mf6.MFSimulation(
            sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
        )
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
        ims = flopy.mf6.ModflowIms(
            sim, **base_ims_kwargs("bicgstab", relax=0.0), preconditioner_type="AMG"
        )
        gwf = build_checker_gwf(sim, name)
        sim.register_ims_package(ims, [gwf.name])
        sim.write_simulation()
        return sim

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=_build,
        check=None,
        compare=None,
        overwrite=False,
        xfail=True,
    )
    test.run()
