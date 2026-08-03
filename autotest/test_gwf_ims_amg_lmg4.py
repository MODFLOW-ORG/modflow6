"""
Integration test for the AMG preconditioner on a large three-dimensional
confined groundwater flow model with a complex heterogeneous
hydraulic-conductivity field.

This test is inspired by Problem 4 from the Mehl and Hill (2001) MODFLOW-2000
LMG solver documentation.  Problem 4 is described as a large, three-dimensional
model with a complex heterogeneous hydraulic-conductivity field (O'Brien, G.M.
and D'Agnese, F.A., U.S. Geological Survey, written communication, 2000).  The
original model represents a regional basin-and-range groundwater system.
This test builds a synthetic stand-in that keeps the features of that problem
that matter to the solver:

  - 10 geologic layers alternating between high-K aquifer units and low-K
    confining units, with hydraulic conductivities spanning six orders of
    magnitude (0.0001 to 10 m/d)
  - Lateral K variability within each layer (alluvial-fan, basin-centre, and
    distal-fine facies zones) controlled by column position
  - Strong vertical anisotropy in confining units (K33/K11 as low as 0.001)
  - Steady-state confined flow with constant-head boundaries on the west and
    east faces representing mountain-front recharge and basin discharge
  - Grid: nlay=10, nrow=30, ncol=40 (12,000 active cells)

The AMG hierarchy depth for N=12,000 is ceil(log2(12000/30)) ≈ 9; the test
uses PRECONDITIONER_LEVELS = 10 to give the hierarchy one level of headroom.

Two AMG configurations are tested against a MILUT reference (PRECONDITIONER_LEVELS=5):

  - CG  + adaptive damping (RELAXATION_FACTOR 0, SPD)  [lmg4_amg_cg]
  - BiCGSTAB + fixed damping (RELAXATION_FACTOR 1)      [lmg4_amg_bgs]

References
----------
Mehl, S.W. and Hill, M.C., 2001, MODFLOW-2000, The U.S. Geological Survey
  Modular Ground-Water Model--Documentation of the Local Multigrid (LMG)
  Package: U.S. Geological Survey Open-File Report 01-177, 19 p.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "ims"
cmp_prefix = "mf6"  # subdirectory used for the ILU reference run

cases = [
    "lmg4_amg_cg",  # CG  + AMG + adaptive omega (RELAXATION_FACTOR 0)
    "lmg4_amg_bgs",  # BiCGSTAB + AMG + fixed omega (RELAXATION_FACTOR 1)
]
linear_accel = ["cg", "bicgstab"]
relax_factors = [0.0, 1.0]

# ---------------------------------------------------------------------------
# Grid and geology
# ---------------------------------------------------------------------------

# 10-layer confined model — nrow x ncol = 30 x 40, cell size 1 km
nlay, nrow, ncol = 10, 30, 40
delr, delc = 1000.0, 1000.0  # m

# Layer bottom elevations (m) — 10 geologic units
# alternating aquifer (A) and confining (C) units
top = 0.0
botm = [
    -50.0,  # layer  1 base: valley fill / alluvium          (A,  50 m)
    -100.0,  # layer  2 base: lacustrine clay                 (C,  50 m)
    -250.0,  # layer  3 base: upper carbonate aquifer         (A, 150 m)
    -300.0,  # layer  4 base: volcanic tuff                   (C,  50 m)
    -550.0,  # layer  5 base: main carbonate aquifer          (A, 250 m)
    -650.0,  # layer  6 base: crystalline semi-confining unit (C, 100 m)
    -800.0,  # layer  7 base: volcanic aquifer                (A, 150 m)
    -870.0,  # layer  8 base: mudstone                        (C,  70 m)
    -1050.0,  # layer  9 base: deep carbonate aquifer          (A, 180 m)
    -1200.0,  # layer 10 base: deep basement / low-K           (C, 150 m)
]

# Baseline horizontal K per layer (m/d)
# Aquifer and confining units span ~6 orders of magnitude, simulating the
# large property contrasts in Basin-and-Range regional flow systems.
k_baseline = np.array(
    [
        1.0,  # layer  1: valley fill (A)
        0.001,  # layer  2: lacustrine clay (C)
        5.0,  # layer  3: upper carbonate (A)
        0.0001,  # layer  4: volcanic tuff (C)
        10.0,  # layer  5: main carbonate (A)
        0.01,  # layer  6: crystalline semi-confining (C)
        2.0,  # layer  7: volcanic aquifer (A)
        0.0001,  # layer  8: mudstone (C)
        5.0,  # layer  9: deep carbonate (A)
        0.001,  # layer 10: deep basement (low-K)
    ]
)

# Vertical anisotropy ratio K33/K11 per layer
k33_ratio = np.array(
    [
        0.10,  # layer  1: valley fill
        0.01,  # layer  2: lacustrine clay
        0.10,  # layer  3: upper carbonate
        0.001,  # layer  4: volcanic tuff (strongly anisotropic)
        0.10,  # layer  5: main carbonate
        0.05,  # layer  6: crystalline semi-confining
        0.10,  # layer  7: volcanic aquifer
        0.001,  # layer  8: mudstone (strongly anisotropic)
        0.10,  # layer  9: deep carbonate
        0.05,  # layer 10: deep basement
    ]
)

# Regional head boundary conditions (m)
chd_west = 200.0  # mountain front / recharge area
chd_east = 50.0  # basin discharge (playas / springs)
strt = 125.0  # uniform initial head

# AMG hierarchy depth: ceil(log2(12000/30)) ≈ 9 → use 10 for headroom
amg_levels = 10
# MILUT/ILUT fill level for the ILU reference run
ilu_levels = 5

# Tolerance for AMG-vs-ILU head comparison (m)
htol = 1.0e-6


# ---------------------------------------------------------------------------
# K-field construction
# ---------------------------------------------------------------------------


def _build_k_field():
    """Build a heterogeneous (nlay, nrow, ncol) K array.

    K varies across columns to mimic a west-to-east change in sediment type:
    coarse, high-K deposits near the mountain front (cols 0 to ncol/5-1,
    K multiplier 3), average deposits in the basin centre (cols ncol/5 to
    4*ncol/5-1, K multiplier 1), and fine, low-K deposits at the far edge
    (cols 4*ncol/5 to ncol-1, K multiplier 0.2).
    """
    cols = np.arange(ncol)
    facies = np.where(
        cols < ncol // 5,
        3.0,
        np.where(cols >= 4 * ncol // 5, 0.2, 1.0),
    )

    k = np.empty((nlay, nrow, ncol))
    k33 = np.empty((nlay, nrow, ncol))
    for ilay in range(nlay):
        k[ilay, :, :] = k_baseline[ilay] * facies[np.newaxis, :]
        k33[ilay, :, :] = k[ilay, :, :] * k33_ratio[ilay]
    return k, k33


# ---------------------------------------------------------------------------
# Model-building helpers
# ---------------------------------------------------------------------------


def _build_gwf(sim, name):
    """Add a steady-state confined GWF model with west/east CHD boundaries."""
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    k, k33 = _build_k_field()
    flopy.mf6.ModflowGwfnpf(gwf, k=k, k33=k33)

    # Constant-head boundaries on the west face (recharge) and east face
    # (discharge) of all ten layers
    spd = [[(ilay, irow, 0), chd_west] for ilay in range(nlay) for irow in range(nrow)]
    spd += [
        [(ilay, irow, ncol - 1), chd_east]
        for ilay in range(nlay)
        for irow in range(nrow)
    ]
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return gwf


def _ims_kwargs(accel, relax, levels, inner_maximum=500):
    return dict(
        print_option="ALL",
        outer_dvclose=1.0e-8,
        outer_maximum=500,
        inner_dvclose=1.0e-12,
        inner_maximum=inner_maximum,
        rcloserecord="1.0 strict",
        linear_acceleration=accel,
        preconditioner_levels=levels,
        relaxation_factor=relax,
    )


def _build_amg_sim(name, ws, accel, relax):
    """Build a GWF simulation using the AMG preconditioner."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim, **_ims_kwargs(accel, relax, amg_levels), preconditioner_type="AMG"
    )
    gwf = _build_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()
    return sim


def _build_ilu_sim(name, ws, accel, relax):
    """Build a GWF simulation using the ILU preconditioner (reference)."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(sim, **_ims_kwargs(accel, relax, ilu_levels))
    gwf = _build_gwf(sim, name)
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
    cmp_ws = Path(test.workspace) / cmp_prefix

    amg_sim = _build_amg_sim(name, test.workspace, accel, relax)
    ilu_sim = _build_ilu_sim(name, cmp_ws, accel, relax)
    return amg_sim, ilu_sim


def check_output(idx, test):
    name = cases[idx]

    hds_amg = flopy.utils.HeadFile(Path(test.workspace) / f"{name}.hds").get_alldata()

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
        overwrite=False,
        compare=None,
    )
    test.run()
