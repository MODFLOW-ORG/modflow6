"""
Integration test for the AMG preconditioner on a large, two-dimensional
unconfined groundwater flow model that uses the Newton-Raphson formulation
and a complex, heterogeneous hydraulic-conductivity (K) field.

The setup is inspired by Problem 2 from the Mehl and Hill (2001) MODFLOW-2000
LMG solver documentation, a large 2D model with a very complicated K field
(Konikow, L., and Hornberger, G., U.S. Geological Survey, written
communication, 2000).

We use an unconfined (Newton-Raphson) formulation instead of confined because:

  1. The Newton Jacobian is slightly nonsymmetric, so it needs BiCGSTAB, and
     that makes AMG more competitive against ILU.
  2. Transmissivity T = K * b(h) changes with the water-table height b(h),
     giving variation across many scales that AMG handles well.
  3. An unconfined sand-and-gravel aquifer matches the reference site.

The test builds a single-layer (2D areal) unconfined GWF model with:

  - A spatially correlated, log-normal K field built with an FFT-based
    random-field generator (see _build_k_field).
  - Two correlation scales: a regional one (range 20 cells / 500 m, 70 % of
    the variance) and a local one (range 6 cells / 150 m, 30 %).
  - Mean log10(K) = 1.0 (K_mean = 10 m/d), std = 0.65 in log10 space, so K
    spans about 4 orders of magnitude.
  - Horizontal anisotropy (K22 = 0.3 * K11) for extra K contrast between the
    x and y directions.
  - Grid: nlay=1, nrow=80, ncol=100 (8,000 active cells; dx=dy=25 m;
    domain 2,000 m x 2,500 m).
  - Unconfined flow; west CHD = 9 m, east CHD = 2 m (aquifer top = 10 m,
    bottom = 0 m); saturated thickness ranges from 2 to 9 m, so transmissivity
    spans more than 5 orders of magnitude.
  - Initial heads interpolated linearly between the two CHD values. Starting
    far from the solution makes the first-step residual huge, so AMG cannot
    give a useful correction on the first Newton iteration.

Both AMG cases use BiCGSTAB + RELAXATION_FACTOR = 1 (adaptive omega does not
work with BiCGSTAB; see ImsLinear.f90).
AMG hierarchy depth: ceil(log2(8000/30)) ~ 9; PRECONDITIONER_LEVELS = 9.

Two Newton-Raphson AMG configurations are compared against a MILUT reference
(PRECONDITIONER_LEVELS = 5):
  - NEWTON + BiCGSTAB + AMG                              [lmg2_amg_nwt]
  - NEWTON + BiCGSTAB + AMG + PRECONDITIONER_STRENGTH_THRESHOLD 0.10
                                                         [lmg2_amg_nwt_s]

The second case checks whether limiting AMG coarsening to strongly connected
cell pairs improves the coarse grids for this heterogeneous K field, where
weak connections between low-K and high-K cells can hurt accuracy.

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
from amg_test_utils import count_inner_iters
from framework import TestFramework

paktest = "ims"
cmp_prefix = "mf6"  # subdirectory used for the ILU reference run

cases = [
    "lmg2_amg_nwt",  # NEWTON + BiCGSTAB + AMG (default sthresh)
    "lmg2_amg_nwt_s",  # NEWTON + BiCGSTAB + AMG + STRENGTH_THRESHOLD 0.10
]
# PRECONDITIONER_STRENGTH_THRESHOLD for each case (None = default 0.25)
sthresh_values = [None, 0.10]

# ---------------------------------------------------------------------------
# Grid and aquifer geometry
# ---------------------------------------------------------------------------

# Single-layer (2D areal) unconfined model — domain 2,000 m x 2,500 m
nlay, nrow, ncol = 1, 80, 100
delr, delc = 25.0, 25.0  # m
top = 10.0  # aquifer top elevation (m)
botm = [0.0]  # aquifer bottom (m); 10 m thick

# Constant-head boundaries — water table ranges from 2 to 9 m, so saturated
# thickness b(h) ranges from 2 to 9 m and T = K*b varies by > 5 orders of
# magnitude across the domain when combined with the K heterogeneity.
chd_west = 9.0  # m (near-top; mountain-front recharge)
chd_east = 2.0  # m (near-bottom; basin discharge)

# Initial heads: linear gradient from west CHD to east CHD.  A constant
# initial head far from the solution produces a huge initial residual that
# makes the AMG preconditioner ineffective on the first Newton step.
strt = np.linspace(chd_west, chd_east, ncol)  # shape (ncol,); broadcast in _build_gwf

# K22 = k22_ratio * K11: 3.3:1 horizontal anisotropy
k22_ratio = 0.30

# AMG hierarchy depth: ceil(log2(8000/30)) ≈ 9
amg_levels = 9
# MILUT fill level for the ILU reference run
ilu_levels = 5

# Head tolerance for AMG-vs-ILU comparison.  Both solvers converge to
# outer_dvclose = 1e-6 m, so the solutions should agree to within ~1e-4 m.
htol = 1.0e-4


# ---------------------------------------------------------------------------
# K-field generation — FFT-based random field
# ---------------------------------------------------------------------------


def _spectral_unit_field(nrow, ncol, range_cells, seed):
    """Generate a smooth, random 2D field (zero mean, unit variance).

    Uses an FFT-based method in which nearby cells are correlated over a
    distance of about range_cells, following an exponential covariance
    C(h) = exp(-h/a).  The grid is padded by one correlation length on each
    side to avoid wrap-around edge effects, then cropped back to
    (nrow, ncol).
    """
    rng = np.random.default_rng(seed)

    # Pad by one correlation length to reduce periodic-boundary artifacts
    pad = max(int(2 * range_cells), 20)
    Nr, Nc = nrow + 2 * pad, ncol + 2 * pad

    # Radial frequency magnitude at each grid point (cycles per cell)
    fr = np.fft.fftfreq(Nr)
    fc = np.fft.fftfreq(Nc)
    Fr, Fc = np.meshgrid(fr, fc, indexing="ij")
    f_mag = np.sqrt(Fr**2 + Fc**2)

    # 2D power spectral density for the isotropic exponential covariance
    # C(h) = exp(-h/a)  <==>  P(f) = 2*pi*a^2 / (1 + (2*pi*a*f)^2)^(3/2)
    a = float(range_cells)
    psd = 2.0 * np.pi * a**2 / (1.0 + (2.0 * np.pi * a * f_mag) ** 2) ** 1.5
    psd[0, 0] = 0.0  # zero mean

    # Scale complex white noise by sqrt(PSD) and transform to spatial domain
    noise = rng.standard_normal((Nr, Nc)) + 1j * rng.standard_normal((Nr, Nc))
    field = np.fft.ifft2(noise * np.sqrt(psd)).real

    # Crop padding and standardize to zero mean, unit variance
    field = field[pad : pad + nrow, pad : pad + ncol]
    return (field - field.mean()) / field.std()


def _build_k_field(nrow, ncol, seed=42):
    """Build a heterogeneous K11 field by combining two correlation scales.

      - Regional (range = 20 cells / 500 m): 70 % of the log10(K) variance
      - Local    (range =  6 cells / 150 m): 30 % of the log10(K) variance

    Mean log10(K) = 1.0 (K_mean = 10 m/d), std = 0.65 in log10 space, so K
    spans about 4 orders of magnitude.  Combined with the 4.5:1 range in
    saturated thickness, transmissivity varies by more than 5 orders of
    magnitude.
    """
    z_regional = _spectral_unit_field(nrow, ncol, range_cells=20, seed=seed)
    z_local = _spectral_unit_field(nrow, ncol, range_cells=6, seed=seed + 1)

    # Weights for the two scales: w1² + w2² = 1
    w1, w2 = np.sqrt(0.70), np.sqrt(0.30)
    z = w1 * z_regional + w2 * z_local
    z = (z - z.mean()) / z.std()  # re-standardize after combination

    log10k_mean, log10k_std = 1.0, 0.65
    return 10.0 ** (log10k_mean + log10k_std * z)


# ---------------------------------------------------------------------------
# Model-building helpers
# ---------------------------------------------------------------------------


def _build_gwf(sim, name):
    """Add a 2D unconfined Newton-Raphson GWF model with a complex K field."""
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )

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
    # strt is shape (ncol,); broadcast to (nlay, nrow, ncol)
    strt_arr = np.broadcast_to(strt, (nlay, nrow, ncol)).copy()
    flopy.mf6.ModflowGwfic(gwf, strt=strt_arr)

    k11 = _build_k_field(nrow, ncol)
    k22 = k11 * k22_ratio  # K in y-direction: 3.3:1 horizontal anisotropy

    # icelltype=1: convertible (unconfined); Newton handles the nonlinearity
    # in effective conductance as the water table rises and falls.
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k11, k22=k22)

    # Constant-head boundaries on the west and east faces
    spd = [[(0, irow, 0), chd_west] for irow in range(nrow)]
    spd += [[(0, irow, ncol - 1), chd_east] for irow in range(nrow)]
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return gwf


def _ims_kwargs(levels, name, inner_maximum=500):
    # All Newton cases use BiCGSTAB (nonsymmetric Jacobian) with fixed omega
    # (RELAXATION_FACTOR = 1); adaptive omega is incompatible with BiCGSTAB.
    return dict(
        print_option="ALL",
        outer_dvclose=1.0e-6,
        outer_maximum=100,
        inner_dvclose=1.0e-12,
        inner_maximum=inner_maximum,
        rcloserecord="1.0e-6 strict",
        linear_acceleration="bicgstab",
        preconditioner_levels=levels,
        relaxation_factor=1.0,
        csv_outer_output_filerecord=f"{name}_outer.csv",
    )


def _build_amg_sim(name, ws, sthresh=None):
    """Build and write an unconfined Newton GWF simulation using AMG."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(
        sim,
        **_ims_kwargs(amg_levels, name),
        preconditioner_type="AMG",
        preconditioner_strength_threshold=sthresh,
    )
    gwf = _build_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()
    return sim


def _build_ilu_sim(name, ws):
    """Build and write an unconfined Newton GWF simulation using ILU (reference)."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims = flopy.mf6.ModflowIms(sim, **_ims_kwargs(ilu_levels, name))
    gwf = _build_gwf(sim, name)
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()
    return sim


# ---------------------------------------------------------------------------
# TestFramework interface
# ---------------------------------------------------------------------------


def build_models(idx, test):
    name = cases[idx]
    sthresh = sthresh_values[idx]
    cmp_ws = Path(test.workspace) / cmp_prefix

    amg_sim = _build_amg_sim(name, test.workspace, sthresh=sthresh)
    ilu_sim = _build_ilu_sim(name, cmp_ws)
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

    amg_iters = count_inner_iters(Path(test.workspace) / f"{name}_outer.csv")
    ilu_iters = count_inner_iters(
        Path(test.workspace) / cmp_prefix / f"{name}_outer.csv"
    )
    assert amg_iters < ilu_iters, (
        f"AMG total inner iterations ({amg_iters}) >= ILU ({ilu_iters})"
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
