"""
Generate AMG vs ILU(0) convergence figure for the SuppTechInfo AMG chapter.

Replicates the spirit of Mehl and Hill (2001) Problem 1 (Figure 1) using the
MODFLOW 6 API (libmf6) for direct per-iteration head access.

Problem setup
-------------
2D confined steady-state GWF model, 201x201 cells (N = 200 interior intervals),
uniform K = 1, constant-head boundaries of zero on all four sides.  The true
solution is h = 0 everywhere.  Two models are run per preconditioner (four
total), one for each of the two Fourier error modes:

  h_0(i, j) = sin(1 *pi*i/200) * sin(1 *pi*j/200)   (smooth,      k = 1)
  h_0(i, j) = sin(15*pi*i/200) * sin(15*pi*j/200)   (oscillatory, k = 15)

The combined (k = 1 + k = 15) amplitude is not run as a separate model.
Because the governing equation is linear and the two modes are orthogonal, the
combined error amplitude at each iteration is found from the two separate runs
by combining them as a root of the sum of squares:

  a_both(n) = sqrt(a_k1(n)^2 + a_k15(n)^2)

Each call to mf6.solve() performs one outer Picard iteration containing exactly
one inner CG step (INNER_MAXIMUM 1).  With RELAXATION_FACTOR 0 and
LINEAR_ACCELERATION CG the preconditioner applies an A-norm-optimal step
length, so the per-iteration amplitude reduction directly measures how
effectively each preconditioner attenuates the corresponding error mode.

Because the true solution is h = 0, the Fourier mode amplitude of the current
head equals the amplitude of the corresponding error component.

Figure panels (following Mehl and Hill (2001) Figure 1)
-------------------------------------------------------
  A  Initial head along the center row (row 101 of 201, index nrow//2=100)
     showing the smooth (k=1), oscillatory (k=15), and combined (k=1+k=15)
     components and the true solution h=0
  B  ILU(0) preconditioner: error amplitude vs. outer iteration for
     k=1, k=15, and the derived k=1+k=15 combined amplitude
  C  AMG preconditioner (CG, adaptive omega): same quantities as B

Output: doc/SuppTechInfo/Figures/amg-convergence.pdf

References
----------
Mehl, S.W. and Hill, M.C., 2001, MODFLOW-2000, The U.S. Geological Survey
  Modular Ground-Water Model -- Documentation of the Local Multigrid (LMG)
  Package: U.S. Geological Survey Open-File Report 01-177, 19 p.
"""

import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as fstyles
import matplotlib
import matplotlib.patheffects as pe
import numpy as np
from modflowapi import ModflowApi

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
repo_root = Path(__file__).resolve().parents[3]
libmf6 = repo_root / "bin" / "libmf6.dylib"
mf6_exe = repo_root / "bin" / "mf6"
figpth = Path(__file__).resolve().parent.parent / "Figures"
figpth.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Problem geometry
# ---------------------------------------------------------------------------
nlay, nrow, ncol = 1, 201, 201
delr = delc = 1.0
top = 1.0
botm = 0.0
K = 1.0

N1 = nrow - 1  # 200
N2 = ncol - 1  # 200

KX_SMOOTH, KY_SMOOTH = 1, 1  # smooth 2-D mode (k = 1 in both directions)
KX_OSCILL, KY_OSCILL = 15, 15  # oscillatory mode (k = 15, as in Mehl & Hill)

N_ITER = 200  # outer iterations to record per solver

# Build mode arrays
_i = np.arange(nrow)[:, np.newaxis]  # (201, 1)
_j = np.arange(ncol)[np.newaxis, :]  # (1, 201)
smooth_2d = np.sin(KX_SMOOTH * np.pi * _i / N1) * np.sin(KY_SMOOTH * np.pi * _j / N2)
oscill_2d = np.sin(KX_OSCILL * np.pi * _i / N1) * np.sin(KY_OSCILL * np.pi * _j / N2)
ic_2d = smooth_2d + oscill_2d  # (201, 201)

# Two initial conditions; combined amplitude found from these (root-sum-square)
smooth_ic = smooth_2d[np.newaxis, :, :]  # k=1 only   (1, 201, 201)
oscill_ic = oscill_2d[np.newaxis, :, :]  # k=15 only  (1, 201, 201)

# Colors consistent across all panels
C_SMOOTH = "C0"  # blue  — k = 1
C_OSCILL = "C3"  # red   — k = 15
C_BOTH = "C2"  # green — k = 1 + k = 15


# ---------------------------------------------------------------------------
# Mode-amplitude helper
# ---------------------------------------------------------------------------


def mode_amplitude(h_2d, kx, ky):
    """Compute the amplitude of Fourier mode (kx, ky) in h_2d."""
    phi = np.sin(kx * np.pi * _i / N1) * np.sin(ky * np.pi * _j / N2)
    return float(abs(np.sum(h_2d * phi)) / np.sum(phi**2))


# ---------------------------------------------------------------------------
# Model builder
# ---------------------------------------------------------------------------


def build_model(ws, name, use_amg, ic):
    """Write a 2-D steady-state GWF model with CHD = 0 on all four edges."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=str(mf6_exe),
        sim_ws=str(ws),
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=1.0e-30,  # never satisfy dvclose early
        outer_maximum=N_ITER + 100,  # headroom so MODFLOW does not abort
        inner_dvclose=1.0e-30,
        inner_maximum=1,  # exactly one CG step per outer iter
        rcloserecord="1.0e-30 strict",
        linear_acceleration="cg",
        preconditioner_levels=0 if not use_amg else 10,
        relaxation_factor=0.0,  # ILU(0) (not MILU) or AMG adaptive omega
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
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
    flopy.mf6.ModflowGwfic(gwf, strt=ic)
    flopy.mf6.ModflowGwfnpf(gwf, k=K, icelltype=0)

    # CHD = 0 on all four boundary edges
    chd_cells = []
    for j in range(ncol):
        chd_cells.append([(0, 0, j), 0.0])
        chd_cells.append([(0, nrow - 1, j), 0.0])
    for i in range(1, nrow - 1):
        chd_cells.append([(0, i, 0), 0.0])
        chd_cells.append([(0, i, ncol - 1), 0.0])
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_cells})

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST")],
    )
    sim.register_ims_package(ims, [gwf.name])
    sim.write_simulation()

    # Allow simulation to continue past a non-converged time step so
    # finalize_time_step() succeeds after our intentional early stop.
    nam_path = ws / "mfsim.nam"
    txt = nam_path.read_text()
    if "CONTINUE" not in txt.upper():
        txt = txt.replace("BEGIN Options", "BEGIN Options\n  CONTINUE")
        nam_path.write_text(txt)

    if use_amg:
        ims_path = ws / f"{name}.ims"
        txt = ims_path.read_text()
        txt = txt.replace("END LINEAR", "  PRECONDITIONER_TYPE  AMG\nEND LINEAR")
        ims_path.write_text(txt)


# ---------------------------------------------------------------------------
# API-driven solver loop
# ---------------------------------------------------------------------------


def run_and_collect(ws, name):
    """
    Run N_ITER outer iterations via the MODFLOW 6 API, returning a list of
    (nrow x ncol) head arrays: index 0 is the IC, index k is after step k.
    """
    mf6 = ModflowApi(str(libmf6), working_directory=str(ws))
    mf6.initialize()

    head_tag = mf6.get_var_address("X", name.upper())

    dt = mf6.get_time_step()
    mf6.prepare_time_step(dt)
    mf6.prepare_solve()

    h0 = mf6.get_value(head_tag).reshape(nrow, ncol).copy()
    history = [h0]

    for _ in range(N_ITER):
        has_converged = mf6.solve()
        h = mf6.get_value(head_tag).reshape(nrow, ncol).copy()
        history.append(h)
        if has_converged:
            break

    mf6.finalize_solve()
    mf6.finalize_time_step()
    mf6.finalize()

    return history


# ---------------------------------------------------------------------------
# Build, run, collect
# ---------------------------------------------------------------------------

CASES = [
    ("amg_k1", True, smooth_ic),
    ("amg_k15", True, oscill_ic),
    ("ilu_k1", False, smooth_ic),
    ("ilu_k15", False, oscill_ic),
]

print("Building and running models...")
histories = {}
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for name, use_amg, ic in CASES:
        ws = tmp / name
        ws.mkdir()
        build_model(ws, name, use_amg, ic)
        print(f"  Running {name} ({nrow}x{ncol}, up to {N_ITER} iters)...")
        histories[name] = run_and_collect(ws, name)


# Compute per-iteration mode amplitudes.
# The PDE is linear and the two modes are orthogonal, so the combined
# (k=1+k=15) amplitude can be found from the separate runs as a root of the
# sum of squares, without needing a third model per solver.
def get_amps(hist, kx, ky):
    return np.array([mode_amplitude(h, kx, ky) for h in hist])


amg_k1_amp = get_amps(histories["amg_k1"], KX_SMOOTH, KY_SMOOTH)
amg_k15_amp = get_amps(histories["amg_k15"], KX_OSCILL, KY_OSCILL)
amg_both_amp = np.sqrt(amg_k1_amp**2 + amg_k15_amp**2)

ilu_k1_amp = get_amps(histories["ilu_k1"], KX_SMOOTH, KY_SMOOTH)
ilu_k15_amp = get_amps(histories["ilu_k15"], KX_OSCILL, KY_OSCILL)
ilu_both_amp = np.sqrt(ilu_k1_amp**2 + ilu_k15_amp**2)

amg_iters = np.arange(len(histories["amg_k1"]))
ilu_iters = np.arange(len(histories["ilu_k1"]))
xmax = max(len(amg_iters), len(ilu_iters)) - 1

for tag, k1, k15, both in [
    ("AMG", amg_k1_amp, amg_k15_amp, amg_both_amp),
    ("ILU(0)", ilu_k1_amp, ilu_k15_amp, ilu_both_amp),
]:
    print(
        f"  {tag}: k=1 {k1[0]:.3f}->{k1[-1]:.2e}"
        f"  k=15 {k15[0]:.3f}->{k15[-1]:.2e}"
        f"  both {both[0]:.3f}->{both[-1]:.2e}"
    )

# ILU k=1 stagnation analysis
print("\nILU(0) k=1 amplitude at selected iterations:")
for it in [0, 5, 10, 25, 50, 100, 150, 200]:
    idx = min(it, len(ilu_k1_amp) - 1)
    print(f"  iter {it:3d}: {ilu_k1_amp[idx]:.4e}")

# ---------------------------------------------------------------------------
# Figure (following Mehl and Hill 2001, Figure 1 style)
# ---------------------------------------------------------------------------

width = 6.8  # inches, USGS two-column width
dpi = 300

with fstyles.USGSPlot():
    # A centred at top, B and C side-by-side on the bottom.
    # "." marks empty grid cells; equal-width columns keep B and C the same
    # size as A so the three panels share a consistent visual weight.
    fig, axes = plt.subplot_mosaic(
        [["A", "A", "A", "A"], ["B", "B", "C", "C"]],
        figsize=(width, width * 0.62),
        layout="constrained",
    )
    axes["C"].sharey(axes["B"])

    fstyles.heading(axes["A"], idx=0)
    fstyles.heading(axes["B"], idx=1)
    fstyles.heading(axes["C"], idx=2)

    # --- Panel A: initial head components along center row ---
    ax = axes["A"]
    center = nrow // 2
    x = np.arange(ncol)
    halo = [pe.withStroke(linewidth=3.5, foreground="white"), pe.Normal()]
    (h_both,) = ax.plot(
        x,
        ic_2d[center, :],
        color=C_BOTH,
        ls="-",
        lw=3.0,
        label=rf"$k=1 + k={KX_OSCILL}$ (smooth + oscillatory)",
    )
    (h_k15,) = ax.plot(
        x,
        oscill_2d[center, :],
        color=C_OSCILL,
        ls=":",
        lw=1.0,
        label=rf"$k={KX_OSCILL}$ (oscillatory)",
        path_effects=halo,
    )
    (h_k1,) = ax.plot(
        x,
        smooth_2d[center, :],
        color=C_SMOOTH,
        ls="-",
        lw=1.0,
        label=r"$k=1$ (smooth)",
        path_effects=halo,
    )
    h_true = ax.axhline(
        0,
        color="k",
        lw=0.8,
        ls="-",
        label=r"True solution ($h = 0$)",
        path_effects=halo,
    )
    ax.set_xlabel("Column")
    ax.set_ylabel("Head, in meters")
    ax.set_title("Initial head components")
    ax.set_xlim([0, ncol - 1])
    ax.set_ylim([-1.5, 3.2])
    ax.legend(
        handles=[h_true, h_k1, h_k15, h_both],
        fontsize=6.5,
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.99),
        framealpha=0.9,
    )

    # y-axis limits for B and C
    all_vals = np.concatenate(
        [
            amg_k1_amp,
            amg_k15_amp,
            amg_both_amp,
            ilu_k1_amp,
            ilu_k15_amp,
            ilu_both_amp,
        ]
    )
    all_vals = all_vals[all_vals > 0]
    ylo = 1e-10
    yhi = 2.0  # show full range from ~sqrt(2) initial combined amplitude

    # --- Panel B: ILU(0) error-component convergence ---
    ax = axes["B"]
    # Plot k=1+k=15 first (thick) so k=1 renders on top and is clearly visible
    ax.semilogy(ilu_iters, ilu_both_amp, color=C_BOTH, ls="-", lw=3.0)
    ax.semilogy(ilu_iters, ilu_k15_amp, color=C_OSCILL, ls="--", lw=1.2)
    ax.semilogy(ilu_iters, ilu_k1_amp, color=C_SMOOTH, ls="-", lw=1.2)
    ax.set_xlabel("Outer iteration")
    ax.set_ylabel("Error, in meters")
    ax.set_title("ILU(0)")
    ax.set_xlim([0, xmax])
    ax.set_ylim([ylo, yhi])

    # --- Panel C: AMG error-component convergence ---
    ax = axes["C"]
    ax.semilogy(amg_iters, amg_both_amp, color=C_BOTH, ls="-", lw=3.0)
    ax.semilogy(amg_iters, amg_k15_amp, color=C_OSCILL, ls="--", lw=1.2)
    ax.semilogy(amg_iters, amg_k1_amp, color=C_SMOOTH, ls="-", lw=1.2)
    ax.set_xlabel("Outer iteration")
    ax.set_title(r"AMG (CG, adaptive $\omega$)")
    ax.set_xlim([0, xmax])

    fig.align_ylabels([axes["A"], axes["B"]])

    outfile = figpth / "amg-convergence.pdf"
    fig.savefig(outfile, dpi=dpi)
    print(f"\nSaved {outfile}")
