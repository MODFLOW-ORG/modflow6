"""
Generate AMG vs ILU(0) convergence figure for a heterogeneous K field.

Companion to amg-convergence.py.  The problem setup is identical except that
the hydraulic conductivity field is a spatially uncorrelated lognormal random
field (sigma_log = 1.5), which raises the condition number of the system
matrix and amplifies the contrast between the ILU(0) and AMG preconditioners.

Problem setup
-------------
2D confined steady-state GWF model, 201x201 cells (N = 200 interior intervals),
lognormal K field (mean_log = 0, sigma_log = 1.5, seed = 42), constant-head
boundaries of zero on all four sides.  Because there is no flow gradient and
all boundary heads are zero, the true solution is h = 0 everywhere regardless
of the K field.  Two models are run per preconditioner (four total), one for
each of the two Fourier error modes:

  h_0(i, j) = sin(1 *pi*i/200) * sin(1 *pi*j/200)   (smooth,      k = 1)
  h_0(i, j) = sin(15*pi*i/200) * sin(15*pi*j/200)   (oscillatory, k = 15)

The combined (k = 1 + k = 15) amplitude is found from the two separate runs as
a root of the sum of squares; see amg-convergence.py for the full justification.

Each call to mf6.solve() performs one outer Picard iteration containing exactly
one inner CG step (INNER_MAXIMUM 1).  With RELAXATION_FACTOR 0 and
LINEAR_ACCELERATION CG the preconditioner applies an A-norm-optimal step
length, so the per-iteration amplitude reduction directly measures how
effectively each preconditioner attenuates the corresponding error mode.

Figure panels
-------------
  A  ILU(0) preconditioner: error amplitude vs. outer iteration for
     k=1, k=15, and the derived k=1+k=15 combined amplitude; legend shown
  B  AMG preconditioner (CG, adaptive omega): same quantities as A

Output: doc/SuppTechInfo/Figures/amg-convergence-heteroK.pdf

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

N1 = nrow - 1  # 200
N2 = ncol - 1  # 200

KX_SMOOTH, KY_SMOOTH = 1, 1  # smooth 2-D mode (k = 1 in both directions)
KX_OSCILL, KY_OSCILL = 15, 15  # oscillatory mode (k = 15)

N_ITER = 200  # outer iterations to record per solver

# ---------------------------------------------------------------------------
# Lognormal K field (reproducible)
# ---------------------------------------------------------------------------
rng = np.random.default_rng(seed=42)
K_2d = rng.lognormal(mean=0.0, sigma=1.5, size=(nrow, ncol))  # (201, 201)

print(
    f"K field: min={K_2d.min():.3f}  max={K_2d.max():.3f}  "
    f"geometric mean={np.exp(np.log(K_2d).mean()):.3f}  "
    f"sigma_log={np.log(K_2d).std():.3f}"
)

# ---------------------------------------------------------------------------
# Mode arrays and initial conditions
# ---------------------------------------------------------------------------
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


def build_model(ws, name, use_amg, ic, k):
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
    flopy.mf6.ModflowGwfnpf(gwf, k=k, icelltype=0)

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
        build_model(ws, name, use_amg, ic, K_2d)
        print(f"  Running {name} ({nrow}x{ncol}, up to {N_ITER} iters)...")
        histories[name] = run_and_collect(ws, name)


# Compute per-iteration mode amplitudes.
# The PDE is linear and the two modes are orthogonal eigenvectors of the
# uniform-K Laplacian; for heterogeneous K they are only approximately
# orthogonal, but the root-sum-square combination remains a useful summary.
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

# ---------------------------------------------------------------------------
# Figure — 2-panel layout: A = ILU(0), B = AMG
# ---------------------------------------------------------------------------

width = 6.8  # inches, USGS two-column width
dpi = 300

with fstyles.USGSPlot():
    fig, axes = plt.subplot_mosaic(
        [["A", "A", "B", "B"]],
        figsize=(width, width * 0.31),
        layout="constrained",
    )
    axes["B"].sharey(axes["A"])

    fstyles.heading(axes["A"], idx=0)
    fstyles.heading(axes["B"], idx=1)

    # y-axis limits (shared)
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
    yhi = 2.0

    # --- Panel A: ILU(0) error-component convergence ---
    ax = axes["A"]
    (h_both,) = ax.semilogy(
        ilu_iters,
        ilu_both_amp,
        color=C_BOTH,
        ls="-",
        lw=3.0,
        label=rf"$k=1 + k={KX_OSCILL}$ (smooth + oscillatory)",
    )
    (h_k15,) = ax.semilogy(
        ilu_iters,
        ilu_k15_amp,
        color=C_OSCILL,
        ls=":",
        lw=1.2,
        label=rf"$k={KX_OSCILL}$ (oscillatory)",
    )
    (h_k1,) = ax.semilogy(
        ilu_iters, ilu_k1_amp, color=C_SMOOTH, ls="-", lw=1.2, label=r"$k=1$ (smooth)"
    )
    ax.set_xlabel("Outer iteration")
    ax.set_ylabel("Error, in meters")
    ax.set_title("ILU(0)")
    ax.set_xlim([0, xmax])
    ax.set_ylim([ylo, yhi])
    ax.legend(
        handles=[h_k1, h_k15, h_both],
        fontsize=6.5,
        loc="lower left",
        ncol=1,
        framealpha=0.9,
    )

    # --- Panel B: AMG error-component convergence ---
    ax = axes["B"]
    ax.semilogy(amg_iters, amg_both_amp, color=C_BOTH, ls="-", lw=3.0)
    ax.semilogy(amg_iters, amg_k15_amp, color=C_OSCILL, ls=":", lw=1.2)
    ax.semilogy(amg_iters, amg_k1_amp, color=C_SMOOTH, ls="-", lw=1.2)
    ax.set_xlabel("Outer iteration")
    ax.set_title(r"AMG (CG, adaptive $\omega$)")
    ax.set_xlim([0, xmax])

    outfile = figpth / "amg-convergence-heteroK.pdf"
    fig.savefig(outfile, dpi=dpi)
    print(f"\nSaved {outfile}")
