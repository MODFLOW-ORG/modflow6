"""
Figure showing kinematic-wave Courant number and optimal time step for
a rectangular channel as a function of discharge and bed slope.

Panel A: Courant number Cr = c*dt/dx vs Q for five bed slopes,
         evaluated at a reference time step of 10 min (600 s) and
         reach length of 500 m.  Gray band shows recommended range
         0.5 <= Cr <= 2.  For other choices, Cr scales as
         (dt / 10 min) * (500 m / dx).

Panel B: Optimal time step Delta_t_opt = dx / c (minutes) vs Q for
         the same five slopes (dx = 500 m).  Horizontal dashed lines
         mark common model time steps (5, 15, and 60 min).  Where a
         colored curve lies above a reference line, that time step
         gives Cr < 1; where below, Cr > 1.  The acceptable range
         is roughly half to twice the plotted curve value
         (0.5 <= Cr <= 2).

Channel parameters are fixed at B = 10 m and n = 0.035.
"""

import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Channel parameters (fixed)
# ---------------------------------------------------------------------------
ALPHA = 1.0
N = 0.035
B = 10.0
DT_REF = 600.0  # reference time step for Cr panel: 10 min = 600 s
DX_REF = 500.0  # reference reach length: 500 m

# ---------------------------------------------------------------------------
# Hydraulic utilities
# ---------------------------------------------------------------------------


def Q_from_depth(d, S0):
    if d <= 0.0:
        return 0.0
    A = B * d
    P = B  # MF6 uses a wetted perimeter of B, so R = A / P = d
    R = A / P
    return (ALPHA / N) * A * R ** (2.0 / 3.0) * S0**0.5


def depth_from_Q(Q, S0):
    if Q <= 0.0:
        return 0.0
    d_max = 200.0
    while Q_from_depth(d_max, S0) < Q:
        d_max *= 2.0
    return brentq(lambda d: Q_from_depth(d, S0) - Q, 1e-12, d_max)


def celerity_from_Q(Q, S0, dq=1e-5):
    if Q <= 0.0:
        return 0.0
    d1 = depth_from_Q(Q, S0)
    d2 = depth_from_Q(Q + dq, S0)
    dA = B * (d2 - d1)
    return dq / dA if dA > 0.0 else 0.0


# ---------------------------------------------------------------------------
# Slopes and discharge range
# ---------------------------------------------------------------------------
S0_values = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
S0_labels = [
    r"$S_0 = 0.0001$",
    r"$S_0 = 0.0003$",
    r"$S_0 = 0.001$",
    r"$S_0 = 0.003$",
    r"$S_0 = 0.01$",
]

Q_arr = np.logspace(-0.5, 3.0, 300)  # ~0.3 to 1000 m³/s

celerity = {}
courant = {}
dt_opt = {}
for S0 in S0_values:
    c_arr = np.array([celerity_from_Q(Q, S0) for Q in Q_arr])
    celerity[S0] = c_arr
    courant[S0] = c_arr * DT_REF / DX_REF
    dt_opt[S0] = DX_REF / c_arr / 60.0  # minutes

# ---------------------------------------------------------------------------
# Tick formatter: plain integers with comma separator for >= 1000
# ---------------------------------------------------------------------------


def _log_label(x, pos):
    if x >= 1000:
        return f"{int(x):,}"
    if x >= 1:
        return f"{int(x)}"
    return f"{x:g}"


_fmt = ticker.FuncFormatter(_log_label)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
colors = plt.cm.plasma(np.linspace(0.10, 0.85, len(S0_values)))

with styles.USGSPlot():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.2), constrained_layout=True)
    ax_cr, ax_dt = axes

    # -----------------------------------------------------------------------
    # Panel A: Courant number
    # -----------------------------------------------------------------------
    ax_cr.axhspan(0.5, 2.0, color="0.87", zorder=0, label=r"$0.5 \leq C_r \leq 2$")

    for S0, label, color in zip(S0_values, S0_labels, colors):
        ax_cr.plot(Q_arr, courant[S0], color=color, lw=1.2, label=label)

    for cr_val, ls in [(0.5, "--"), (1.0, "-"), (2.0, "--")]:
        ax_cr.axhline(cr_val, color="0.5", lw=0.6, ls=ls, zorder=1)
    ax_cr.text(
        Q_arr[-1] * 0.75,
        0.5,
        r"$C_r = 0.5$",
        va="top",
        ha="right",
        fontsize=6.5,
        color="0.4",
    )
    ax_cr.text(
        Q_arr[-1] * 0.75,
        1.0,
        r"$C_r = 1$",
        va="bottom",
        ha="right",
        fontsize=6.5,
        color="0.4",
    )
    ax_cr.text(
        Q_arr[-1] * 0.75,
        2.0,
        r"$C_r = 2$",
        va="bottom",
        ha="right",
        fontsize=6.5,
        color="0.4",
    )

    ax_cr.set_xscale("log")
    ax_cr.set_yscale("log")
    ax_cr.set_ylim(0.05, 30)
    ax_cr.set_xlabel(r"Discharge, $Q$ (m$^3$ s$^{-1}$)")
    ax_cr.set_ylabel(r"Courant number, $C_r$")
    ax_cr.xaxis.set_major_formatter(_fmt)
    ax_cr.yaxis.set_major_formatter(_fmt)
    ax_cr.tick_params(direction="in", top=True, right=True, which="both")
    styles.heading(ax=ax_cr, letter="A")
    styles.graph_legend(
        ax=ax_cr,
        loc="upper left",
        fontsize=6.5,
        handlelength=1.4,
        borderpad=0.5,
        framealpha=0.9,
        edgecolor="0.7",
    )

    # -----------------------------------------------------------------------
    # Panel B: Optimal time step
    # -----------------------------------------------------------------------
    for S0, label, color in zip(S0_values, S0_labels, colors):
        ax_dt.plot(Q_arr, dt_opt[S0], color=color, lw=1.2, label=label)

    # Horizontal reference lines for common model time steps
    ref_dts = [5, 15, 60]
    ref_labels = ["5 min", "15 min", "60 min"]
    for dt_ref, lbl in zip(ref_dts, ref_labels):
        ax_dt.axhline(dt_ref, color="0.5", lw=0.7, ls="--", zorder=1)
        ax_dt.text(
            Q_arr[0] * 1.35,
            dt_ref,
            lbl,
            va="center",
            ha="left",
            fontsize=6.5,
            color="0.4",
        )

    ax_dt.set_xscale("log")
    ax_dt.set_yscale("log")
    ax_dt.set_xlabel(r"Discharge, $Q$ (m$^3$ s$^{-1}$)")
    ax_dt.set_ylabel(r"Optimal time step, $\Delta t$ (min)")
    ax_dt.xaxis.set_major_formatter(_fmt)
    ax_dt.yaxis.set_major_formatter(_fmt)
    ax_dt.tick_params(direction="in", top=True, right=True, which="both")
    styles.heading(ax=ax_dt, letter="B")
    styles.graph_legend(
        ax=ax_dt,
        loc="upper right",
        fontsize=6.5,
        handlelength=1.4,
        borderpad=0.5,
        framealpha=0.9,
        edgecolor="0.7",
    )

figpth = "../Figures"
fig.savefig(f"{figpth}/SFRKinematicWaveCourant.pdf", dpi=300)
print(f"Saved {figpth}/SFRKinematicWaveCourant.pdf")
print(f"B={B} m  n={N}  DT_ref={DT_REF / 60:.0f} min  DX_ref={DX_REF:.0f} m")
