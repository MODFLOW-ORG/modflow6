"""
Figure comparing MODFLOW 6 kinematic-wave finite-difference scheme
against the exact kinematic-wave solution for a single rectangular-
channel reach at Courant numbers ranging from 0.5 to 2.0.

Governing equations match kinematic_residual / kinematic_storage in
src/Model/ModelUtilities/GwfSfrCommon.f90.
"""

import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Channel and reach parameters
# ---------------------------------------------------------------------------
ALPHA = 1.0  # SI unit conversion (m^(1/3)/s)
N_MANN = 0.03  # Manning roughness (s/m^(1/3))
B = 10.0  # channel width (m)
S0 = 0.001  # bed slope (m/m)
L = 500.0  # reach length (m)
Q0 = 5.0  # base discharge (m^3/s)
THETA = 1.0  # temporal weighting (default MF6 value)

# ---------------------------------------------------------------------------
# Manning's-equation utilities for a rectangular channel
# ---------------------------------------------------------------------------


def Q_from_depth(d):
    if d <= 0.0:
        return 0.0
    A = B * d
    P = B + 2.0 * d
    R = A / P
    return (ALPHA / N_MANN) * A * R ** (2.0 / 3.0) * S0**0.5


def depth_from_Q(Q):
    if Q <= 0.0:
        return 0.0
    d_max = 50.0
    while Q_from_depth(d_max) < Q:
        d_max *= 2.0
    return brentq(lambda d: Q_from_depth(d) - Q, 1e-12, d_max)


def area_from_Q(Q):
    return B * depth_from_Q(Q)


def celerity_from_Q(Q, dq=1e-4):
    if Q <= 0.0:
        return 0.0
    dA = area_from_Q(Q + dq) - area_from_Q(Q)
    return dq / dA if dA > 0.0 else 0.0


# ---------------------------------------------------------------------------
# MODFLOW 6 kinematic-wave residual (GwfSfrCommon.f90)
# ---------------------------------------------------------------------------


def kinematic_residual(
    Qa, Qb, Qc, Qd, Aa, Ab, Ac, Ad, qsrc, length, theta, delt, courant
):
    if courant > 1.0:
        f11 = (Qd - Qc) / length
        f12 = (Ac - Aa) / delt
    else:
        f11 = (theta * (Qd - Qc) + (1.0 - theta) * (Qb - Qa)) / length
        f12 = 0.5 * ((Ad - Ab) + (Ac - Aa)) / delt
    return f11 + f12 - qsrc


# ---------------------------------------------------------------------------
# Newton-Raphson solver (sfr_calc_transient)
# ---------------------------------------------------------------------------


def solve_reach(Qa, Qb, Qc, dt, courant, qsrc=0.0, max_iter=50, dq=1e-4, tol=1e-10):
    Aa = area_from_Q(Qa)
    Ab = area_from_Q(Qb)
    Ac = area_from_Q(Qc)
    Qd = max(Qb, 0.0)
    for _ in range(max_iter):
        Ad = area_from_Q(Qd)
        res = kinematic_residual(
            Qa, Qb, Qc, Qd, Aa, Ab, Ac, Ad, qsrc, L, THETA, dt, courant
        )
        Ad2 = area_from_Q(Qd + dq)
        res2 = kinematic_residual(
            Qa, Qb, Qc, Qd + dq, Aa, Ab, Ac, Ad2, qsrc, L, THETA, dt, courant
        )
        deriv = (res2 - res) / dq
        if abs(deriv) < 1e-30:
            break
        delta = -res / deriv
        Qd = max(Qd + delta, 0.0)
        if abs(delta) < tol and abs(res) < tol:
            break
    return Qd


# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------

c0 = celerity_from_Q(Q0)
T_travel = L / c0
T_wave = 10.0 * T_travel
Q_amp = 0.4 * Q0


def Q_upstream(t):
    return Q0 + Q_amp * np.sin(2.0 * np.pi * t / T_wave)


def Q_exact(t):
    return Q_upstream(t - T_travel)


Cr_values = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
t_end = 2.0 * T_wave + T_travel

# ---------------------------------------------------------------------------
# Run MF6 FD scheme
# ---------------------------------------------------------------------------

results = {}
for Cr in Cr_values:
    dt = Cr * T_travel
    n_steps = int(np.ceil(t_end / dt)) + 1
    times = np.arange(n_steps) * dt
    Qd_arr = np.zeros(n_steps)
    Qd_arr[0] = Q0
    Qa_prev = Q0
    Qb_prev = Q0
    for i in range(1, n_steps):
        Qc = Q_upstream(times[i])
        Qd_arr[i] = solve_reach(Qa_prev, Qb_prev, Qc, dt, Cr)
        Qa_prev = Qc
        Qb_prev = Qd_arr[i]
    results[Cr] = (times, Qd_arr)

# ---------------------------------------------------------------------------
# Reference curves
# ---------------------------------------------------------------------------

t_dense = np.linspace(0.0, t_end, 2000)
Q_up_dense = Q_upstream(t_dense)
Q_ex_dense = np.where(t_dense >= T_travel, Q_exact(t_dense), Q0)
t_norm_dense = t_dense / T_travel

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

colors_low = plt.cm.Blues_r(np.linspace(0.10, 0.65, 3))
colors_high = plt.cm.Reds(np.linspace(0.35, 0.85, 3))

panel_cfg = [
    ([0.5, 0.75, 1.0], colors_low, "A"),
    ([1.25, 1.5, 2.0], colors_high, "B"),
]

with styles.USGSPlot():
    fig, axes = plt.subplots(
        1, 2, figsize=(6.8, 3.2), constrained_layout=True, sharey=True
    )

    for ax, (ax_Cr_values, ax_colors, panel_letter) in zip(axes, panel_cfg):
        ax.plot(
            t_norm_dense,
            Q_up_dense / Q0,
            color="black",
            lw=0.8,
            ls="--",
            label=r"$Q_U$ (upstream)",
            zorder=5,
        )
        ax.plot(
            t_norm_dense,
            Q_ex_dense / Q0,
            color="black",
            lw=1.5,
            ls="-",
            label="Exact",
            zorder=5,
        )
        for Cr, color in zip(ax_Cr_values, ax_colors):
            t_sim, Qd_sim = results[Cr]
            ax.plot(
                t_sim / T_travel,
                Qd_sim / Q0,
                color=color,
                lw=1.0,
                ls="-",
                label=rf"$C_r = {Cr}$",
                zorder=4,
            )

        ax.set_xlim(0, t_end / T_travel)
        ax.set_ylim(0.20, 1.85)
        ax.set_xlabel(r"Normalized time, $t \, / \, T$")
        ax.tick_params(direction="in", top=True, right=True)
        ax.axhline(1.0, color="0.75", lw=0.5, zorder=1)

        styles.heading(ax=ax, letter=panel_letter)
        styles.graph_legend(
            ax=ax,
            ncols=3,
            loc="lower left",
            fontsize=7,
            handlelength=1.4,
            borderpad=0.5,
            framealpha=0.9,
            edgecolor="0.7",
        )

    axes[0].set_ylabel(r"Normalized discharge, $Q \, / \, Q_0$")

figpth = "../Figures"
fig.savefig(f"{figpth}/SFRKinematicWave.pdf", dpi=300)
print(f"Saved {figpth}/SFRKinematicWave.pdf")
print(f"B={B} m  S={S0}  n={N_MANN}  L={L} m  Q0={Q0} m^3/s")
print(f"c0={c0:.4f} m/s  T_travel={T_travel:.2f} s  T_wave={T_wave:.2f} s")
