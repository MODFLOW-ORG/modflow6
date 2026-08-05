"""
Figure comparing MODFLOW 6 kinematic-wave finite-difference scheme
against the exact kinematic-wave solution at Courant numbers ranging
from 0.5 to 2.0.

A two-reach MODFLOW 6 model is run at each Courant number with a
sinusoidal upstream boundary condition.  The downstream outflow (reach 1)
is compared against the exact kinematic-wave solution, which is a pure
translation delayed by two travel times (2T).

Channel: L = 500 m per reach, B = 10 m, S0 = 0.001, n = 0.03, Q0 = 5 m^3/s.
"""

import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Channel and reach parameters
# ---------------------------------------------------------------------------
N_MANN = 0.03
B = 10.0
S0 = 0.001
L_REACH = 500.0
Q0 = 5.0
DEPS = 1.0e-5

MF6_EXE = Path(__file__).resolve().parent.parent.parent.parent / "bin" / "mf6"

# ---------------------------------------------------------------------------
# Manning's-equation utilities for a rectangular channel
# ---------------------------------------------------------------------------


def _Q_from_depth(d):
    if d <= 0.0:
        return 0.0
    A = B * d
    P = B  # MF6 uses a wetted perimeter of B, so R = A / P = d
    return (1.0 / N_MANN) * A * (A / P) ** (2.0 / 3.0) * S0**0.5


def _depth_from_Q(Q):
    if Q <= 0.0:
        return 0.0
    d_max = 50.0
    while _Q_from_depth(d_max) < Q:
        d_max *= 2.0
    return brentq(lambda d: _Q_from_depth(d) - Q, 1.0e-12, d_max)


# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
c0 = DEPS / (B * max(_depth_from_Q(Q0 + DEPS) - _depth_from_Q(Q0), 1.0e-20))
T_travel = L_REACH / c0
T_wave = 10.0 * T_travel
Q_amp = 0.4 * Q0
t_end = 2.0 * T_wave + 2.0 * T_travel  # two full waves plus two travel-time delays


def Q_upstream(t):
    return Q0 + Q_amp * np.sin(2.0 * np.pi * t / T_wave)


# ---------------------------------------------------------------------------
# Two-reach MODFLOW 6 model
# ---------------------------------------------------------------------------
_CONNECTION_DATA = [
    [0, -1],
    [1, 0],
]
_NCON = [len(row) - 1 for row in _CONNECTION_DATA]

NWARM = 20  # warmup periods at Q0; initialises usinflowold before sinusoid


def _build_run(name, ws, Cr):
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)

    dt = Cr * T_travel
    n_main = int(np.ceil(t_end / dt))
    nper_total = NWARM + n_main

    rtp = 20.0
    stage0 = rtp + _depth_from_Q(Q0)

    pak_data = [
        (irno, -1, -1, -1, L_REACH, B, S0, rtp, 1.0, 0.0, N_MANN, nc, 1.0, 0)
        for irno, nc in enumerate(_NCON)
    ]
    initial_stages = [(irno, stage0) for irno in range(2)]

    # Warmup at constant Q0; main hydrograph uses new-time inflow values
    sfr_perioddata = {i: [(0, "inflow", Q0)] for i in range(NWARM)}
    sfr_perioddata.update(
        {
            NWARM + i: [(0, "inflow", float(Q_upstream((i + 1) * dt)))]
            for i in range(n_main)
        }
    )

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=str(ws),
        version="mf6",
        exe_name=str(MF6_EXE),
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=nper_total,
        perioddata=[(dt, 1, 1.0)] * nper_total,
    )
    flopy.mf6.ModflowIms(sim, outer_dvclose=1.0e-8, inner_dvclose=1.0e-9)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=1,
        nrow=1,
        ncol=1,
        delr=L_REACH,
        delc=L_REACH,
        top=20.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0e-4)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1.0e-6, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_flows=True,
        storage=True,
        nreaches=2,
        packagedata=pak_data,
        connectiondata=_CONNECTION_DATA,
        perioddata=sfr_perioddata,
        initialstages=initial_stages,
        pname="sfr-1",
    )

    obs_fname = f"{name}.sfr.obs"
    sfr.obs.initialize(
        filename=obs_fname,
        print_input=True,
        continuous={
            f"{obs_fname}.csv": [
                ("inflow", "ext-inflow", (0,)),
                ("outflow", "ext-outflow", (1,)),
            ]
        },
    )

    sim.write_simulation()
    success, _ = sim.run_simulation(silent=True)
    assert success, f"MODFLOW 6 run failed for {name}"

    obs = np.genfromtxt(ws / f"{obs_fname}.csv", delimiter=",", names=True)
    t_all = obs[obs.dtype.names[0]]
    q_out_all = -obs["OUTFLOW"]

    # Drop warmup periods; reset time origin to start of main hydrograph
    skip = NWARM
    t_s = t_all[skip:] - t_all[skip - 1]
    q_out = q_out_all[skip:]

    # Prepend initial condition at t = 0
    t_s = np.concatenate([[0.0], t_s])
    q_out = np.concatenate([[Q0], q_out])

    return t_s, q_out


# ---------------------------------------------------------------------------
# Run all six Courant numbers
# ---------------------------------------------------------------------------
Cr_values = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
results = {}

for Cr in Cr_values:
    cr_tag = str(Cr).replace(".", "p")
    with tempfile.TemporaryDirectory() as _tmpdir:
        t_s, q_out = _build_run(f"sfr-kw-cr{cr_tag}", Path(_tmpdir), Cr)
    results[Cr] = (t_s / T_travel, q_out / Q0)
    print(f"Cr={Cr}: done ({len(t_s) - 1} main steps)")

# ---------------------------------------------------------------------------
# Reference curves
# ---------------------------------------------------------------------------
t_dense = np.linspace(0.0, t_end, 2000)
Q_up_dense = Q_upstream(t_dense)
Q_ex_dense = np.where(
    t_dense >= 2.0 * T_travel,
    Q_upstream(t_dense - 2.0 * T_travel),
    Q0,
)
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
                t_sim,
                Qd_sim,
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
print(f"B={B} m  S={S0}  n={N_MANN}  L_REACH={L_REACH} m  Q0={Q0} m^3/s")
print(f"c0={c0:.4f} m/s  T_travel={T_travel:.2f} s  T_wave={T_wave:.2f} s")
