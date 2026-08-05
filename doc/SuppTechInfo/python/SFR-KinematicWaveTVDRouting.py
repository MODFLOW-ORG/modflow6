"""
Figure showing that the TVD kinematic-wave scheme (ATS_COURANT) closely
matches the exact downstream hydrograph at Courant numbers where the
standard scheme introduces substantial numerical diffusion.

A two-reach MODFLOW 6 model is used so that the downstream reach has a
single upstream neighbor and the van Leer anti-diffusion correction is
active.  Results for the downstream reach are compared against the exact
kinematic-wave solution (pure translation delayed by two travel times).

Two figures are saved:
  SFRKinematicWaveTVDRouting.pdf -- hydrograph comparison at Cr=0.5 and 0.75
  SFRKinematicWaveTVDBudget.pdf  -- per-step budget discrepancy at Cr=0.5 and 0.75

Channel parameters match SFR-KinematicWave.py (Example 1):
  L = 500 m per reach, B = 10 m, S0 = 0.001, n = 0.03, Q0 = 5 m^3/s per reach.
"""

import re
import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Channel parameters (match SFR-KinematicWave.py, Example 1)
# ---------------------------------------------------------------------------
N_MANN = 0.03
B = 10.0
S0 = 0.001
L_REACH = 500.0
Q0 = 5.0
DEPS = 1.0e-5

MF6_EXE = Path(__file__).resolve().parent.parent.parent.parent / "bin" / "mf6"

# ---------------------------------------------------------------------------
# Hydraulics (rectangular channel)
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


def _area_from_Q(Q):
    return B * max(_depth_from_Q(Q), 0.0)


def _celerity_from_Q(Q):
    if Q <= 0.0:
        return 0.0
    da = _area_from_Q(Q + DEPS) - _area_from_Q(Q)
    return DEPS / da if da > 1.0e-20 else 0.0


# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
c0 = _celerity_from_Q(Q0)
T_travel = L_REACH / c0
T_wave = 10.0 * T_travel
Q_amp = 0.4 * Q0
t_end = T_wave + 2.0 * T_travel  # full wave plus two travel times


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

NWARM = 20  # warmup periods at Q0 before sinusoidal hydrograph


def _build_run(name, ws, use_tvd, Cr):
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)

    dt = Cr * T_travel
    n_main = int(np.ceil(t_end / dt))
    nper_total = NWARM + n_main

    rtp = 20.0
    d0 = _depth_from_Q(Q0)
    stage0 = rtp + d0

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

    sfr_kwargs = dict(
        print_flows=True,
        storage=True,
        nreaches=2,
        packagedata=pak_data,
        connectiondata=_CONNECTION_DATA,
        perioddata=sfr_perioddata,
        initialstages=initial_stages,
        pname="sfr-1",
    )
    if use_tvd:
        sfr_kwargs["ats_courant"] = 1.0
        sfr_kwargs["maximum_depth_change"] = 1.0e-9

    sfr = flopy.mf6.ModflowGwfsfr(gwf, **sfr_kwargs)

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

    listing = (ws / f"{name}.lst").read_text()
    sections = re.split(r"SFR-1 BUDGET FOR ENTIRE MODEL", listing)[1:]
    disc_all = []
    for sec in sections:
        m = re.search(r"PERCENT DISCREPANCY\s*=\s*([0-9.Ee+\-]+)", sec)
        if m:
            disc_all.append(abs(float(m.group(1))))
    disc = np.asarray(disc_all[skip:])

    return t_s, q_out, disc


# ---------------------------------------------------------------------------
# Run all four cases
# ---------------------------------------------------------------------------
Cr_values = [0.5, 0.75]
std_results = {}
tvd_results = {}
std_budget = {}
tvd_budget = {}

for Cr in Cr_values:
    with tempfile.TemporaryDirectory() as _tmpdir:
        _tmp = Path(_tmpdir)
        t_std, q_std, d_std = _build_run(
            f"sfr-tvd-std-cr{int(Cr * 100)}", _tmp / "std", False, Cr
        )
        t_tvd, q_tvd, d_tvd = _build_run(
            f"sfr-tvd-tvd-cr{int(Cr * 100)}", _tmp / "tvd", True, Cr
        )
    std_results[Cr] = (t_std / T_travel, q_std / Q0)
    tvd_results[Cr] = (t_tvd / T_travel, q_tvd / Q0)
    std_budget[Cr] = (t_std / T_travel, d_std)
    tvd_budget[Cr] = (t_tvd / T_travel, d_tvd)

# Dense curves for exact solution and upstream inflow
t_dense = np.linspace(0.0, t_end, 4000)
Q_up_dense = Q_upstream(t_dense)
Q_ex_dense = np.where(
    t_dense >= 2.0 * T_travel,
    Q_upstream(t_dense - 2.0 * T_travel),
    Q0,
)
t_norm_dense = t_dense / T_travel

EPS_FLOOR = 1.0e-15
figpth = "../Figures"
panel_labels = ["A", "B"]

# ---------------------------------------------------------------------------
# Figure 1: Routing hydrographs
# ---------------------------------------------------------------------------
with styles.USGSPlot():
    fig, axes = plt.subplots(
        1, 2, figsize=(6.8, 3.2), constrained_layout=True, sharey=True
    )

    for ax, Cr, letter in zip(axes, Cr_values, panel_labels):
        t_norm_std, Q_std_norm = std_results[Cr]
        t_norm_tvd, Q_tvd_norm = tvd_results[Cr]

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
        ax.plot(
            t_norm_std,
            Q_std_norm,
            color="C1",
            lw=1.2,
            ls="-",
            label=rf"Standard ($C_r = {Cr}$)",
            zorder=4,
        )
        ax.plot(
            t_norm_tvd,
            Q_tvd_norm,
            color="C0",
            lw=1.2,
            ls="-",
            label=rf"TVD ($C_r = {Cr}$)",
            zorder=3,
        )

        ax.set_xlim(0, t_end / T_travel)
        ax.set_ylim(0.20, 1.85)
        ax.set_xlabel(r"Normalized time, $t \, / \, T$")
        ax.axhline(1.0, color="0.75", lw=0.5, zorder=1)
        ax.tick_params(direction="in", top=True, right=True)
        styles.heading(ax=ax, letter=letter)
        styles.graph_legend(
            ax=ax,
            ncols=2,
            loc="lower left",
            fontsize=7,
            handlelength=1.4,
            borderpad=0.5,
            framealpha=0.9,
            edgecolor="0.7",
        )

    axes[0].set_ylabel(r"Normalized discharge, $Q \, / \, Q_0$")

fig.savefig(f"{figpth}/SFRKinematicWaveTVDRouting.pdf", dpi=300)
print(f"Saved {figpth}/SFRKinematicWaveTVDRouting.pdf")

# ---------------------------------------------------------------------------
# Figure 2: Budget discrepancy
# ---------------------------------------------------------------------------
with styles.USGSPlot():
    fig2, axes2 = plt.subplots(
        1, 2, figsize=(6.8, 3.2), constrained_layout=True, sharey=True
    )
    t_end_norm = t_end / T_travel

    for ax, Cr, letter in zip(axes2, Cr_values, panel_labels):
        t_std_b, d_std_b = std_budget[Cr]
        t_tvd_b, d_tvd_b = tvd_budget[Cr]

        d_std_plot = np.where(d_std_b > EPS_FLOOR, d_std_b, EPS_FLOOR)
        d_tvd_plot = np.where(d_tvd_b > EPS_FLOOR, d_tvd_b, EPS_FLOOR)

        if len(t_std_b) > 0:
            ax.plot(
                t_std_b,
                d_std_plot,
                color="C1",
                lw=1.0,
                ls="-",
                label=rf"Standard ($C_r = {Cr}$)",
                zorder=4,
            )
        if len(t_tvd_b) > 0:
            ax.plot(
                t_tvd_b,
                d_tvd_plot,
                color="C0",
                lw=1.5,
                ls="-",
                label=rf"TVD ($C_r = {Cr}$)",
                zorder=5,
            )
            i_ann = len(t_tvd_b) // 3
            ax.annotate(
                "TVD: 0.00%\nevery step",
                xy=(t_tvd_b[i_ann], EPS_FLOOR),
                xytext=(t_tvd_b[i_ann], 1.0e-9),
                ha="center",
                va="bottom",
                fontsize=5.5,
                color="C0",
                arrowprops=dict(arrowstyle="-", color="C0", lw=0.6),
            )

        ax.axhline(1.0, color="0.45", lw=0.7, ls="--", zorder=1)
        ax.text(
            t_end_norm * 0.97,
            1.0,
            "1%",
            ha="right",
            va="bottom",
            fontsize=6.5,
            color="0.4",
        )
        ax.set_yscale("log")
        ax.set_ylim(1.0e-16, 100.0)
        ax.set_xlim(0, t_end_norm)
        ax.set_xlabel(r"Normalized time, $t \, / \, T$")
        ax.tick_params(direction="in", top=True, right=True, which="both")
        styles.heading(ax=ax, letter=letter)
        styles.graph_legend(
            ax=ax,
            ncols=1,
            loc="center right",
            fontsize=7,
            handlelength=1.4,
            borderpad=0.5,
            framealpha=0.9,
            edgecolor="0.7",
        )

    axes2[0].set_ylabel("Percent Budget Discrepancy (%)")

fig2.savefig(f"{figpth}/SFRKinematicWaveTVDBudget.pdf", dpi=300)
print(f"Saved {figpth}/SFRKinematicWaveTVDBudget.pdf")

for Cr in Cr_values:
    _, d_std = std_budget[Cr]
    _, d_tvd = tvd_budget[Cr]
    if len(d_std) > 0:
        print(f"Cr={Cr}: Standard max disc: {d_std.max():.3e}%")
    if len(d_tvd) > 0:
        print(f"Cr={Cr}: TVD      max disc: {d_tvd.max():.3e}%")

print(f"c0={c0:.4f} m/s  T_travel={T_travel:.2f} s")
