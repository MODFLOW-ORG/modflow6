"""
Space-time grid diagram for the kinematic-wave four-point finite-difference
stencil used in the MODFLOW 6 SFR Package.

Variable naming follows GwfSfrCommon.f90:
  qa = Q_U^{n-1}   upstream discharge, previous time step
  qb = Q_D^{n-1}   downstream discharge, previous time step
  qc = Q_U^n       upstream discharge, current time step  (known inflow)
  qd = Q_D^n       downstream discharge, current time step (unknown)

Panel A: Cr <= 1  -- four-point weighted stencil
Panel B: Cr >  1  -- upstream-only storage stencil

Node labels are centered directly above (qc, qd) or below (qa, qb) each
node so that they are unambiguously aligned with their node and do not
compete with the Δt annotation on the left margin.
"""

import flopy.plot.styles as styles
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------
X_U, X_D = 0.0, 1.0  # upstream / downstream x positions
T_OLD, T_NEW = 0.0, 1.0  # previous / current time-step positions

NODE_KW = dict(s=80, zorder=5)
KNOWN_COLOR = "0.3"  # filled nodes (known values)
UNKNOWN_COLOR = "C3"  # highlighted node (unknown qd)
INACTIVE_COLOR = "0.75"  # nodes not used in Cr > 1 stencil

LABEL_FS = 8.5  # node-label font size
ANNOT_FS = 7.5  # Δx / Δt annotation font size
EDGE_LW = 1.2
LBL_OFFSET = 0.12  # vertical distance from node centre to label anchor

# ---------------------------------------------------------------------------
# Node positions and labels
# ---------------------------------------------------------------------------
nodes = {
    "qa": (X_U, T_OLD),
    "qb": (X_D, T_OLD),
    "qc": (X_U, T_NEW),
    "qd": (X_D, T_NEW),
}

labels = {
    "qa": r"$Q_a = Q_U^{\,n-1}$",
    "qb": r"$Q_b = Q_D^{\,n-1}$",
    "qc": r"$Q_c = Q_U^{\,n}$",
    "qd": r"$Q_d = Q_D^{\,n}$",
}


def draw_panel(ax, cr_case, title_letter):
    """Draw one panel for the given Cr case."""

    # -----------------------------------------------------------------------
    # Grid lines (light)
    # -----------------------------------------------------------------------
    for x in [X_U, X_D]:
        ax.axvline(x, color="0.82", lw=0.6, zorder=0)
    for t in [T_OLD, T_NEW]:
        ax.axhline(t, color="0.82", lw=0.6, zorder=0)

    # -----------------------------------------------------------------------
    # Stencil edges
    # -----------------------------------------------------------------------
    if cr_case == "low":
        edges = [("qa", "qb"), ("qc", "qd"), ("qa", "qc"), ("qb", "qd")]
        active_nodes = {"qa", "qb", "qc", "qd"}
    else:
        edges = [("qc", "qd"), ("qa", "qc")]
        active_nodes = {"qa", "qc", "qd"}

    for n1, n2 in edges:
        x1, t1 = nodes[n1]
        x2, t2 = nodes[n2]
        ax.plot(
            [x1, x2],
            [t1, t2],
            color="0.35",
            lw=EDGE_LW,
            solid_capstyle="round",
            zorder=2,
        )

    # -----------------------------------------------------------------------
    # Nodes
    # -----------------------------------------------------------------------
    for key, (x, t) in nodes.items():
        if key == "qd":
            ax.scatter(
                x,
                t,
                s=NODE_KW["s"],
                facecolors=UNKNOWN_COLOR,
                edgecolors=UNKNOWN_COLOR,
                zorder=NODE_KW["zorder"],
                linewidths=1.5,
            )
        elif key in active_nodes:
            ax.scatter(
                x,
                t,
                s=NODE_KW["s"],
                facecolors=KNOWN_COLOR,
                edgecolors=KNOWN_COLOR,
                zorder=NODE_KW["zorder"],
            )
        else:
            ax.scatter(
                x,
                t,
                s=NODE_KW["s"],
                facecolors="white",
                edgecolors=INACTIVE_COLOR,
                zorder=NODE_KW["zorder"],
                linewidths=1.2,
            )

    # -----------------------------------------------------------------------
    # Node labels — centered above (qc, qd) or below (qa, qb) each node
    # -----------------------------------------------------------------------
    for key, (x, t) in nodes.items():
        txt = labels[key]
        col = (
            UNKNOWN_COLOR
            if key == "qd"
            else (
                INACTIVE_COLOR if (cr_case == "high" and key == "qb") else KNOWN_COLOR
            )
        )
        if t == T_OLD:  # bottom row: label below the node
            ax.text(
                x,
                t - LBL_OFFSET,
                txt,
                ha="center",
                va="top",
                fontsize=LABEL_FS,
                color=col,
            )
        else:  # top row: label above the node
            ax.text(
                x,
                t + LBL_OFFSET,
                txt,
                ha="center",
                va="bottom",
                fontsize=LABEL_FS,
                color=col,
            )

    # -----------------------------------------------------------------------
    # Δx and Δt annotations
    # -----------------------------------------------------------------------
    # Δx: double-headed arrow below the node labels (further down than before)
    y_arr = T_OLD - 0.32
    ax.annotate(
        "",
        xy=(X_D, y_arr),
        xytext=(X_U, y_arr),
        arrowprops=dict(arrowstyle="<->", color="0.4", lw=0.9),
    )
    ax.text(
        0.5 * (X_U + X_D),
        y_arr - 0.07,
        r"$\Delta x$",
        ha="center",
        va="top",
        fontsize=ANNOT_FS,
        color="0.4",
    )

    # Δt: double-headed arrow to the left of the upstream column
    x_arr = X_U - 0.28
    ax.annotate(
        "",
        xy=(x_arr, T_NEW),
        xytext=(x_arr, T_OLD),
        arrowprops=dict(arrowstyle="<->", color="0.4", lw=0.9),
    )
    ax.text(
        x_arr - 0.05,
        0.5 * (T_OLD + T_NEW),
        r"$\Delta t$",
        ha="right",
        va="center",
        fontsize=ANNOT_FS,
        color="0.4",
    )

    # -----------------------------------------------------------------------
    # Axis labels and limits
    # -----------------------------------------------------------------------
    ax.set_xlim(X_U - 0.48, X_D + 0.48)
    ax.set_ylim(T_OLD - 0.58, T_NEW + 0.35)
    ax.set_xticks([X_U, X_D])
    ax.set_xticklabels([r"$x_U$", r"$x_D$"], fontsize=LABEL_FS)
    ax.set_yticks([T_OLD, T_NEW])
    ax.set_yticklabels([r"$t^{\,n-1}$", r"$t^{\,n}$"], fontsize=LABEL_FS)
    ax.tick_params(direction="in", top=False, right=False, length=0)
    ax.set_xlabel(r"Distance, $x$", labelpad=4)
    ax.set_ylabel(r"Time, $t$", labelpad=4)

    # -----------------------------------------------------------------------
    # Subtitle (Cr case)
    # -----------------------------------------------------------------------
    subtitle = (
        r"$C_r \leq 1$: four-point weighted stencil"
        if cr_case == "low"
        else r"$C_r > 1$: upstream-only storage stencil"
    )
    ax.set_title(subtitle, fontsize=7.5, pad=4, color="0.3")

    styles.heading(ax=ax, letter=title_letter)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
with styles.USGSPlot():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.4), constrained_layout=True)
    draw_panel(axes[0], "low", "A")
    draw_panel(axes[1], "high", "B")

    handles = [
        mpatches.Patch(color=KNOWN_COLOR, label="Known discharge"),
        mpatches.Patch(color=UNKNOWN_COLOR, label=r"Unknown $Q_d = Q_D^n$"),
        mpatches.Patch(color=INACTIVE_COLOR, label="Not used in stencil"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncols=3,
        fontsize=7,
        bbox_to_anchor=(0.5, -0.08),
        framealpha=0.9,
        edgecolor="0.7",
    )

figpth = "../Figures"
fig.savefig(f"{figpth}/SFRKinematicWaveStencil.pdf", dpi=300, bbox_inches="tight")
print(f"Saved {figpth}/SFRKinematicWaveStencil.pdf")
