"""
Conceptual figure showing what happens to solute in a cell whose water table
moves, without the stranded mass option.

Panel A shows a falling water table. The sorbate of the interval that drains is
released to the water that remains, because the sorbed mass of a cell is scaled
by the saturation, and the solute in the water that is retained against drainage
is not carried into the next time step, because the mobile water volume is taken
as the saturation times the porosity, so mass is lost.

Panel B shows the water table rising again. The interval is saturated once more,
as though it had held water at the concentration of the previous time step, so
mass is created.
"""

from pathlib import Path

import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrow, Rectangle

WATER = "#9ecae1"
DRAINED = "#efe3c8"
SOLID = "0.35"
SEED = 20260812  # the scatter of sorbate is the same every time the figure is built
NDOT = 60
ARROW = "#d62728"

XL, XR, W = 0.8, 5.6, 2.6  # left cell x, right cell x, cell width
YB, YT = 1.6, 7.6  # cell bottom and top
HI, LO = 6.4, 3.4  # high and low water table


def draw_cell(ax, x, wt, label):
    """One cell with its water table, and sorbate on the solid."""
    ax.add_patch(
        Rectangle((x, YB), W, wt - YB, facecolor=WATER, edgecolor="none", zorder=2)
    )
    ax.add_patch(
        Rectangle((x, wt), W, YT - wt, facecolor=DRAINED, edgecolor="none", zorder=2)
    )
    ax.add_patch(
        Rectangle(
            (x, YB), W, YT - YB, facecolor="none", edgecolor="black", lw=1.0, zorder=4
        )
    )
    ax.plot([x, x + W], [wt, wt], color="black", lw=1.0, zorder=5)

    # sorbate held on the solid aquifer material, scattered over the whole cell
    rng = np.random.default_rng(SEED)
    xs = rng.uniform(x + 0.18, x + W - 0.18, NDOT)
    ys = rng.uniform(YB + 0.18, YT - 0.18, NDOT)
    ax.plot(xs, ys, marker="o", ms=2.2, color=SOLID, ls="none", zorder=6)

    ax.text(x + W / 2.0, YT + 0.25, label, ha="center", va="bottom", fontsize=7)


def annotate(ax, text, xy, xytext, color=ARROW):
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        fontsize=6.5,
        ha="center",
        color=color,
        arrowprops=dict(arrowstyle="->", color=color, lw=0.9, shrinkA=0, shrinkB=2),
        zorder=8,
    )


def draw_change(ax, x, y0, y1, facecolor):
    """The interval whose saturation changes over the time step."""
    ax.add_patch(
        Rectangle(
            (x, y0),
            W,
            y1 - y0,
            facecolor=facecolor,
            edgecolor="black",
            lw=0.8,
            hatch="////",
            zorder=3,
        )
    )


with styles.USGSPlot():
    fig, axes = plt.subplots(ncols=2, figsize=(6.8, 3.7), constrained_layout=True)

    # -- A, the water table falls and the saturated volume decreases
    ax = axes[0]
    draw_cell(ax, XL, HI, "before")
    draw_cell(ax, XR, LO, "after")
    draw_change(ax, XR, LO, HI, DRAINED)
    ax.add_patch(
        FancyArrow(
            XL + W + 0.3,
            4.6,
            1.4,
            0.0,
            width=0.05,
            head_width=0.28,
            head_length=0.28,
            color="black",
            length_includes_head=True,
            zorder=7,
        )
    )
    annotate(
        ax,
        "sorbate of the drained interval is released\n"
        "to the water that remains, and $C$ rises",
        (XR + 0.55, LO - 0.55),
        (XR + W / 2.0, YB - 1.5),
    )
    annotate(
        ax,
        "solute in the water retained in the\n"
        "drained interval is not carried into the\n"
        "next time step, and mass is lost",
        (XR + W - 0.45, (LO + HI) / 2.0 + 0.5),
        (XR + W / 2.0 + 1.5, YT + 2.6),
    )
    styles.heading(ax=ax, letter="A")

    # -- B, the water table rises and the saturated volume increases
    ax = axes[1]
    draw_cell(ax, XL, LO, "before")
    draw_cell(ax, XR, HI, "after")
    draw_change(ax, XR, LO, HI, WATER)
    ax.add_patch(
        FancyArrow(
            XL + W + 0.3,
            4.6,
            1.4,
            0.0,
            width=0.05,
            head_width=0.28,
            head_length=0.28,
            color="black",
            length_includes_head=True,
            zorder=7,
        )
    )
    annotate(
        ax,
        "the interval is saturated again as though\n"
        "it had held water at the concentration of\n"
        "the previous time step, and mass is created",
        (XR + W - 0.45, (LO + HI) / 2.0 + 0.5),
        (XR + W / 2.0 + 1.5, YT + 2.6),
    )
    styles.heading(ax=ax, letter="B")

    for ax in axes:
        ax.set_xlim(0.2, 8.6)
        ax.set_ylim(-1.9, 11.2)
        ax.set_axis_off()

    # -- explanation
    keys = [
        Rectangle((0, 0), 1, 1, facecolor=WATER, edgecolor="black", lw=0.8),
        Rectangle((0, 0), 1, 1, facecolor=DRAINED, edgecolor="black", lw=0.8),
        Rectangle(
            (0, 0), 1, 1, facecolor="white", edgecolor="black", lw=0.8, hatch="////"
        ),
        plt.Line2D([], [], marker="o", ms=3.0, color=SOLID, ls="none"),
    ]
    labels = [
        "saturated part of the cell",
        "drained part of the cell",
        "interval whose saturation changes",
        "sorbate held on the solid aquifer material",
    ]
    leg = fig.legend(
        keys,
        labels,
        loc="lower center",
        ncols=2,
        fontsize=6.5,
        frameon=False,
        handlelength=1.6,
        title="EXPLANATION",
    )
    leg.get_title().set_fontsize(7.5)
    leg.get_title().set_fontweight("bold")

figpth = Path(__file__).resolve().parent.parent / "Figures"
fig.savefig(figpth / "MSTStrandedMassConcept.pdf", dpi=300)
print(f"Saved {figpth / 'MSTStrandedMassConcept.pdf'}")
