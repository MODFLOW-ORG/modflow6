"""Generate the conceptual figure that illustrates the two-fluid vertical flow
restriction used by the Seawater Intrusion (SWI) Package chapter of the
MODFLOW 6 Supplemental Technical Information document.

Two cells are stacked and share a horizontal face. In the sharp-interface
conceptualization each cell contains freshwater over saltwater, so the shown
configuration is buoyantly unstable: the upper cell's saltwater lies directly
above the lower cell's freshwater at the face. The vertical connection carries
freshwater flow between the two freshwater zones and saltwater flow between the
two saltwater zones. The restriction follows the same rule as the vertical
leakage of Huyakorn, Wu, and Park (1996): same-fluid vertical flow is allowed,
but the lighter fluid may not move down through the heavier fluid and the
heavier fluid may not move up through the lighter fluid.

  A) Buoyant flow (allowed): freshwater rises through the face into the upper
     freshwater and saltwater sinks through the face into the lower saltwater;
     the two flows bypass one another.
  B) Anti-buoyant flow (restricted): freshwater cannot descend through the
     saltwater and saltwater cannot rise through the freshwater.

Written to ../Figures/SWIVerticalRestriction.pdf.
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from flopy.plot import styles

figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")

C_FRESH = "#4a90d9"
C_SALT = "#3aa6a0"
C_FRESH_DK = "#1f5fa8"
C_SALT_DK = "#1c6b66"
C_BLOCK = "#c0392b"

XL, XR = 0.0, 3.0  # cell horizontal extent
YB, YF, YT = 0.0, 3.0, 6.0  # lower-cell bottom, shared face, upper-cell top
XC = 0.5 * (XL + XR)
ZU = 4.5  # interface in the upper cell (fresh above, salt below to the face)
ZL = 1.5  # interface in the lower cell (fresh above to the face, salt below)


def rect(ax, y0, y1, color):
    ax.add_patch(
        mpatches.Rectangle(
            (XL, y0), XR - XL, y1 - y0, fc=color, ec="none", alpha=0.34, zorder=1
        )
    )


def draw_cells(ax):
    # each cell: freshwater over saltwater (both fluids in both cells)
    rect(ax, ZU, YT, C_FRESH)  # upper cell freshwater
    rect(ax, YF, ZU, C_SALT)  # upper cell saltwater (sits at the face)
    rect(ax, ZL, YF, C_FRESH)  # lower cell freshwater (sits at the face)
    rect(ax, YB, ZL, C_SALT)  # lower cell saltwater
    # cell outlines
    for y0, y1 in ((YF, YT), (YB, YF)):
        ax.add_patch(
            mpatches.Rectangle(
                (XL, y0), XR - XL, y1 - y0, fc="none", ec="0.3", lw=1.2, zorder=4
            )
        )
    # interface lines and the shared face
    for z in (ZU, ZL):
        ax.plot([XL, XR], [z, z], color=C_SALT_DK, lw=1.6, zorder=5)
    ax.plot([XL, XR], [YF, YF], color="0.15", lw=2.4, zorder=5)
    ax.text(
        XR + 0.12,
        YF,
        "shared\nface",
        ha="left",
        va="center",
        fontsize=7.5,
        color="0.15",
    )
    # fluid labels to the left of each band (clear of the arrows)
    for y, txt, col in (
        (0.5 * (ZU + YT), "fresh", C_FRESH_DK),
        (0.5 * (YF + ZU), "salt", C_SALT_DK),
        (0.5 * (ZL + YF), "fresh", C_FRESH_DK),
        (0.5 * (YB + ZL), "salt", C_SALT_DK),
    ):
        ax.text(XL - 0.15, y, txt, ha="right", va="center", fontsize=7.5, color=col)


def arrow(ax, x, y0, y1, color, restricted=False):
    """flow arrow from y0 to y1; if restricted, draw it faded/dashed and add a
    red cross at the face where the crossing is blocked."""
    if restricted:
        ax.annotate(
            "",
            xy=(x, y1),
            xytext=(x, y0),
            arrowprops=dict(
                arrowstyle="->", color=color, lw=1.4, ls=(0, (3, 2)), alpha=0.55
            ),
        )
        ax.plot(x, YF, marker="x", ms=11, mew=2.6, color=C_BLOCK, zorder=8)
    else:
        ax.annotate(
            "",
            xy=(x, y1),
            xytext=(x, y0),
            arrowprops=dict(arrowstyle="->", color=color, lw=2.2),
        )


def base(ax, letter, title, note):
    styles.heading(ax=ax, letter=letter, heading=title)
    ax.text(XC, YB - 0.75, note, ha="center", va="top", fontsize=7.5, color="0.25")
    ax.set_xlim(XL - 1.3, XR + 1.6)
    ax.set_ylim(YB - 1.7, YT + 0.5)
    ax.set_aspect("equal", "box")
    ax.axis("off")


def panel_A(ax):
    draw_cells(ax)
    # freshwater rises (lower fresh -> up through salt -> upper fresh)
    arrow(ax, XC - 0.7, 2.2, 5.2, C_FRESH_DK)
    ax.text(
        XC - 0.95,
        5.35,
        "freshwater\nrises",
        ha="center",
        va="bottom",
        fontsize=8,
        color=C_FRESH_DK,
    )
    # saltwater sinks (upper salt -> down through fresh -> lower salt)
    arrow(ax, XC + 0.7, 3.8, 0.8, C_SALT_DK)
    ax.text(
        XC + 0.95,
        0.65,
        "saltwater\nsinks",
        ha="center",
        va="top",
        fontsize=8,
        color=C_SALT_DK,
    )
    base(
        ax,
        "A",
        "buoyant flow: allowed",
        "each fluid crosses in its buoyant direction; the two flows bypass",
    )


def panel_B(ax):
    draw_cells(ax)
    # freshwater cannot descend through the saltwater
    arrow(ax, XC - 0.7, 5.2, 2.2, C_FRESH_DK, restricted=True)
    ax.text(
        XC - 0.95,
        5.35,
        "freshwater\ndown",
        ha="center",
        va="bottom",
        fontsize=8,
        color=C_FRESH_DK,
    )
    # saltwater cannot rise through the freshwater
    arrow(ax, XC + 0.7, 0.8, 3.8, C_SALT_DK, restricted=True)
    ax.text(
        XC + 0.95,
        0.65,
        "saltwater\nup",
        ha="center",
        va="top",
        fontsize=8,
        color=C_SALT_DK,
    )
    base(
        ax,
        "B",
        "anti-buoyant flow: restricted",
        "freshwater cannot sink through salt; salt cannot rise through fresh",
    )


def make_figure():
    with styles.USGSPlot():
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.8))
        panel_A(axes[0])
        panel_B(axes[1])
        fig.tight_layout()

    fpth = os.path.join(figpth, "SWIVerticalRestriction.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    make_figure()


if __name__ == "__main__":
    main()
