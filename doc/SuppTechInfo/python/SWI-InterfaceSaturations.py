"""Generate the conceptual figure that illustrates the sharp-interface geometry
and the bottom-referenced saturations used by the Seawater Intrusion (SWI)
Package chapter of the MODFLOW 6 Supplemental Technical Information document.

The figure shows a single aquifer cell on the left (freshwater above the
interface, saltwater below) and, on the right, the three saturations as vertical
"measuring bars":

  S^w = S(h^f)      total-water column, measured from the cell bottom Z_B
  S^s = S(zeta)     salt column, measured from the cell bottom Z_B
  S^f = S^w - S^s   freshwater fraction, bounded below by the moving interface

S^w and S^s rise from the datum (Z_B), so both fit the standard MODFLOW 6
storage kernels; S^f floats above the datum (its base is the interface zeta),
which is why the freshwater storage is assembled as the difference of the two
bottom-referenced columns rather than by passing S^f to the kernels directly.

The figure is written to ../Figures/SWIInterfaceSaturations.pdf.
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from flopy.plot import styles

figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")

# colors (kept muted, consistent with the other schematic chapters)
C_FRESH = "#4a90d9"
C_SALT = "#3aa6a0"

# cell geometry (schematic elevations)
ZB = 0.0  # cell bottom (datum)
ZT = 5.0  # cell top
HF = 3.7  # freshwater head (water table)
ZETA = 1.8  # interface elevation


def elev_labels(ax, x):
    """right-justified elevation labels just left of x."""
    for y, txt, col in (
        (ZT, r"$Z_T$", "black"),
        (HF, r"$h^f$", C_FRESH),
        (ZETA, r"$\zeta$", C_SALT),
        (ZB, r"$Z_B$", "black"),
    ):
        ax.text(x - 0.12, y, txt, ha="right", va="center", fontsize=10, color=col)


def bar(ax, x, ybot, ytop, color, floating=False):
    """draw a saturation measuring bar from ybot to ytop at horizontal position x.
    A floating bar (base above the datum) is drawn with a dashed outline and a
    dotted stem down to the datum to emphasize that its base is not Z_B."""
    w = 0.5
    if floating:
        ax.add_patch(
            mpatches.Rectangle(
                (x - w / 2, ybot),
                w,
                ytop - ybot,
                fc="none",
                ec="#c0392b",
                lw=1.3,
                ls=(0, (4, 2)),
                zorder=3,
            )
        )
        ax.plot([x, x], [ZB, ybot], color="#c0392b", lw=0.8, ls=(0, (1, 2)), zorder=2)
    else:
        ax.add_patch(
            mpatches.Rectangle(
                (x - w / 2, ybot),
                w,
                ytop - ybot,
                fc=color,
                ec="0.3",
                lw=0.8,
                alpha=0.85,
                zorder=3,
            )
        )


def make_figure():
    with styles.USGSPlot():
        fig, ax = plt.subplots(1, 1, figsize=(7.4, 4.2))

        # -- datum line across the whole figure
        ax.plot([-0.2, 8.8], [ZB, ZB], color="0.5", lw=0.8, ls=(0, (5, 3)), zorder=0)
        ax.text(
            8.8, ZB - 0.3, "datum $Z_B$", ha="right", va="top", fontsize=8, color="0.5"
        )

        # -- the aquifer cell (left)
        xL, xR = 0.6, 3.0
        # saltwater zone (Z_B -> zeta)
        ax.add_patch(
            mpatches.Rectangle(
                (xL, ZB), xR - xL, ZETA - ZB, fc=C_SALT, ec="none", alpha=0.35, zorder=1
            )
        )
        # freshwater zone (zeta -> h^f)
        ax.add_patch(
            mpatches.Rectangle(
                (xL, ZETA),
                xR - xL,
                HF - ZETA,
                fc=C_FRESH,
                ec="none",
                alpha=0.30,
                zorder=1,
            )
        )
        # cell outline
        ax.add_patch(
            mpatches.Rectangle(
                (xL, ZB), xR - xL, ZT - ZB, fc="none", ec="0.3", lw=1.2, zorder=4
            )
        )
        # water table and interface surfaces
        ax.plot([xL, xR], [HF, HF], color=C_FRESH, lw=2.0, zorder=5)
        ax.plot([xL, xR], [ZETA, ZETA], color=C_SALT, lw=2.0, zorder=5)
        # zone labels
        ax.text(
            0.5 * (xL + xR),
            0.5 * (ZT + HF),
            "unsaturated",
            ha="center",
            va="center",
            fontsize=8,
            color="0.4",
        )
        ax.text(
            0.5 * (xL + xR),
            0.5 * (HF + ZETA),
            "freshwater",
            ha="center",
            va="center",
            fontsize=9,
            color="#1f5fa8",
        )
        ax.text(
            0.5 * (xL + xR),
            0.5 * (ZETA + ZB),
            "saltwater",
            ha="center",
            va="center",
            fontsize=9,
            color="#1c6b66",
        )
        elev_labels(ax, xL)

        # -- the three saturation measuring bars (right)
        xw, xs, xf = 5.0, 6.4, 7.8
        bar(ax, xw, ZB, HF, C_FRESH)
        bar(ax, xs, ZB, ZETA, C_SALT)
        bar(ax, xf, ZETA, HF, None, floating=True)

        # "from bottom" / "floating base" annotations
        for x, txt in ((xw, "from $Z_B$"), (xs, "from $Z_B$")):
            ax.text(x, -0.55, txt, ha="center", va="top", fontsize=7, color="0.35")
        ax.text(
            xf,
            -0.55,
            "floating base",
            ha="center",
            va="top",
            fontsize=7,
            color="#c0392b",
        )

        # bar equations
        ax.text(
            xw,
            HF + 0.28,
            r"$S^w=\mathcal{S}(h^f)$",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#1f5fa8",
        )
        ax.text(
            xs,
            ZETA + 0.28,
            r"$S^s=\mathcal{S}(\zeta)$",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#1c6b66",
        )
        ax.text(
            xf,
            HF + 0.28,
            r"$S^f=S^w-S^s$",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#c0392b",
        )

        ax.set_xlim(-0.3, 9.0)
        ax.set_ylim(-1.4, 6.1)
        ax.set_aspect("equal", "box")
        ax.axis("off")
        fig.tight_layout()

    fpth = os.path.join(figpth, "SWIInterfaceSaturations.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    make_figure()


if __name__ == "__main__":
    main()
