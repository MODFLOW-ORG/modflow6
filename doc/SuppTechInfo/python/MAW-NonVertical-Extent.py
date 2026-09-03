"""Generate the figure that illustrates the horizontal reach of a non-vertical
multi-aquifer well connection in the Multi-Aquifer Well Package Non-Vertical
Connections chapter of the MODFLOW 6 Supplemental Technical Information
document.

Both panels are drawn to the same scale for cells 20 long and 10 thick, with a
screen that spans the cell thickness:

  A) a moderate tilt, where the screen stays within the connected cell
  B) a strong tilt, where the screen reaches across several cells and a
     connection to each of them is required

The figure is written to ../Figures/MAWNonVerticalExtent.pdf.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import styles

figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")

DX, DZ = 20.0, 10.0  # cell length and thickness
RW = 0.25  # well radius (exaggerated)
NCELL = 6  # cells drawn in both panels


def screen_length(dz, rw, omega_deg):
    """in-cell screen length from the vertical screen extent and tilt angle."""
    o = np.radians(omega_deg)
    return (dz - 2.0 * rw * np.sin(o)) / np.cos(o)


def dim_line(ax, p0, p1, text, off, **kwargs):
    """draw a double-headed dimension line between p0 and p1 with a label."""
    ax.annotate(
        "",
        xy=p1,
        xytext=p0,
        arrowprops=dict(arrowstyle="<->", color="black", lw=0.8, shrinkA=0, shrinkB=0),
    )
    mid = 0.5 * (np.array(p0) + np.array(p1)) + np.array(off)
    ax.text(mid[0], mid[1], text, **kwargs)


def draw_panel(ax, omega_deg, letter, title):
    lw_in = screen_length(DZ, RW, omega_deg)
    o = np.radians(omega_deg)
    run = lw_in * np.sin(o)  # horizontal reach of the screen
    reached = int(np.ceil(run / DX))  # cells the screen passes through

    # the screen enters at the left edge of the first cell; a cell the screen
    # reaches but is not connected to is shaded
    for i in range(NCELL):
        ax.add_patch(
            plt.Rectangle(
                (i * DX, 0.0),
                DX,
                DZ,
                facecolor="white" if i == 0 else ("0.90" if i < reached else "0.97"),
                edgecolor="black",
                lw=0.8,
                zorder=1,
            )
        )
    ax.text(0.5 * DX, DZ + 1.0, "connected\ncell", ha="center", va="bottom", size=6)
    if reached > 1:
        ax.text(
            0.5 * (1 + reached) * DX,
            DZ + 1.0,
            "reached by the screen, but not connected",
            ha="center",
            va="bottom",
            size=6,
        )

    # the screen, drawn as a band of width equal to the well diameter
    x0, z0 = 0.0, DZ
    x1, z1 = run, DZ - lw_in * np.cos(o)
    perp = np.array([np.cos(o), -np.sin(o)]) * RW
    band = np.array(
        [[x0, z0] + perp, [x1, z1] + perp, [x1, z1] - perp, [x0, z0] - perp]
    )
    ax.add_patch(plt.Polygon(band, closed=True, facecolor="red", ec="red", zorder=3))

    # the screen elevations, labeled clear of the cells
    xr = NCELL * DX
    for z, label in ((z0, "$TOPSCRN$"), (z1, "$BOTSCRN$")):
        ax.plot([0.0, xr], [z, z], "k--", lw=0.6, zorder=2)
        ax.text(xr + 0.06 * DX, z, label, va="center", size=7, style="italic")
    ax.annotate(
        "$L_w$",
        xy=(0.45 * (x0 + x1), 0.5 * (z0 + z1)),
        xytext=(0, 9),
        textcoords="offset points",
        ha="center",
        size=8,
    )

    # the tilt angle, measured from vertical at the top of the screen
    ax.plot([x0, x0], [z0, z0 - 0.5 * DZ], "k:", lw=0.6, zorder=2)
    ax.annotate(r"$\omega$", xy=(x0 + 0.04 * DX, z0 - 0.34 * DZ), size=8)

    # the horizontal reach
    dim_line(
        ax,
        (x0, -0.30 * DZ),
        (x1, -0.30 * DZ),
        r"$L_w \sin \omega$",
        (0.0, -0.14 * DZ),
        ha="center",
        va="top",
        size=8,
    )

    ax.set_xlim(-0.05 * DX, (NCELL + 0.6) * DX)
    ax.set_ylim(-0.75 * DZ, DZ + 0.55 * DZ)
    ax.set_aspect("equal")
    ax.axis("off")
    styles.heading(ax=ax, letter=letter, heading=title)


def make_figure():
    with styles.USGSMap():
        fig, axes = plt.subplots(nrows=2, figsize=(7.0, 3.6))
        draw_panel(axes[0], 60.0, "A", r"$\omega = 60^{\circ}$, screen within the cell")
        draw_panel(axes[1], 85.0, "B", r"$\omega = 85^{\circ}$, screen beyond the cell")
        fig.tight_layout()

    fpth = os.path.join(figpth, "MAWNonVerticalExtent.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    make_figure()


if __name__ == "__main__":
    main()
