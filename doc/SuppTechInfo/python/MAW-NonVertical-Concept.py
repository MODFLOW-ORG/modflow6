"""Generate the conceptual figure that illustrates how the in-cell screen length
(and the saturated conductance) is calculated for non-vertical multi-aquifer
well connections in the Multi-Aquifer Well Package Non-Vertical Connections
chapter of the MODFLOW 6 Supplemental Technical Information document.

The figure has four panels that contrast the radial (THIEM/SKIN/CUMULATIVE) and
cylindrical (MEAN) conductance equations for a slanted and a horizontal
connection:

  A) THIEM/SKIN, slanted connection
  B) MEAN, slanted connection
  C) THIEM/SKIN, horizontal connection
  D) MEAN, horizontal connection

The figure is written to ../Figures/MAWNonVerticalConcept.pdf.
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import styles

figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")

# schematic geometry (exaggerated for clarity)
W, H = 10.0, 10.0  # cell width and height
xc, zc = 5.0, 5.0  # cell center
rw = 0.7  # well radius (exaggerated)
rs = 1.6  # skin (filter pack) radius (exaggerated)


def rect(center, length, width, omega_deg):
    """corner coordinates of a rectangle of the given length (along the screen
    axis) and width (perpendicular), with the axis tilted omega degrees from
    vertical."""
    o = np.radians(omega_deg)
    a = np.array([np.sin(o), np.cos(o)])  # axis direction
    p = np.array([np.cos(o), -np.sin(o)])  # perpendicular direction
    c = np.array(center)
    return np.array(
        [
            c + 0.5 * length * a + 0.5 * width * p,
            c + 0.5 * length * a - 0.5 * width * p,
            c - 0.5 * length * a - 0.5 * width * p,
            c - 0.5 * length * a + 0.5 * width * p,
        ]
    )


def dim_line(ax, p0, p1, text, off, ha="center", va="center", rotation=0):
    """draw a double-headed dimension line between p0 and p1 with a label."""
    ax.annotate(
        "",
        xy=p1,
        xytext=p0,
        arrowprops=dict(arrowstyle="<->", color="black", lw=0.8, shrinkA=0, shrinkB=0),
    )
    mid = 0.5 * (np.array(p0) + np.array(p1)) + np.array(off)
    ax.text(mid[0], mid[1], text, ha=ha, va=va, fontsize=8, rotation=rotation)


def radius_dim(ax, base, direction, length, text, color, label_pad=0.28):
    """draw a single-headed radius arrow from base along direction and label it
    just beyond the arrow tip so the label is clear of the borehole fill."""
    base = np.array(base)
    direction = np.array(direction)
    tip = base + direction * length
    ax.annotate(
        "",
        xy=tip,
        xytext=base,
        arrowprops=dict(arrowstyle="->", color=color, lw=0.8, shrinkA=0, shrinkB=0),
    )
    lab = tip + direction * label_pad
    ax.text(lab[0], lab[1], text, fontsize=8, ha="center", va="center", color=color)


def underbrace(ax, x_left, x_right, y, label, color="0.3"):
    """draw an underbrace spanning [x_left, x_right] at height y with a label."""
    xm = 0.5 * (x_left + x_right)
    ax.plot([x_left, x_right], [y, y], color=color, lw=0.8)
    ax.plot([x_left, x_left], [y, y + 0.12], color=color, lw=0.8)
    ax.plot([x_right, x_right], [y, y + 0.12], color=color, lw=0.8)
    ax.plot([xm, xm], [y, y - 0.12], color=color, lw=0.8)
    ax.text(xm, y - 0.45, label, ha="center", va="top", fontsize=7, color=color)


def draw_radial_eq(ax, cx, cy):
    """draw the THIEM/SKIN/CUMULATIVE conductance equation as a fraction with
    the Thiem and skin terms in the denominator identified by underbraces."""
    # C = and fraction bar
    ax.text(cx - 4.6, cy + 0.15, r"$C = $", ha="right", va="center", fontsize=9)
    ax.plot([cx - 4.3, cx + 4.2], [cy + 0.15, cy + 0.15], color="black", lw=0.8)
    # numerator
    ax.text(cx, cy + 0.55, r"$2 \pi K_h\, L_w$", ha="center", va="center", fontsize=9)
    # denominator: Thiem term + skin term
    xt1 = cx - 2.4
    xt2 = cx + 1.9
    ax.text(xt1, cy - 0.65, r"$\ln(r_0/r_w)$", ha="center", va="center", fontsize=9)
    ax.text(cx - 0.7, cy - 0.65, r"$+$", ha="center", va="center", fontsize=9)
    ax.text(
        xt2,
        cy - 0.65,
        r"$\left(\dfrac{K_h}{K_s}-1\right)\ln(r_s/r_w)$",
        ha="center",
        va="center",
        fontsize=9,
    )
    # underbraces identifying the two denominator components
    underbrace(ax, xt1 - 1.0, xt1 + 1.0, cy - 1.35, "Thiem")
    underbrace(ax, xt2 - 2.0, xt2 + 2.0, cy - 1.35, "skin")


def draw_mean_eq(ax, cx, cy):
    """draw the MEAN conductance equation as a fraction with the numerator and
    denominator clearly separated from the fraction bar."""
    ax.text(cx - 2.5, cy + 0.15, r"$C = $", ha="right", va="center", fontsize=9)
    ax.plot([cx - 2.2, cx + 2.2], [cy + 0.15, cy + 0.15], color="black", lw=0.8)
    ax.text(
        cx,
        cy + 0.6,
        r"$\pi (r_w + r_s) K_s\, L_w$",
        ha="center",
        va="center",
        fontsize=9,
    )
    ax.text(cx, cy - 0.3, r"$r_s - r_w$", ha="center", va="center", fontsize=9)


def draw_panel(ax, omega_deg, mean, lw_in, letter, title):
    # cell
    ax.add_patch(
        mpatches.Rectangle((0, 0), W, H, fc="0.92", ec="0.4", lw=1.0, zorder=0)
    )

    o = np.radians(omega_deg)
    a = np.array([np.sin(o), np.cos(o)])
    horizontal = omega_deg >= 90.0

    # surrounding filter pack (skin) and the well screen (borehole); the skin
    # is shown on all panels because the skin radius and conductivity appear in
    # both the radial (SKIN/CUMULATIVE) and MEAN conductance equations
    ax.add_patch(
        mpatches.Polygon(
            rect((xc, zc), lw_in, 2 * rs, omega_deg),
            closed=True,
            fc="0.75",
            ec="0.5",
            lw=0.6,
            zorder=1,
        )
    )
    ax.add_patch(
        mpatches.Polygon(
            rect((xc, zc), lw_in, 2 * rw, omega_deg),
            closed=True,
            fc="red",
            ec="darkred",
            lw=0.6,
            zorder=2,
        )
    )

    # screen-axis end points
    top = np.array([xc, zc]) + 0.5 * lw_in * a
    bot = np.array([xc, zc]) - 0.5 * lw_in * a

    # dashed well centerline along the screen axis (the radius arrows start here)
    ax.plot(
        [bot[0], top[0]],
        [bot[1], top[1]],
        color="0.25",
        lw=0.7,
        ls=(0, (5, 3)),
        zorder=3,
    )

    # for a slanted connection, show the vertical band occupied by the borehole
    # (scrn_top - scrn_bot) as a vertical dimension at the left side of the cell;
    # for a horizontal connection this band is just the borehole diameter and is
    # represented instead by the well radius (below)
    if not horizontal:
        # vertical extent of the borehole (top to bottom of the pipe), taken
        # from the highest and lowest corners of the slanted screen
        corners = rect((xc, zc), lw_in, 2 * rw, omega_deg)
        itop = int(np.argmax(corners[:, 1]))
        ibot = int(np.argmin(corners[:, 1]))
        ztop, zbot = corners[itop, 1], corners[ibot, 1]
        xdim = 1.2
        # witness lines connecting the screen top and bottom to the dimension line
        ax.plot(
            [xdim, corners[itop, 0]], [ztop, ztop], color="0.5", lw=0.5, ls=(0, (3, 2))
        )
        ax.plot(
            [xdim, corners[ibot, 0]], [zbot, zbot], color="0.5", lw=0.5, ls=(0, (3, 2))
        )
        dim_line(
            ax,
            (xdim, zbot),
            (xdim, ztop),
            r"$z_t - z_b$",
            off=(-0.35, 0),
            ha="right",
            rotation=90,
        )

    # in-cell screen length L_w along the screen axis
    noff = np.array([np.cos(o), -np.sin(o)]) * (rw + 0.9)
    dim_line(
        ax,
        bot + noff,
        top + noff,
        r"$L_w$",
        off=(0.55, 0.0) if not horizontal else (0.0, -0.6),
        ha="left" if not horizontal else "center",
    )

    # well radius (and skin radius for MEAN), drawn perpendicular to the screen
    # axis from the axis to the edge of the borehole (and filter pack)
    pneg = np.array([-np.cos(o), np.sin(o)])
    if horizontal:
        rbase = np.array([xc - 0.32 * lw_in, zc])
    else:
        rbase = np.array([xc, zc]) - 0.18 * lw_in * a
    # offset the rw and rs arrows along the screen axis so both are visible
    radius_dim(ax, rbase + 1.0 * a, pneg, rs, r"$r_s$", color="0.3")
    radius_dim(ax, rbase - 1.0 * a, pneg, rw, r"$r_w$", color="darkred")

    # tilt angle omega between vertical and the screen axis
    if not horizontal:
        ax.plot([xc, xc], [zc, zc + 0.5 * lw_in], color="0.4", lw=0.6, ls=":")
        arc = mpatches.Arc(
            (xc, zc), 3.0, 3.0, angle=0, theta1=90 - omega_deg, theta2=90, color="0.2"
        )
        ax.add_patch(arc)
        ax.text(xc - 0.55, zc + 1.5, r"$\omega$", fontsize=9, color="0.2")

    # labels: how L_w is obtained and the conductance equation
    if horizontal:
        lw_txt = r"$L_w$ = conn$\_$length (specified)"
    else:
        lw_txt = r"$L_w = (z_t - z_b - 2 r_w \sin\omega)\,/\,\cos\omega$"
    ax.text(0.5 * W, -0.7, lw_txt, ha="center", va="top", fontsize=8)

    if mean:
        draw_mean_eq(ax, cx=0.5 * W, cy=-2.3)
        pp_txt = r"filter pack: $K_s$"
    else:
        # THIEM / SKIN / CUMULATIVE: Thiem term plus skin term in the denominator
        draw_radial_eq(ax, cx=0.5 * W, cy=-2.3)
        pp_txt = r"aquifer: $K_h$, $r_0$;  skin: $K_s$"

    ax.text(0.5 * W, H + 0.4, pp_txt, ha="center", va="bottom", fontsize=8)

    styles.heading(ax=ax, letter=letter, heading=title)
    ax.set_xlim(-0.5, W + 0.5)
    ax.set_ylim(-4.8, H + 1.6)
    ax.set_aspect("equal", "box")
    ax.axis("off")


def make_figure():
    with styles.USGSPlot():
        fig, axes = plt.subplots(2, 2, figsize=(7.5, 8.5))

        # slanted in-cell length from the radius-aware formula
        omega = 40.0
        zt_zb = 6.0
        lw_slant = (zt_zb - 2 * rw * np.sin(np.radians(omega))) / np.cos(
            np.radians(omega)
        )
        lw_horiz = 7.5  # specified connection length

        draw_panel(
            axes[0, 0],
            omega,
            mean=False,
            lw_in=lw_slant,
            letter="A",
            title="THIEM / SKIN, slanted",
        )
        draw_panel(
            axes[0, 1],
            omega,
            mean=True,
            lw_in=lw_slant,
            letter="B",
            title="MEAN, slanted",
        )
        draw_panel(
            axes[1, 0],
            90.0,
            mean=False,
            lw_in=lw_horiz,
            letter="C",
            title="THIEM / SKIN, horizontal",
        )
        draw_panel(
            axes[1, 1],
            90.0,
            mean=True,
            lw_in=lw_horiz,
            letter="D",
            title="MEAN, horizontal",
        )

        fig.tight_layout()

    fpth = os.path.join(figpth, "MAWNonVerticalConcept.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    make_figure()


if __name__ == "__main__":
    main()
