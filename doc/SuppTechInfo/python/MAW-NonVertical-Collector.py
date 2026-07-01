"""Generate the radial collector (Ranney) well figure for the Multi-Aquifer
Well Package Non-Vertical Connections chapter of the MODFLOW 6 Supplemental
Technical Information document.

The figure shows a single multi-aquifer well with a vertical caisson and four
horizontal laterals (non-vertical connections) and the simulated head field.
The MAW and CHD boundary cells are drawn with ``flopy.plot.PlotMapView.plot_bc``
and the figure is styled with ``flopy.plot.styles``.

The figure is written to ../Figures/MAWNonVerticalCollector.pdf.
"""

import os
import shutil

import flopy
import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import styles

# -- paths
ws = os.path.join(os.path.dirname(__file__), "temp", "maw-collector")
figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")
mf6exe = os.path.join(os.path.dirname(__file__), "..", "..", "..", "bin", "mf6")

# -- grid and aquifer properties
name = "collector"
nlay, nrow, ncol = 1, 41, 41
delr = delc = 25.0
top = 50.0
botm = 0.0
hk = 25.0
strt = 50.0

# -- multi-aquifer well properties
radius = 0.5
sradius = 1.0
hks = 25.0
mawrate = -25000.0
arm = 12  # number of cells each lateral extends from the central caisson
center = (nrow // 2, ncol // 2)


def lateral_cells():
    """cells traversed by the four horizontal laterals (excluding the caisson)."""
    cells = []
    ic, jc = center
    for d in range(1, arm + 1):
        cells.append((0, ic - d, jc))  # north
        cells.append((0, ic + d, jc))  # south
        cells.append((0, ic, jc - d))  # west
        cells.append((0, ic, jc + d))  # east
    return cells


def build_model():
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6exe, sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-9,
        inner_dvclose=1e-9,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk, save_flows=True)

    # constant heads around the perimeter
    chdspd = [
        [(0, i, j), strt]
        for i in range(nrow)
        for j in range(ncol)
        if i in (0, nrow - 1) or j in (0, ncol - 1)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd, pname="CHD")

    # multi-aquifer well: a vertical caisson connection plus horizontal laterals
    laterals = lateral_cells()
    ncon = 1 + len(laterals)
    packagedata = [[0, radius, botm, strt, "MEAN", ncon]]

    # a thin vertical band represents the horizontal screen within a cell
    thin = 2.0 * radius
    sc_top = 0.5 * (top + botm) + 0.5 * thin
    sc_bot = 0.5 * (top + botm) - 0.5 * thin

    connectiondata = []
    angledata = []
    # caisson connection (vertical, spans the full cell thickness)
    connectiondata.append([0, 0, (0, center[0], center[1]), top, botm, hks, sradius])
    # lateral connections (horizontal: 90 degrees, in-cell length = cell width)
    for icon, cellid in enumerate(laterals, start=1):
        connectiondata.append([0, icon, cellid, sc_top, sc_bot, hks, sradius])
        angledata.append([0, icon, 90.0, delr])

    flopy.mf6.ModflowGwfmaw(
        gwf,
        non_vertical_wells=True,
        save_flows=True,
        print_head=True,
        packagedata=packagedata,
        connectiondata=connectiondata,
        angledata=angledata,
        perioddata={0: [[0, "rate", mawrate]]},
        pname="MAW",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim, gwf


def make_figure(gwf):
    head = gwf.output.head().get_data()

    with styles.USGSMap():
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 6.0), tight_layout=True)

        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
        hmin = float(np.nanmin(head))
        ca = pmv.plot_array(head, cmap="viridis_r", alpha=0.85)
        levels = np.linspace(hmin, strt, 9)
        cs = pmv.contour_array(head, levels=levels, colors="white", linewidths=0.75)
        ax.clabel(cs, fmt="%.1f", fontsize=8)
        pmv.plot_grid(lw=0.2, color="0.7")
        pmv.plot_bc("CHD", color="0.4")
        pmv.plot_bc("MAW", color="red")

        # mark the central caisson
        xc = gwf.modelgrid.xcellcenters[center]
        yc = gwf.modelgrid.ycellcenters[center]
        ax.plot(xc, yc, marker="o", mfc="white", mec="black", ms=7, mew=1.0, zorder=5)

        ax.set_aspect("equal", "box")
        styles.xlabel(ax=ax, label="x, in meters")
        styles.ylabel(ax=ax, label="y, in meters")

        cbar = fig.colorbar(ca, ax=ax, shrink=0.7)
        cbar.ax.set_title("Head, m", pad=10, fontdict={"size": 8})

        # legend proxies for the boundary conditions
        handles = [
            plt.Line2D([], [], color="red", lw=6, label="MAW lateral cells"),
            plt.Line2D([], [], color="0.4", lw=6, label="Constant-head cells"),
            plt.Line2D(
                [],
                [],
                color="black",
                marker="o",
                mfc="white",
                ls="",
                ms=7,
                label="MAW caisson",
            ),
        ]
        styles.graph_legend(ax=ax, handles=handles, loc="upper right")

    fpth = os.path.join(figpth, "MAWNonVerticalCollector.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    if os.path.isdir(ws):
        shutil.rmtree(ws)
    os.makedirs(ws, exist_ok=True)
    sim, gwf = build_model()
    sim.write_simulation()
    success, buff = sim.run_simulation(silent=True)
    if not success:
        raise RuntimeError("MODFLOW 6 did not terminate normally")
    make_figure(gwf)


if __name__ == "__main__":
    main()
