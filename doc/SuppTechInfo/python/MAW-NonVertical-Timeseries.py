"""Generate the radial collector (Ranney) well head time-series figure for the
Multi-Aquifer Well Package Non-Vertical Connections chapter of the MODFLOW 6
Supplemental Technical Information document.

A transient simulation (pumping followed by recovery) of the radial collector
well is run two ways: (1) with the horizontal laterals represented as
non-vertical connections (the length correction is applied), and (2) with the
same lateral cells represented as thin vertical screens (no length correction).

The figure (styled with ``flopy.plot.styles.USGSPlot``) shows the simulated
well head with and without the length correction. Representing the laterals as
thin vertical screens (no length correction) greatly reduces the connection
conductance and produces substantially more drawdown in the well for the same
pumping rate.

The figure is written to ../Figures/MAWNonVerticalTimeseries.pdf.
"""

import os
import shutil

import flopy
import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import styles

# -- paths
basews = os.path.join(os.path.dirname(__file__), "temp")
figpth = os.path.join(os.path.dirname(__file__), "..", "Figures")
mf6exe = os.path.join(os.path.dirname(__file__), "..", "..", "..", "bin", "mf6")

# -- grid and aquifer properties
nlay, nrow, ncol = 1, 41, 41
delr = delc = 25.0
top = 50.0
botm = 0.0
hk = 25.0
ss = 1.0e-4
strt = 50.0

# -- multi-aquifer well properties
radius = 0.5
sradius = 1.0
hks = 25.0
pumprate = -40000.0
arm = 12
center = (nrow // 2, ncol // 2)
caisson_cell = (0, center[0], center[1])

# -- transient stress periods: pumping then recovery
perioddata_tdis = [(5.0, 25, 1.2), (5.0, 25, 1.2)]


def lateral_cells():
    cells = []
    ic, jc = center
    for d in range(1, arm + 1):
        cells.append((0, ic - d, jc))
        cells.append((0, ic + d, jc))
        cells.append((0, ic, jc - d))
        cells.append((0, ic, jc + d))
    return cells


def build_model(corrected):
    name = "corrected" if corrected else "uncorrected"
    ws = os.path.join(basews, f"ts-{name}")
    if os.path.isdir(ws):
        shutil.rmtree(ws)
    os.makedirs(ws, exist_ok=True)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6exe, sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(perioddata_tdis),
        perioddata=perioddata_tdis,
    )
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
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=ss, sy=0.0, transient={0: True})

    chdspd = [
        [(0, i, j), strt]
        for i in range(nrow)
        for j in range(ncol)
        if i in (0, nrow - 1) or j in (0, ncol - 1)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd, pname="CHD")

    laterals = lateral_cells()
    ncon = 1 + len(laterals)
    packagedata = [[0, radius, botm, strt, "MEAN", ncon]]

    thin = 2.0 * radius
    sc_top = 0.5 * (top + botm) + 0.5 * thin
    sc_bot = 0.5 * (top + botm) - 0.5 * thin

    connectiondata = [[0, 0, caisson_cell, top, botm, hks, sradius]]
    angledata = []
    for icon, cellid in enumerate(laterals, start=1):
        connectiondata.append([0, icon, cellid, sc_top, sc_bot, hks, sradius])
        angledata.append([0, icon, 90.0, delr])

    maw_kwargs = dict(
        no_well_storage=True,
        save_flows=True,
        print_head=True,
        observations={f"{name}.maw.obs.csv": [("hwell", "head", 1)]},
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [[0, "rate", pumprate]], 1: [[0, "rate", 0.0]]},
        pname="MAW",
    )
    if corrected:
        # represent the laterals as horizontal non-vertical connections
        maw_kwargs["non_vertical_wells"] = True
        maw_kwargs["angledata"] = angledata
    # when not corrected the same lateral cells are connected with thin
    # vertical screens (no ANGLEDATA), so the in-cell screen length is the thin
    # vertical band rather than the lateral length

    flopy.mf6.ModflowGwfmaw(gwf, **maw_kwargs)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim, gwf, ws, name


def run(corrected):
    sim, _, ws, name = build_model(corrected)
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    if not success:
        raise RuntimeError(f"MODFLOW 6 ({name}) did not terminate normally")
    # well head time series from the MAW observation file
    obs = np.genfromtxt(
        os.path.join(ws, f"{name}.maw.obs.csv"), delimiter=",", names=True
    )
    return obs["time"], obs["HWELL"]


def make_figure():
    totim, hwell = run(corrected=True)
    _, hwell_u = run(corrected=False)
    tpump = perioddata_tdis[0][0]

    with styles.USGSPlot():
        fig, ax = plt.subplots(1, 1, figsize=(5.0, 4.0), tight_layout=True)

        # shade the pumping (red) and recovery (blue) portions of the simulation
        ax.axvspan(0.0, tpump, color="red", alpha=0.5, lw=0.0, zorder=0)
        ax.axvspan(tpump, totim[-1], color="blue", alpha=0.5, lw=0.0, zorder=0)

        ax.plot(totim, hwell, color="black", lw=1.5, label="With length correction")
        ax.plot(
            totim,
            hwell_u,
            color="black",
            lw=1.5,
            ls="--",
            label="Without length correction",
        )
        ax.set_xlim(0, totim[-1])
        styles.xlabel(ax=ax, label="Time, in days")
        styles.ylabel(ax=ax, label="MAW well head, in meters")
        styles.add_text(
            ax=ax, text="pumping", x=0.02, y=0.96, ha="left", va="top", bold=False
        )
        styles.add_text(
            ax=ax, text="recovery", x=0.98, y=0.04, ha="right", va="bottom", bold=False
        )
        styles.graph_legend(ax=ax, loc="upper center", bbox_to_anchor=(0.5, -0.18))

    fpth = os.path.join(figpth, "MAWNonVerticalTimeseries.pdf")
    fig.savefig(fpth)
    print(f"saved {fpth}")
    plt.close(fig)


def main():
    os.makedirs(basews, exist_ok=True)
    make_figure()


if __name__ == "__main__":
    main()
