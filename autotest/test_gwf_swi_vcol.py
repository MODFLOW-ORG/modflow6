"""
Vertical-column test of SWI vertical flow, single-fluid (freshwater-only) case.

A 1-D vertical column of confined cells is driven by a freshwater well at the
top: injection for the first stress period pushes the flat interface down, and
extraction for the second period pulls it back up. Because the interface storage
is drainable (SY) and the column is closed, the interface elevation moves at a
rate set purely by the mass balance, dz/dt = Q / (SY * A). With Q = 0.02,
SY = 0.2 and a unit cell area the interface marches 0.1 m per step -- down
through two layer faces and then back up to its start.

In the single-fluid case the vertical flow is NOT buoyancy-restricted: the
interface is a diagnostic derived from a continuous freshwater potential, so the
vertical conductance is the standard NPF conductance. Restricting it would
decouple the deep cells below the interface and leave the head field
underdetermined (the column would not converge on flow reversal, and would need
an artificial bottom boundary). With unrestricted vertical flow the closed
column is well posed and needs no bottom boundary. A high K keeps the (small)
vertical gradient that feeds the
interface-cell storage negligible, so the interface stays essentially flat and
the per-cell interface elevation zeta = -alphaf*hf tracks it.

Run with --plot to write a figure of the interface elevation versus time
(interface marching down then up through the layer faces).
"""

import pathlib as pl

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["swi-vcol"]

nlay, nrow, ncol = 5, 1, 1
delr = delc = 1.0
top = 0.0
botm = [-1.0, -2.0, -3.0, -4.0, -5.0]
sy = 0.2
alphaf = 1000.0 / (1025.0 - 1000.0)  # 40.0

nstp = 20
q = 0.02  # injection/extraction rate -> dz = q / (sy * A) = 0.1 m per step
zeta0 = -0.5  # initial (flat) interface elevation
dz_expected = q / (sy * (delr * delc))  # 0.1 m per step


def build_models(idx, test):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    # two periods (inject then extract), nstp steps each, dt = 1
    flopy.mf6.ModflowTdis(
        sim, nper=2, perioddata=[(nstp, nstp, 1.0), (nstp, nstp, 1.0)]
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=100,
        inner_maximum=100,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-9,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
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
    flopy.mf6.ModflowGwfic(gwf, strt=-zeta0 / alphaf)  # zeta = -0.5 initially
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=100.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=sy, transient=True)
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    # freshwater well in the top cell: inject then extract
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={0: [[(0, 0, 0), q]], 1: [[(0, 0, 0), -q]]},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budgetcsv_filerecord=name + ".bud.csv",
        saverecord=[("HEAD", "ALL")],
    )
    return sim, None


def _zeta(test):
    """per-cell interface elevation over time: zeta = -alphaf*hf."""
    gwf = test.sims[0].gwf[0]
    hobj = gwf.output.head()
    times = hobj.times
    z = np.array([-alphaf * hobj.get_data(totim=t).flatten() for t in times])
    return np.array(times), z  # z shape (ntime, nlay)


def plot_output(idx, test):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    ws = pl.Path(test.workspace)
    # per-cell (clamped) interface: each cell splits into freshwater (blue,
    # above zeta) and saltwater (red, below zeta) using the .zta output
    zobj = flopy.utils.HeadFile(ws / "mymodel.zta", text="zeta")
    times = zobj.times
    cell_top = np.array([top] + list(botm[:-1]))
    cell_bot = np.array(botm)

    snap = np.linspace(0, len(times) - 1, 9).round().astype(int)
    fig = plt.figure(figsize=(12, 3.2))
    for j, ti in enumerate(snap):
        ax = fig.add_subplot(1, len(snap), j + 1)
        zeta = zobj.get_data(totim=times[ti]).flatten()
        for k in range(nlay):
            z = min(max(zeta[k], cell_bot[k]), cell_top[k])
            ax.add_patch(
                Rectangle(
                    (0, cell_bot[k]),
                    1.0,
                    z - cell_bot[k],
                    facecolor="red",
                    edgecolor="k",
                    lw=0.5,
                )
            )
            ax.add_patch(
                Rectangle(
                    (0, z),
                    1.0,
                    cell_top[k] - z,
                    facecolor="blue",
                    edgecolor="k",
                    lw=0.5,
                )
            )
        ax.set_ylim(botm[-1], top)
        ax.set_xlim(0, 1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_title(f"t={times[ti]:.0f}", fontsize=9)
        if j > 0:
            ax.set_yticklabels([])
    fig.suptitle("SWI vertical column (single fluid): fresh (blue) / salt (red)")
    fig.tight_layout()
    fig.savefig(ws / "interface.png", dpi=150)
    plt.close("all")


def check_output(idx, test):
    ws = pl.Path(test.workspace)
    times, z = _zeta(test)
    zeta = z[:, 0]  # top cell

    # -- interface stays essentially flat (the small vertical gradient that feeds
    #    the interface-cell storage is negligible at this K)
    spread = float(np.ptp(z, axis=1).max())
    print(f"max interface spread across column = {spread:.4f}")
    assert spread < 0.1, f"interface should stay flat, spread={spread:.4f}"

    for j, (t, zz) in enumerate(zip(times, zeta)):
        print(f"step {j:2d} t={t:6.1f} zeta={zz:8.4f}")

    # -- interface descends through period 1, ascends through period 2
    zeta_p1, zeta_p2 = zeta[:nstp], zeta[nstp:]
    assert np.all(np.diff(zeta_p1) < 0.0), "interface should descend while injecting"
    assert np.all(np.diff(zeta_p2) > 0.0), "interface should ascend while extracting"

    # -- average rate is set by the mass balance (0.1 m/step), and the interface
    #    descends ~2 m so it crosses the -1 and -2 layer faces
    assert abs(np.diff(zeta_p1).mean() + dz_expected) < 1.0e-3, (
        f"mean descent should be {dz_expected} m/step, got {np.diff(zeta_p1).mean()}"
    )
    assert zeta_p1.min() < -2.0, "interface should cross the -1 and -2 faces"

    # -- interface returns to its starting elevation (symmetric inject/extract)
    assert abs(zeta[-1] - zeta0) < 1.0e-3, (
        f"interface should return to {zeta0}, got {zeta[-1]}"
    )

    # -- the model budget must conserve mass
    df = pd.read_csv(ws / "mymodel.bud.csv")
    imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
    peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
    tol = 1.0e-6 + 1.0e-4 * peak
    print(f"max|imbalance|={imbalance:.3e} peak={peak:.3e} tol={tol:.3e}")
    assert imbalance < tol, f"budget does not conserve mass ({imbalance:.3e})"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
    )
    test.run()
