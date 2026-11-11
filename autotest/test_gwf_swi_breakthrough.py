"""
Premature-breakthrough test for the two-fluid SWI vertical flow.

A 1-D vertical column of two coupled models (freshwater over saltwater) is driven by
recharge on the freshwater top; the freshwater front descends as a sharp,
storage-controlled front. The physical requirement -- the point of the vertical buoyancy
restriction -- is that freshwater does NOT show up in a cell before the interface in
the OVERLYING cell has reached that cell's bottom. In other words the profile of
freshwater saturated cell fraction must be a sharp descending front: a cell may only
acquire freshwater once the cell above it is full. Freshwater "leaking" ahead of the
front (a deep cell wetting while a shallower cell is not yet full) is premature
breakthrough -- an artifact of head-driven downward freshwater flow into the
saltwater zone, which the buoyancy restriction suppresses.

Without the restriction, a small amount of freshwater moves one cell ahead of the
front (Sf~0.04 while the overlying cell is ~0.46); with the restriction the cell below
the front stays dry until the front cell fills. This test asserts the restricted
(correct) behavior and would fail if the vertical flow were unrestricted. It also
checks that the front descends at the mass-balance rate and that both budgets
conserve mass.

Run with --plot for the per-cell fresh(blue)/salt(red) filmstrip.
"""

import pathlib as pl

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["swi-breakthrough"]

nlay, nrow, ncol = 10, 1, 1
delr = delc = 1.0
top = 0.0
botm = [-float(k + 1) for k in range(nlay)]  # 0 .. -10, 1 m cells
sy = 0.2
recharge = 0.02
alphaf = 1000.0 / (1025.0 - 1000.0)  # 40.0
alphas = 1025.0 / (1025.0 - 1000.0)  # 41.0
zeta0 = -0.5  # initial thin freshwater lens in the top cell

nstp = 20
perlen = 40.0  # -> mass-balance front depth = recharge*perlen/sy = 4 m below zeta0

cell_top = np.array([top] + botm[:-1])
cell_bot = np.array(botm)


def build_gwf(sim, is_salt):
    name = "saltwater" if is_salt else "freshwater"
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0 if is_salt else -zeta0 / alphaf)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=sy, transient=True)
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    if is_salt:
        # saltwater exits the bottom as the interface descends
        flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[[nlay - 1, 0, 0, 0.0, 10.0]])
    else:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budgetcsv_filerecord=name + ".bud.csv",
        saverecord=[("HEAD", "ALL")],
    )
    return gwf


def build_models(idx, test):
    sim = flopy.mf6.MFSimulation(
        sim_name="mymodel", sim_ws=test.workspace, exe_name="mf6"
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(perlen, nstp, 1.0)])
    # COMPLEX solver settings: the freshwater model has no head boundary
    # (recharge only); its head level is set through the exchange with the
    # saltwater model, which has the GHB at the bottom
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        complexity="complex",
        outer_maximum=200,
        inner_maximum=200,
        outer_dvclose=1.0e-7,
        inner_dvclose=1.0e-8,
        linear_acceleration="bicgstab",
    )
    build_gwf(sim, False)
    build_gwf(sim, True)
    flopy.mf6.ModflowSwiswi(
        sim, exgtype="SWI6-SWI6", exgmnamea="freshwater", exgmnameb="saltwater"
    )
    return sim, None


def _saturation(test):
    """per-cell freshwater saturated cell fraction over time, from the coupled heads."""
    ws = pl.Path(test.workspace)
    hf = flopy.utils.HeadFile(ws / "freshwater.hds")
    hs = flopy.utils.HeadFile(ws / "saltwater.hds")
    times = hf.times
    sf = []
    for t in times:
        zeta = (
            -alphaf * hf.get_data(totim=t).flatten()
            + alphas * hs.get_data(totim=t).flatten()
        )
        sf.append(np.clip((cell_top - zeta) / (cell_top - cell_bot), 0.0, 1.0))
    return np.array(times), np.array(sf)  # (ntime, nlay)


def plot_output(idx, test):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    ws = pl.Path(test.workspace)
    times, sf = _saturation(test)
    snap = np.linspace(0, len(times) - 1, 9).round().astype(int)
    fig = plt.figure(figsize=(12, 4))
    for j, ti in enumerate(snap):
        ax = fig.add_subplot(1, len(snap), j + 1)
        for k in range(nlay):
            zc = cell_top[k] - sf[ti, k] * (cell_top[k] - cell_bot[k])
            ax.add_patch(
                Rectangle(
                    (0, cell_bot[k]),
                    1.0,
                    zc - cell_bot[k],
                    facecolor="red",
                    edgecolor="k",
                    lw=0.4,
                )
            )
            ax.add_patch(
                Rectangle(
                    (0, zc),
                    1.0,
                    cell_top[k] - zc,
                    facecolor="blue",
                    edgecolor="k",
                    lw=0.4,
                )
            )
        ax.set_ylim(botm[-1], top)
        ax.set_xlim(0, 1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_title(f"t={times[ti]:.0f}", fontsize=9)
        if j > 0:
            ax.set_yticklabels([])
    fig.suptitle("SWI premature-breakthrough test: sharp descending front")
    fig.tight_layout()
    fig.savefig(ws / "breakthrough.png", dpi=150)
    plt.close("all")


def check_output(idx, test):
    ws = pl.Path(test.workspace)
    times, sf = _saturation(test)

    # -- NO PREMATURE BREAKTHROUGH: a cell may only contain freshwater once the cell
    #    above it is full. Flag any cell that is wet while its overlying cell is
    #    not (freshwater leaked ahead of the storage front).
    wet, full = 0.02, 0.95
    worst = 0.0
    for it, t in enumerate(times):
        for k in range(1, nlay):
            if sf[it, k] > wet and sf[it, k - 1] < full:
                worst = max(worst, sf[it, k])
                print(
                    f"  premature at t={t:.1f} layer {k + 1}: "
                    f"Sf={sf[it, k]:.3f} while overlying Sf={sf[it, k - 1]:.3f}"
                )
    assert worst == 0.0, (
        f"premature breakthrough: freshwater reached a cell (max Sf={worst:.3f}) "
        "before its overlying cell filled"
    )

    # -- sanity: the front descends and reaches the mass-balance depth
    front = np.array(
        [cell_bot[np.where(s > 0.5)[0][-1]] if np.any(s > 0.5) else top for s in sf]
    )
    mb_front = zeta0 - recharge * perlen / sy  # -4.5
    print(f"final front={front[-1]:.1f}  mass-balance front={mb_front:.1f}")
    assert np.all(np.diff(front) <= 1.0e-6), "front should descend monotonically"
    assert abs(front[-1] - mb_front) <= 1.0, "front should track mass balance"

    # -- both budgets conserve mass
    for nam in ("freshwater", "saltwater"):
        df = pd.read_csv(ws / f"{nam}.bud.csv")
        imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
        peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
        tol = 1.0e-6 + 1.0e-4 * peak
        print(f"{nam}: max|imbalance|={imbalance:.3e} tol={tol:.3e}")
        assert imbalance < tol, f"{nam} budget does not conserve mass"


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
