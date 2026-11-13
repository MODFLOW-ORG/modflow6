"""
Vertical-column test of SWI vertical flow, two-fluid case (coupled freshwater and
saltwater models via the SWI-SWI exchange).

The model is a 1-D vertical column with a freshwater well at
the top of the freshwater model and a fixed saltwater head at the bottom of the
saltwater model. Injecting freshwater at the top pushes the interface down and
displaces saltwater out the bottom boundary; extracting pulls the interface back
up and draws saltwater in. Over the two stress periods the interface marches down
through two layer faces and back up.

Unlike the single-fluid case, the vertical flow here IS buoyancy-restricted: each
fluid moves freely in its buoyant direction and is restricted in the other by the
fluid present at the shared face, so freshwater and saltwater vertical exchange
is controlled by the interface position as it crosses each face.

Checks that the interface descends then ascends, stays flat, crosses the -1 and
-2 faces, that the freshwater and saltwater models agree on the interface
(zeta_f == zeta_s), and that both model budgets conserve mass. Run with --plot to
write the interface-vs-time figure.
"""

import pathlib as pl

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["swi-vcol2"]

nlay, nrow, ncol = 5, 1, 1
delr = delc = 1.0
top = 0.0
botm = [-1.0, -2.0, -3.0, -4.0, -5.0]
sy = 0.2
alphaf = 1000.0 / (1025.0 - 1000.0)  # 40.0
alphas = 1025.0 / (1025.0 - 1000.0)  # 41.0

nstp = 20
q = 0.02
zeta0 = -0.5  # initial (flat) interface elevation


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
    # freshwater hf = 0.0125 -> zeta = -0.5; saltwater hs = 0
    strt = 0.0 if is_salt else -zeta0 / alphaf
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=100.0, save_specific_discharge=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=sy, transient=True)
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    if is_salt:
        # fixed saltwater head at the bottom cell: saltwater exits/enters here
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(nlay - 1, 0, 0), 0.0]])
    else:
        # freshwater well in the top cell: inject then extract
        flopy.mf6.ModflowGwfwel(
            gwf, stress_period_data={0: [[(0, 0, 0), q]], 1: [[(0, 0, 0), -q]]}
        )
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
    flopy.mf6.ModflowTdis(
        sim, nper=2, perioddata=[(nstp, nstp, 1.0), (nstp, nstp, 1.0)]
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        no_ptcrecord=True,
        under_relaxation="DBD",
        under_relaxation_gamma=0.1,
        under_relaxation_theta=0.7,
        under_relaxation_kappa=0.07,
        under_relaxation_momentum=0.0,
        outer_maximum=500,
        inner_maximum=600,
        outer_dvclose=1.0e-7,
        inner_dvclose=1.0e-8,
        linear_acceleration="bicgstab",
    )
    build_gwf(sim, False)
    build_gwf(sim, True)
    flopy.mf6.ModflowSwiswi(
        sim,
        exgtype="SWI6-SWI6",
        exgmnamea="freshwater",
        exgmnameb="saltwater",
    )
    sim.register_ims_package(ims, ["freshwater", "saltwater"])
    return sim, None


def _interface(test):
    """unclamped interface from heads: zeta = -alphaf*hf + alphas*hs."""
    ws = pl.Path(test.workspace)
    hf = flopy.utils.HeadFile(ws / "freshwater.hds")
    hs = flopy.utils.HeadFile(ws / "saltwater.hds")
    times = hf.times
    z = np.array(
        [
            -alphaf * hf.get_data(totim=t).flatten()
            + alphas * hs.get_data(totim=t).flatten()
            for t in times
        ]
    )
    return np.array(times), z  # (ntime, nlay)


def _zeta_output(test, model):
    """clamped per-cell zeta as written to the .zta output file."""
    ws = pl.Path(test.workspace)
    zobj = flopy.utils.HeadFile(ws / f"{model}.zta", text="zeta")
    return np.array([zobj.get_data(totim=t).flatten() for t in zobj.times])


def plot_output(idx, test):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    ws = pl.Path(test.workspace)
    # per-cell (clamped) interface: each cell splits into freshwater (blue,
    # above zeta) and saltwater (red, below zeta) using the .zta output
    zobj = flopy.utils.HeadFile(ws / "freshwater.zta", text="zeta")
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
    fig.suptitle("SWI vertical column (two fluid): fresh (blue) / salt (red)")
    fig.tight_layout()
    fig.savefig(ws / "interface.png", dpi=150)
    plt.close("all")


def check_output(idx, test):
    ws = pl.Path(test.workspace)
    times, zf = _interface(test)
    zeta = zf[:, 0]  # top cell (flat interface)

    # -- freshwater and saltwater models agree on the (clamped) zeta output
    zof = _zeta_output(test, "freshwater")
    zos = _zeta_output(test, "saltwater")
    assert np.allclose(zof, zos, atol=1.0e-6), "fresh and salt zeta must agree"

    # -- interface stays essentially flat
    spread = float(np.ptp(zf, axis=1).max())
    print(f"max interface spread across column = {spread:.4f}")
    assert spread < 0.1, f"interface should stay flat, spread={spread:.4f}"

    for j, (t, zz) in enumerate(zip(times, zeta)):
        print(f"step {j:2d} t={t:6.1f} zeta={zz:8.4f}")

    # -- interface descends through period 1, ascends through period 2
    zeta_p1, zeta_p2 = zeta[:nstp], zeta[nstp:]
    assert np.all(np.diff(zeta_p1) < 0.0), "interface should descend while injecting"
    assert np.all(np.diff(zeta_p2) > 0.0), "interface should ascend while extracting"

    # -- interface descends ~2 m (crossing the -1 and -2 faces) and returns near
    #    its start (small asymmetry from the two-fluid coupling is acceptable)
    assert zeta_p1.min() < -2.0, "interface should cross the -1 and -2 faces"
    assert abs(zeta[-1] - zeta0) < 0.1, (
        f"interface should return near {zeta0}, got {zeta[-1]}"
    )

    # -- quantitative vertical-flow rate: injecting freshwater volume q*t at the
    #    top displaces the sharp interface downward as saltwater is pushed out the
    #    bottom, and every unit of injected volume sweeps 1/(area*sy) of interface
    #    descent. So the cumulative descent over the injection period must equal
    #    the analytical q/(area*sy)*perlen. This checks the buoyancy-restricted
    #    vertical saltwater conductance quantitatively (a wrong vertical
    #    conductance would move the interface at the wrong speed while still
    #    conserving mass). The small excess is the interface smoothing as zeta
    #    crosses the two layer faces.
    area = delr * delc
    perlen = float(nstp)  # perioddata uses perlen == nstp with tsmult 1.0
    descent = zeta0 - float(zeta_p1[-1])
    descent_expected = q / (area * sy) * perlen
    print(f"interface descent: observed={descent:.4f} expected={descent_expected:.4f}")
    assert abs(descent - descent_expected) < 0.05 * descent_expected, (
        f"interface descent {descent:.4f} should match q/(area*sy)*perlen "
        f"= {descent_expected:.4f} (vertical saltwater flow rate)"
    )

    # -- both model budgets must conserve mass
    for nam in ("freshwater", "saltwater"):
        df = pd.read_csv(ws / f"{nam}.bud.csv")
        imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
        peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
        tol = 1.0e-6 + 1.0e-4 * peak
        print(f"{nam}: max|imbalance|={imbalance:.3e} peak={peak:.3e} tol={tol:.3e}")
        assert imbalance < tol, f"{nam} budget does not conserve mass ({imbalance:.3e})"


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
