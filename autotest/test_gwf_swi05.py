"""
SWI Package test -- single-fluid (freshwater-only) sharp-interface model of a
flat strip-island aquifer bounded by the sea on both sides.

Conceptual model: a 1-row, 100-column aquifer with a flat top at +50 m and
bottom at -50 m. The sea is represented by a GHB cell at each end of the row
(head fixed at the saltwater head), and freshwater recharge is applied
everywhere. The SWI Package runs in single-fluid mode with a fixed saltwater
head, and the Newton formulation is used. The two cases differ only in the
saltwater head (sea level): 0.0 and 40.0; the latter places much of the island
below sea level.

What check_output verifies: the simulated interface elevation zeta obeys the
Ghyben-Herzberg relation zeta = clip(-40*head + 41*saltwater_head, botm, top) at
every cell (alphaf=40, alphas=41), and the model's volumetric budget closes
(percent discrepancy ~ 0). The zeta output is limited to the cell, so the
clipped relation is what is compared.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = [
    "swi05a",
    "swi05b",
]
saltwater_head = [0.0, 40.0]

Lx = 10000.0  # x position of right edge of model domain
delr = 100.0
delc = 1.0
nlay = 1
nrow = 1
ncol = int(Lx / delr)
top_island = 50.0
botm_aquifer = -50.0
recharge = 0.0001
hydraulic_conductivity = 100.0
ghb_cond_fact = 1.0


def build_models(idx, test):
    ws = test.workspace
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=ws,
        exe_name="mf6",
        memory_print_option="all",
        print_input=True,
    )
    tdis = flopy.mf6.ModflowTdis(sim)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        no_ptcrecord=True,
        outer_maximum=500,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions=True,
    )
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top_island,
        botm=botm_aquifer,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        # alternative_cell_averaging="amt-hmk",  # "harmonic",
        icelltype=1,
        k=hydraulic_conductivity,
    )
    h0 = saltwater_head[idx]
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=zeta_file,
        saltwater_head=h0,
    )
    cghb = hydraulic_conductivity * ghb_cond_fact * delr * delc / 1
    ghb_list = []
    for j in [0, ncol - 1]:
        freshwater_head = h0
        ghb_list.append([0, 0, j, freshwater_head, cghb])

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_list,
    )
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    budget_file = name + ".bud"
    head_file = name + ".hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        budgetcsv_filerecord=name + ".bud.csv",
        head_filerecord=head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    ws = test.workspace
    sim = test.sims[0]
    gwf = sim.gwf[0]
    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    volume_fresh = (head - zeta) * delr
    volume_fresh = volume_fresh.sum()
    pxs = flopy.plot.PlotCrossSection(gwf, line={"row": 0})
    ax = pxs.ax

    colors = ["cyan", "red"]
    pxs.plot_fill_between(zeta, head=head, colors=colors, ax=ax, edgecolors="none")

    pxs.plot_grid()
    title = (
        f"saltwater head = {saltwater_head[idx]}; Fresh volume = {volume_fresh:0.2f}"
    )
    ax.set_title(title)
    ax.set_ylim(botm_aquifer, top_island)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


def check_output(idx, test):
    sim = test.sims[0]
    gwf = sim.gwf[0]
    ws = sim.sim_path
    h0 = saltwater_head[idx]

    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    print(f"head={head}")
    print(f"zeta={zeta}")

    # The interface must satisfy the Ghyben-Herzberg relation (alphaf=40,
    # alphas=41) with the interface elevation constrained to the cell
    # (botm <= zeta <= top).
    zeta_expected = np.clip(-40.0 * head + 41.0 * h0, botm_aquifer, top_island)
    assert np.allclose(zeta, zeta_expected, atol=1.0e-4), (
        f"zeta does not satisfy the Ghyben-Herzberg relation; "
        f"max diff {np.max(np.abs(zeta - zeta_expected))}"
    )

    # the model budget must conserve mass (absolute imbalance relative to peak
    # flow rather than the per-step PERCENT_DIFFERENCE)
    df = pd.read_csv(ws / f"{gwf.name}.bud.csv")
    imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
    peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
    tol = 1.0e-6 + 1.0e-4 * peak
    print(f"{gwf.name}: max|imbalance|={imbalance:.3e} peak={peak:.3e} tol={tol:.3e}")
    assert imbalance < tol, f"budget does not conserve mass (imbalance={imbalance:.3e})"


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
