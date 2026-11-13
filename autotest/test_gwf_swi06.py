"""
SWI Package test -- single-fluid (freshwater-only) sharp-interface model of a
flat strip-island aquifer under a transient, rising sea level.

Conceptual model: the same flat strip-island geometry as swi05 (1 row, 100
columns, top +50 m, bottom -50 m, a GHB cell at each end, recharge everywhere),
but run transient over 10 stress periods of 3,650 days each (100 years total).
Sea level rises linearly from 0 to 40 m. The time-varying saltwater head is
supplied through the TVG (Time-Varying Grid) Utility of the SWI Package as a
SALTWATER_HEAD auxiliary variable, while the GHB heads track the same rising sea
level. The Newton formulation is used.

What check_output verifies: at every output time the simulated interface
elevation zeta obeys the Ghyben-Herzberg relation
zeta = clip(-40*head + 41*saltwater_head(t), botm, top) using that period's sea
level (alphaf=40, alphas=41), and the model's volumetric budget closes for all
periods. The zeta output is limited to the cell, so the clipped relation is
what is compared.
This case also exercises the TVG Utility and transient interface movement.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = [
    "swi06a",
]

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
perlen = 10 * 365.0
nper = 10
perioddata = nper * [(perlen, 1, 1.0)]
sealevel_start = 0.0
sealevel_stop = 40.0


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
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
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
        icelltype=1,
        k=hydraulic_conductivity,
    )
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=zeta_file,
        tvg6_filerecord=f"{name}.swi.tvg",
    )

    # initialize the TVG Utility with the saltwater head
    sl = np.linspace(sealevel_start, sealevel_stop, nper)
    aux = {}
    for kper in range(nper):
        # this list has 2 entries, one for hsalt and one for myauxvar
        aux[kper] = [sl[kper], np.ones((nlay, nrow, ncol))]
    swi.tvg.initialize(auxiliary=["saltwater_head", "myauxvar"], aux=aux)

    cghb = hydraulic_conductivity * ghb_cond_fact * delr * delc / 1
    ghbspd = {}
    for kper in range(nper):
        ghb_list = []
        for j in [0, ncol - 1]:
            freshwater_head = sl[kper]
            ghb_list.append([0, 0, j, freshwater_head, cghb])
        ghbspd[kper] = ghb_list

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
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

    sim = test.sims[0]
    ws = test.workspace
    gwf = sim.gwf[0]
    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    volume_fresh = (head - zeta) * delr
    volume_fresh = volume_fresh.sum()
    pxs = flopy.plot.PlotCrossSection(gwf, line={"row": 0})
    ax = pxs.ax

    zobj = gwf.swi.output.zeta()
    times = zobj.times
    for t in times:
        zeta = zobj.get_data(totim=t)
        pxs.plot_surface(zeta)

    pxs.plot_grid()
    title = "Rising Sea Level"
    ax.set_title(title)
    ax.set_ylim(botm_aquifer, top_island)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


def check_output(idx, test):
    sim = test.sims[0]
    gwf = sim.gwf[0]
    ws = sim.sim_path

    # time-varying saltwater head (sea level) imposed through the TVG Utility
    sl = np.linspace(sealevel_start, sealevel_stop, nper)
    hobj = gwf.output.head()
    zobj = gwf.swi.output.zeta()
    times = hobj.times
    assert len(times) == nper, f"expected {nper} output times, got {len(times)}"

    # The interface must satisfy the Ghyben-Herzberg relation (alphaf=40,
    # alphas=41) at every time, using the period's saltwater head, with the
    # interface elevation constrained to the cell (botm <= zeta <= top).
    for i, t in enumerate(times):
        head = hobj.get_data(totim=t).flatten()
        zeta = zobj.get_data(totim=t).flatten()
        zeta_expected = np.clip(-40.0 * head + 41.0 * sl[i], botm_aquifer, top_island)
        adiff = np.max(np.abs(zeta - zeta_expected))
        print(f"period {i}: hsalt={sl[i]:.3f} max|zeta-exp|={adiff:.2e}")
        assert np.allclose(zeta, zeta_expected, atol=1.0e-4), (
            f"period {i}: zeta does not satisfy the Ghyben-Herzberg relation; "
            f"max diff {adiff}"
        )

    # the model budget must conserve mass for all periods (absolute imbalance
    # relative to peak flow rather than the per-step PERCENT_DIFFERENCE)
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
