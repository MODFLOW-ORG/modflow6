"""
Test of the seawater intrusion (SWI) package for steady-state,
one layer, confined, unconfined, and unconfined newton.

"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "swi01a",  # confined
    "swi01b",  # unconfined
    "swi01c",  # unconfined newton
]

top = [0.0, 50.0, 50.0]
newtonoptions = [False, False, "newton"]
icelltype = [0, 1, 1]

Lx = 10000  # meters
delr, delc = 100.0, 1.0
ncol = int(Lx / delr) + 1
nlay = 1
nrow = 1
botm = -400.0
recharge = 0.001
k = 10.0
h0 = 0.0
h1 = h0


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
        outer_maximum=500,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions=newtonoptions[idx],
    )
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top[idx],
        botm=botm,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        alternative_cell_averaging="amt-hmk",  # "harmonic",
        icelltype=icelltype[idx],
        k=k,
    )
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=zeta_file)
    cghb = 10.0 * delr * delc / 10.0
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=[[0, 0, 0, h0, cghb], [0, 0, ncol - 1, h1, cghb]],
    )
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    budget_file = name + ".bud"
    head_file = name + ".hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord=head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def make_cross_section_plot(sim, idx, title):
    import matplotlib.pyplot as plt

    ws = pl.Path(sim.sim_path)
    gwf = sim.gwf[0]
    x = gwf.modelgrid.xcellcenters.flatten()
    botm = gwf.dis.botm.array.flatten()
    fpth = ws / "mymodel.zta"
    head = gwf.output.head().get_data().flatten()
    zeta = flopy.utils.HeadFile(fpth, text="zeta").get_data().flatten()
    if idx == 0:
        tp = np.zeros(head.shape)
    else:
        tp = head
    pxs = flopy.plot.PlotCrossSection(gwf, line={"row": 0})
    ax = pxs.ax
    ax.plot(x, head, "k-")
    ax.plot(x, zeta, "k--")
    ax.fill_between(x, tp, zeta, color="cyan")
    ax.fill_between(x, zeta, botm, color="red")
    ax.set_title(title)
    ax.set_ylim(-400, 50)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


# analytical solutions based on Fetter 1972 paper
def analytical_head_unconfined(x, recharge, a, k, alphaf):
    return np.sqrt(recharge * (a**2 - x**2) / k / (1 + alphaf))


def analytical_head_confined(x, recharge, a, k, alphaf):
    return np.sqrt(recharge * (a**2 - x**2) / k / alphaf)


def get_heads_for_comparison(sim, idx):
    gwf = sim.gwf[0]
    ncol = gwf.dis.ncol.get_data()
    delr = gwf.dis.delr.get_data()[0]
    hk = gwf.npf.k.get_data()[0][0, 0]
    recharge = gwf.get_package("RCH").recharge.get_data()[0][
        0, 0
    ]  # period 0 and cell 0, 0
    jstart = int((ncol - 1) / 2)
    xc = gwf.modelgrid.xcellcenters.flatten()
    Lx = xc[-1] - xc[0]
    xoffset = Lx / 2 + delr / 2
    xanalytical = xc[jstart:] - xoffset
    if idx == 0:
        hanalytical = analytical_head_confined(xanalytical, recharge, Lx / 2.0, hk, 40)
    else:
        hanalytical = analytical_head_unconfined(
            xanalytical, recharge, Lx / 2.0, hk, 40
        )
    hmodel = sim.gwf[0].output.head().get_data().flatten()[jstart:]
    return xanalytical, hanalytical, hmodel


def make_head_plot(sim, idx, title):
    import matplotlib.pyplot as plt

    ws = pl.Path(sim.sim_path)
    xanalytical, hanalytical, hmodel = get_heads_for_comparison(sim, idx)
    ax = plt.subplot(1, 1, 1)
    ax.plot(xanalytical, hanalytical, "b-", label="analytical")
    ax.plot(xanalytical, hmodel, "bo", mfc="none", label="numerical")
    ax.legend()
    ax.set_title(title)
    plt.savefig(ws / "head.png")
    plt.close("all")


def plot_output(idx, test):
    title = {
        0: "Case 1a -- Confined",
        1: "Case 1b -- Unconfined",
        2: "Case 1c -- Unconfined Newton",
    }
    sim = test.sims[0]
    make_cross_section_plot(sim, idx, title[idx])
    make_head_plot(sim, idx, title[idx])


def check_output(idx, test):
    # get the flopy sim object
    sim = test.sims[0]

    xanalytical, hanalytical, hmodel = get_heads_for_comparison(sim, idx)
    diff = hmodel - hanalytical
    print(f"\nSimulation idx {idx}")
    for j in range(hmodel.shape[0]):
        print(j, hanalytical[j], hmodel[j], diff[j])

    if idx in [0, 1]:
        atol = 0.003
        assert np.allclose(diff[:-1], 0, atol=atol), (
            f"Difference between model and analytical solution {diff} "
            f"greater than atol {atol}"
        )
    elif idx == 2:
        # newton
        assert np.all(diff[:-1] < 0.0), (
            f"Difference between model and analytical solution {diff} "
            "should be less than zero for newton"
        )
        assert np.all(diff[:-1] > -0.5), (
            f"Difference between model and analytical solution {diff} "
            "should be greater than tol for newton"
        )


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
