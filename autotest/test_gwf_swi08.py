"""
Test of the seawater intrusion (SWI) package for coupled
freshwater and saltwater model using the SWI-SWI exchange.
Use a simple 1-layer model with 21 columns and 1 row.
Case a is for confined and case b is for unconfined. There
are two transient stress periods.  The first stress period
has fresh groundwater recharge, which causes a freshwater
bubble to form.  There is no freshwater recharge for the
second stress period, which causes the freshwater bubble to
shrink. The test checks that the zeta computed by the SWI
package is correct.

"""

import pathlib as pl

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = [
    "swi08a-conf",
    "swi08b-unconf",
]

ncol = 21
nlay = 1
nrow = 1
delr = np.array([10.0] + 19 * [100.0] + [10.0])
delc = 1.0
botm = -80.0
recharge = {0: 0.0075, 1: 0.0}
k_fw = 10.0
k_sw = 10.0
h0 = 0.0
icelltype = [0, 1]
iconvert = [0, 1]
top = [0.0, 10.0]
newtonoptions = "NEWTON"
ss = 0.0
sy = 0.2


def build_gwf_model(idx, sim, is_saltwater):
    if is_saltwater:
        name = "saltwater"
    else:
        name = "freshwater"

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions=newtonoptions,
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
    strt = 0.0 if is_saltwater else 0.001
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=icelltype[idx],
        k=k_sw if is_saltwater else k_fw,
    )
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=iconvert[idx], ss=ss, sy=sy)
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=zeta_file,
    )
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[0, 0, 0, h0], [0, 0, ncol - 1, h0]],
    )
    if not is_saltwater:
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

    return gwf


def build_models(idx, test):
    ws = test.workspace
    sim_name = "mymodel"
    sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        sim_ws=ws,
        exe_name="mf6",
        memory_print_option="all",
    )

    # transient tdis
    nper = 2
    perioddata = [(80000.0, 100, 1.0), (10000.0, 10, 1.0)]
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)

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
        outer_dvclose=1e-7,
        inner_dvclose=1e-8,
        linear_acceleration="bicgstab",
    )

    gwf_freshwater = build_gwf_model(idx, sim, False)
    gwf_saltwater = build_gwf_model(idx, sim, True)

    swiswi = flopy.mf6.ModflowSwiswi(
        sim,
        exgtype="SWI6-SWI6",
        exgmnamea="freshwater",
        exgmnameb="saltwater",
    )
    sim.register_ims_package(ims, [gwf_freshwater.name, gwf_saltwater.name])

    return sim, None


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    ws = test.workspace
    sim = test.sims[0]
    gwf = sim.gwf[0]
    x = gwf.modelgrid.xcellcenters.flatten()
    fpth = pl.Path(ws) / f"{gwf.name}.zta"
    head = gwf.output.head().get_data().flatten()
    zobj = flopy.utils.HeadFile(fpth, text="zeta")
    times = zobj.times

    ax = plt.subplot(1, 1, 1)
    flopy.plot.PlotCrossSection(gwf, line={"row": 0}, ax=ax)
    for t in times:
        zeta = zobj.get_data(totim=t).flatten()
        ax.plot(x, zeta, "k-")
    ax.plot(x, head, "b-")

    ax.set_ylim(-100, 10.0)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


def check_output(idx, test):
    # get the flopy sim object
    sim = test.sims[0]

    alphaf = 1000.0 / (1025.0 - 1000.0)
    alphas = 1025.0 / (1025.0 - 1000.0)

    # fresh gwf model
    ws = pl.Path(sim.sim_path)
    gwf_fresh = sim.gwf[0]
    fpth = pl.Path(ws) / f"{gwf_fresh.name}.zta"
    head_f = gwf_fresh.output.head().get_data().flatten()
    zeta = flopy.utils.HeadFile(fpth, text="zeta").get_data().flatten()

    # salt gwf model
    gwf_salt = sim.gwf[1]
    fpth = pl.Path(ws) / f"{gwf_salt.name}.zta"
    head_s = gwf_salt.output.head().get_data().flatten()
    zeta_s = flopy.utils.HeadFile(fpth, text="zeta").get_data().flatten()

    zeta_answer = -alphaf * head_f + alphas * head_s
    for j in range(head_f.shape[0]):
        print(j, head_f[j], head_s[j], zeta[j], zeta_s[j], zeta_answer[j])
    assert np.allclose(zeta, zeta_answer), f"zeta is not right {zeta} /= {zeta_answer}"
    assert np.allclose(zeta, zeta_s), (
        f"salt zeta not equal fresh zeta {zeta_s} /= {zeta}"
    )

    # both model budgets must conserve mass. Check the absolute imbalance
    # (TOTAL_IN - TOTAL_OUT) relative to the peak flow rather than the per-step
    # PERCENT_DIFFERENCE, which is meaningless at the many near-zero-flow steps
    # of the quasi-static saltwater model.
    for nam in (gwf_fresh.name, gwf_salt.name):
        df = pd.read_csv(ws / f"{nam}.bud.csv")
        imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
        peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
        tol = 1.0e-6 + 1.0e-4 * peak
        print(f"{nam}: max|imbalance|={imbalance:.3e} peak={peak:.3e} tol={tol:.3e}")
        assert imbalance < tol, (
            f"{nam} budget does not conserve mass (imbalance={imbalance:.3e})"
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
