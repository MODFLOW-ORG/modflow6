"""
Weighting-bias test for the SWI horizontal conductance: steady freshwater lens
beneath a strip island, single-fluid mode, compared with the exact Dupuit
solution.

A 1-row, 100-column unconfined aquifer (top +20, bottom -400) receives uniform
recharge N and is bounded by constant freshwater head h_e at the two end cells
(sea level is zero, so the saltwater head is zero and zeta = -alphaf*h). The
interface stays well above the aquifer bottom everywhere, so the freshwater
thickness is b = (1 + alphaf) h and the steady Dupuit balance

    d/dx [ K (1 + alphaf) h dh/dx ] = -N

has the exact solution h^2 = h_e^2 + N (L^2 - x^2) / (K (1 + alphaf)), with x
measured from the island center and L the distance to the boundary cell
centers. Because h^2 is quadratic in x, a scheme that uses the arithmetic mean
of the two cell thicknesses reproduces this solution exactly at the cell
centers; upstream weighting of the thickness overestimates the conductance on
the thinning side and biases the heads low.

Three cases solve the same problem:
  island-upstream : default SWI (Newton, upstream-weighted fluid-zone fraction)
  island-amtlmk   : DEV_NO_UPSTREAM_WEIGHTING with ALTERNATIVE_CELL_AVERAGING
                    AMT-LMK (arithmetic-mean thickness)
  island-hmean    : DEV_NO_UPSTREAM_WEIGHTING with the default harmonic mean

Checks: every case converges and conserves mass, zeta obeys Ghyben-Herzberg,
the AMT-LMK case matches the exact solution to solver tolerance, and the
upstream case is measurably worse than AMT-LMK (the bias this option exists to
quantify). The harmonic-mean error is reported for information.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["island-upstream", "island-amtlmk", "island-hmean"]
no_upstream = [False, True, True]
cell_averaging = [None, "AMT-LMK", None]

Lx = 10000.0
delr = 100.0
delc = 1.0
ncol = int(Lx / delr)
top = 20.0
botm = -400.0
k = 10.0
recharge = 0.001
h_edge = 0.5
alphaf = 40.0
alphas = 41.0


def analytical_head():
    x = (np.arange(ncol) + 0.5) * delr - 0.5 * Lx
    L = 0.5 * Lx - 0.5 * delr
    h2 = h_edge**2 + recharge * (L**2 - x**2) / (k * (1.0 + alphaf))
    return np.sqrt(h2)


def build_models(idx, test):
    ws = test.workspace
    name = "island"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        no_ptcrecord=True,
        outer_maximum=500,
        inner_maximum=300,
        outer_dvclose=1.0e-9,
        inner_dvclose=1.0e-10,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=h_edge)
    npf_kwargs = {}
    if cell_averaging[idx] is not None:
        npf_kwargs["alternative_cell_averaging"] = cell_averaging[idx]
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, **npf_kwargs)
    swi_kwargs = {}
    if no_upstream[idx]:
        swi_kwargs["dev_no_upstream_weighting"] = True
    flopy.mf6.ModflowGwfswi(
        gwf, zeta_filerecord=f"{name}.zta", saltwater_head=0.0, **swi_kwargs
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[0, 0, 0, h_edge], [0, 0, ncol - 1, h_edge]]
    )
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.bud",
        budgetcsv_filerecord=f"{name}.bud.csv",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, None


def head_error(test):
    sim = test.sims[0]
    gwf = sim.gwf[0]
    head = gwf.output.head().get_data().flatten()
    return head, np.abs(head - analytical_head()).max()


def check_output(idx, test):
    sim = test.sims[0]
    gwf = sim.gwf[0]
    ws = sim.sim_path
    head, err = head_error(test)
    zeta = gwf.swi.output.zeta().get_data().flatten()
    print(f"{cases[idx]}: max |h - h_exact| = {err:.4e} m")

    # Ghyben-Herzberg with a prescribed saltwater head of zero
    zeta_expected = np.clip(-alphaf * head, botm, top)
    assert np.allclose(zeta, zeta_expected, atol=1.0e-4), "zeta violates GH"
    assert zeta.min() > botm + 1.0, "interface reached the bottom; lens not simple"

    # mass balance
    df = pd.read_csv(ws / f"{gwf.name}.bud.csv")
    imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
    peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
    assert imbalance < 1.0e-6 + 1.0e-5 * peak, f"budget imbalance {imbalance:.3e}"

    if cases[idx] == "island-amtlmk":
        # arithmetic-mean thickness reproduces the quadratic h^2 exactly
        assert err < 1.0e-3, f"AMT-LMK head error {err:.3e} exceeds 1e-3 m"
    elif cases[idx] == "island-upstream":
        # upstream weighting biases the lens; it must be visibly worse than the
        # unbiased scheme (which is exact), and the bias must be toward low heads
        assert err > 1.0e-2, f"upstream error {err:.3e} unexpectedly small"
        assert (head - analytical_head()).max() < 1.0e-3, (
            "upstream weighting should bias heads low, not high"
        )


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
