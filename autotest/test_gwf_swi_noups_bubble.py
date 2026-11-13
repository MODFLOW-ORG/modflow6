"""
Weighting-bias test for the two-fluid SWI formulation: growing and shrinking
freshwater bubble (the swi08 unconfined problem) solved with and without
upstream weighting of the fluid-zone conductance.

A 1-layer, 21-column strip island (bottom -80, unconfined top +10) is
initially all saltwater with both heads at sea level. Recharge in the first
stress period grows a freshwater bubble until the interface reaches the
aquifer bottom and the lens becomes steady; the recharge is removed in the
second period and the bubble shrinks. Freshwater and saltwater are separate
GWF models coupled by a SWI-SWI exchange.

Three cases:
  bubble-upstream : default SWI (Newton, upstream-weighted fluid-zone fraction,
                    cross-storage and cross-flow exchange terms)
  bubble-amtlmk   : DEV_NO_UPSTREAM_WEIGHTING in both models, AMT-LMK averaging
  bubble-hmean    : DEV_NO_UPSTREAM_WEIGHTING in both models, harmonic mean

Without upstream weighting the flow term is Picard and the exchange drops its
cross-flow terms; the storage Newton terms and the cross-storage exchange
terms are retained. All cases use the two-fluid solver recipe (DBD
under-relaxation plus backtracking); without backtracking the AMT-LMK case
diverges when the saltwater zone pinches out at the bottom.

What the test establishes:
  * the AMT-LMK quasi-Newton two-fluid scheme converges and conserves mass
    for both fluids, at a cost of several times the outer iterations of the
    upstream scheme (totals are printed);
  * at the end of period 1 the lens is steady with the interface on the
    aquifer bottom over the interior, and in the near-shore cell the AMT-LMK
    interface is deeper than the upstream one, i.e. closer to the Dupuit
    position (upstream weighting takes the thick interior cell's freshwater
    fraction for the shoreward face, overstating the conductance and
    pulling the interface up); the difference during the shrinking phase is
    reported, as it is a weighting bias with no exact solution to judge it by;
  * the harmonic mean cannot grow a fluid zone from zero thickness (the
    conductance into an all-saltwater cell is zero), so that case is expected
    to fail at the first time step and is marked xfail.
For the no-upstream cases the upstream-weighted reference is built and run
in a sub-workspace so the comparison is self-contained.
"""

import pathlib as pl

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["bubble-upstream", "bubble-amtlmk", "bubble-hmean"]
no_upstream = [False, True, True]
cell_averaging = [None, "AMT-LMK", None]
xfail = [False, False, True]

ncol = 21
delr = 100.0
delc = 1.0
top = 10.0
botm = -80.0
k = 10.0
recharge = {0: 0.0075, 1: 0.0}
sy = 0.2
alphaf = 40.0
alphas = 41.0


def build_gwf(sim, is_salt, noups, avg):
    name = "saltwater" if is_salt else "freshwater"
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0 if is_salt else 0.001)
    npf_kwargs = {}
    if avg is not None:
        npf_kwargs["alternative_cell_averaging"] = avg
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, **npf_kwargs)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=sy)
    swi_kwargs = {}
    if noups:
        swi_kwargs["dev_no_upstream_weighting"] = True
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=f"{name}.zta", **swi_kwargs)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[0, 0, 0, 0.0], [0, 0, ncol - 1, 0.0]]
    )
    if not is_salt:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.bud",
        budgetcsv_filerecord=f"{name}.bud.csv",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def build_sim(ws, noups, avg):
    sim = flopy.mf6.MFSimulation(sim_name="bubble", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim, nper=2, perioddata=[(80000.0, 100, 1.0), (10000.0, 10, 1.0)]
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        csv_outer_output_filerecord="outer.csv",
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
        backtracking_number=20,
        backtracking_tolerance=1.05,
        backtracking_reduction_factor=0.1,
        backtracking_residual_limit=0.002,
        linear_acceleration="bicgstab",
    )
    fresh = build_gwf(sim, False, noups, avg)
    salt = build_gwf(sim, True, noups, avg)
    flopy.mf6.ModflowSwiswi(
        sim, exgtype="SWI6-SWI6", exgmnamea=fresh.name, exgmnameb=salt.name
    )
    sim.register_ims_package(ims, [fresh.name, salt.name])
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, no_upstream[idx], cell_averaging[idx]), None


def read_results(ws):
    """Return (zeta history [ntimes, ncol], fresh head, salt head, outer iters)."""
    ws = pl.Path(ws)
    zobj = flopy.utils.HeadFile(ws / "freshwater.zta", text="zeta")
    zeta = np.array([zobj.get_data(totim=t).flatten() for t in zobj.times])
    hf = flopy.utils.HeadFile(ws / "freshwater.hds").get_data().flatten()
    hs = flopy.utils.HeadFile(ws / "saltwater.hds").get_data().flatten()
    zs = flopy.utils.HeadFile(ws / "saltwater.zta", text="zeta").get_data()
    nouter = len(pd.read_csv(ws / "outer.csv"))
    return zeta, hf, hs, zs.flatten(), nouter


def check_budgets(ws):
    for nam in ("freshwater", "saltwater"):
        df = pd.read_csv(pl.Path(ws) / f"{nam}.bud.csv")
        imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
        peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
        tol = 1.0e-6 + 1.0e-4 * peak
        assert imbalance < tol, f"{nam} budget imbalance {imbalance:.3e} > {tol:.3e}"


def check_output(idx, test):
    if xfail[idx]:
        # the run is expected to fail at the first time step; there is no
        # output to check beyond the framework's failure check
        return
    ws = pl.Path(test.sims[0].sim_path)
    zeta, hf, hs, zs, nouter = read_results(ws)
    print(f"{cases[idx]}: total outer iterations = {nouter}")

    # interface from the two heads, identical in both models; the zeta output
    # is limited to the cell, so compare with the clipped relation
    zgh = np.clip(-alphaf * hf + alphas * hs, botm, top)
    assert np.allclose(zeta[-1], zgh), "zeta violates GH"
    assert np.allclose(zeta[-1], zs), "fresh and salt zeta differ"
    check_budgets(ws)

    # the bubble must have grown: by the end of period 1 the interface is at
    # the bottom over the island interior and above it near the shores
    z_end1 = zeta[99]
    assert z_end1[ncol // 2] < botm + 1.0e-3, "lens did not reach the bottom"
    assert z_end1[1] > botm + 1.0, "no toe near the shore"

    if no_upstream[idx]:
        # run the upstream-weighted reference in a sub-workspace
        ref_ws = ws / "ref-upstream"
        ref = build_sim(ref_ws, False, None)
        ref.set_sim_path(ref_ws)
        ref.exe_name = test.targets["mf6"]
        ref.write_simulation(silent=True)
        success, _ = ref.run_simulation(silent=True)
        assert success, "upstream reference did not run"
        zeta_ref, _, _, _, nouter_ref = read_results(ref_ws)
        print(f"  upstream reference: total outer iterations = {nouter_ref}")
        dz = np.abs(zeta - zeta_ref)
        print(f"  max |zeta - zeta_upstream| at end of period 1 = {dz[99].max():.3f} m")
        print(f"  max |zeta - zeta_upstream| at end of period 2 = {dz[-1].max():.3f} m")
        print(
            f"  center-cell zeta at end: this {zeta[-1][ncol // 2]:.2f} m, "
            f"upstream {zeta_ref[-1][ncol // 2]:.2f} m"
        )
        # steady lens, near-shore cell: the averaged scheme is less biased
        # toward a shallow interface than upstream weighting
        assert z_end1[1] < zeta_ref[99][1] - 1.0, (
            "averaged thickness should give a deeper near-shore interface than "
            f"upstream ({z_end1[1]:.2f} vs {zeta_ref[99][1]:.2f})"
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
        xfail=xfail[idx],
    )
    test.run()
