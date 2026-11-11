"""
Checks that the SWI Package gives the same heads as a plain NPF simulation when
the interface (zeta) lies well below the model bottom. With no saltwater in the
domain SWI must have no effect, and the test exercises the horizontal and
vertical SWI conductance in the all-freshwater limit -- the freshwater zone is
the full cell (S^s = 0) and vertical flow is unrestricted.

To stress the conductance values, the model is confined and multi-layer with
heterogeneous horizontal (k) and vertical (k33) hydraulic conductivity, and the
stresses (recharge + top-corner drains, plus an injection well in the bottom
layer) produce both upward and downward vertical flow. The single-fluid SWI
package uses the default saltwater head of zero, so zeta = -alphaf*hf is on the
order of -300, far below the model bottom. The same model is run with and
without SWI and the heads must match; a self-check confirms the flow field
actually contains both vertical flow directions.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["swi-npf"]

nlay, nrow, ncol = 3, 1, 10

# heterogeneous, reproducible conductivity fields (independent k and k33)
_rng = np.random.default_rng(2024)
k = _rng.uniform(1.0, 100.0, size=(nlay, nrow, ncol))
k33 = _rng.uniform(0.01, 5.0, size=(nlay, nrow, ncol))


def build_sim(ws, exe, with_swi):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim)  # single steady-state period
    flopy.mf6.ModflowIms(
        sim,
        print_option="all",
        linear_acceleration="bicgstab",
        outer_maximum=100,
        inner_maximum=100,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=1.0,
        top=10.0,
        botm=[5.0, 0.0, -5.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=k, k33=k33, save_specific_discharge=True
    )
    # -- stresses that produce both downward and upward vertical flow:
    #    recharge on the top layer + drains at the top corners (downward), and an
    #    injection well in the bottom layer (upward above the well)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.02)
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(nlay - 1, 0, ncol // 2), 100.0]])
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 8.0], [(0, 0, ncol - 1), 8.0]]
    )
    if with_swi:
        # default saltwater head is zero -> zeta = -alphaf*hf << model bottom
        flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budget_filerecord=name + ".bud",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, test.targets["mf6"], with_swi=True), None


def check_output(idx, test):
    gwf_swi = test.sims[0].gwf[0]

    # -- self-check: the flow field must have both upward and downward vertical
    #    flow, otherwise this is not a meaningful test of the vertical conductance
    spdis = gwf_swi.output.budget().get_data(text="DATA-SPDIS")[0]
    qz = spdis["qz"]
    print(f"qz range: [{qz.min():.3e}, {qz.max():.3e}]")
    assert qz.min() < 0.0 and qz.max() > 0.0, (
        "test should exercise both upward and downward vertical flow "
        f"(qz range [{qz.min():.3e}, {qz.max():.3e}])"
    )

    # -- build and run the equivalent NPF-only (no SWI) model in a sibling dir
    npf_ws = test.workspace / "npf"
    sim_npf = build_sim(npf_ws, test.targets["mf6"], with_swi=False)
    sim_npf.write_simulation(silent=True)
    success, _ = sim_npf.run_simulation(silent=True)
    assert success, "NPF reference model failed to run"

    # -- heads must match: SWI is a no-op when zeta is below the model bottom
    h_swi = gwf_swi.output.head().get_data()
    h_npf = sim_npf.gwf[0].output.head().get_data()
    dmax = np.abs(h_swi - h_npf).max()
    print(f"max |h_swi - h_npf| = {dmax:.3e}")
    assert np.allclose(h_swi, h_npf, atol=1.0e-9, rtol=0.0), (
        f"SWI must reproduce NPF when zeta is below the model bottom "
        f"(max head difference {dmax:.3e})"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
