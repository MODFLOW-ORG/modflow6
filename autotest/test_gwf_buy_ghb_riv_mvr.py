"""BUY must be included in the GHB and RIV water provided to MVR.

Two cells: a CHD cell supplies water that discharges through a GHB or RIV in
the second cell, and MVR routes a fraction of it to an SFR reach. The water
leaving the second cell must equal the flow into it from the CHD cell, and MVR
must receive the specified fraction of that flow, with and without a density
contrast between the aquifer and the boundary.

When the second cell is also a constant-head cell, the boundary has no flow,
so it must not provide any water to MVR.
"""

import itertools

import flopy
import numpy as np
import pytest
from framework import TestFramework

# kind, concentration, mover fraction, boundary in a constant-head cell
settings = [
    (*s, False) for s in itertools.product(["GHB", "RIV"], [0.0, 16.0], [0.4, 1.0])
] + [("GHB", 0.0, 1.0, True), ("RIV", 0.0, 1.0, True)]
cases = [f"buyghbrivmvr{i}" for i in range(len(settings))]


def build_models(idx, test):
    kind, concentration, fraction, constant_head = settings[idx]
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="flow", save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(ims, ["flow"])
    grid = dict(nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfdis(gwf, **grid)
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    chd = [[(0, 0, 0), 2.0, concentration]]
    if constant_head:
        chd.append([(0, 0, 1), 1.5, concentration])
    flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD",
        auxiliary=["CONCENTRATION"],
        stress_period_data=chd,
    )
    # Freshwater boundary in the second cell; the RIV stays above its bottom.
    if kind == "GHB":
        flopy.mf6.ModflowGwfghb(
            gwf,
            pname="BND",
            mover=True,
            auxiliary=["CONCENTRATION"],
            stress_period_data=[[(0, 0, 1), 1.0, 10.0, 0.0]],
        )
    else:
        flopy.mf6.ModflowGwfriv(
            gwf,
            pname="BND",
            mover=True,
            auxiliary=["CONCENTRATION"],
            stress_period_data=[[(0, 0, 1), 1.0, 10.0, 0.5, 0.0]],
        )
    flopy.mf6.ModflowGwfbuy(gwf, packagedata=[(0, 1.3, 0.0, "trans", "CONCENTRATION")])
    # SFR reaches with zero streambed conductivity only receive mover water.
    flopy.mf6.ModflowGwfsfr(
        gwf,
        pname="SFR",
        mover=True,
        nreaches=2,
        packagedata=[
            (0, (0, 0, 0), 1.0, 1.0, 0.001, 2.0, 1.0, 0.0, 0.04, 1, 1.0, 0),
            (1, (0, 0, 1), 1.0, 1.0, 0.001, 2.0, 1.0, 0.0, 0.04, 1, 1.0, 0),
        ],
        connectiondata=[(0, -1), (1, 0)],
        budget_filerecord="sfr.cbc",
    )
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=1,
        maxpackages=2,
        packages=[("BND",), ("SFR",)],
        perioddata=[("BND", 0, "SFR", 0, "FACTOR", fraction)],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="flow.hds",
        budget_filerecord="flow.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imst = flopy.mf6.ModflowIms(
        sim, filename="trans.ims", linear_acceleration="BICGSTAB"
    )
    sim.register_ims_package(imst, ["trans"])
    flopy.mf6.ModflowGwtdis(gwt, **grid)
    flopy.mf6.ModflowGwtic(gwt, strt=concentration)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[("CHD", "AUX", "CONCENTRATION"), ("BND", "AUX", "CONCENTRATION")],
    )
    flopy.mf6.ModflowGwtsft(
        gwt,
        pname="SFT",
        flow_package_name="SFR",
        packagedata=[(0, concentration), (1, concentration)],
    )
    flopy.mf6.ModflowGwtmvt(gwt)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="trans.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea="flow", exgmnameb="trans")
    return sim, None


def check_output(idx, test):
    kind, concentration, fraction, constant_head = settings[idx]
    ws = test.workspace
    h = flopy.utils.HeadFile(ws / "flow.hds", precision="double").get_data().ravel()
    budget = flopy.utils.CellBudgetFile(ws / "flow.cbc", precision="double")
    qbnd = budget.get_data(text=kind)[0]["q"].item()
    qtomvr = budget.get_data(text=f"{kind}-TO-MVR")[0]["q"].item()
    # CSR row 1 is [diagonal, neighbor 0]; positive is flow into cell 1.
    qface = budget.get_data(text="FLOW-JA-FACE")[0].ravel()[3]
    sfr = flopy.utils.CellBudgetFile(ws / "sfr.cbc", precision="double")
    qfrommvr = sfr.get_data(text="FROM-MVR")[0]["q"].sum()

    assert qface > 0.0
    if constant_head:
        # A boundary in a constant-head cell has no flow and provides no water.
        np.testing.assert_allclose([qbnd, qtomvr, qfrommvr], 0.0, rtol=0.0, atol=1e-10)
        return
    # All water entering cell 1 leaves through the boundary.
    np.testing.assert_allclose(qbnd + qtomvr, -qface, rtol=0.0, atol=1e-8)
    # MVR receives the specified fraction of the boundary discharge.
    np.testing.assert_allclose(qtomvr, -fraction * qface, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(qfrommvr, fraction * qface, rtol=0.0, atol=1e-8)
    # With a density contrast, the discharge differs from cond * (h - stage),
    # so the checks above detect a mover flow without the BUY terms.
    if concentration > 0.0:
        assert abs(10.0 * (h[1] - 1.0) - qface) > 1e-4


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
