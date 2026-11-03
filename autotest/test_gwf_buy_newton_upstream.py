"""
Test the upstream-saturation treatment of the BUY package with the Newton
formulation.

With the Newton formulation the conductance of a horizontal connection is
the saturated conductance multiplied by the saturation of the upstream cell.
The buoyancy flow across a connection between cells of different density is
generally not zero when the heads of the two cells are equal, so the upstream
cell must be selected from the direction of the total (head plus buoyancy)
flow rather than from the head difference.  Otherwise the flow, and
therefore the residual, jumps when the head ordering of the two cells
changes, and the solution can cycle between two states without converging.

The model has three cells:

    A  thin, dry, fresh cell (top 21, bottom 14)
    B  thick, saline cell (top 21, bottom -60), constant head 0.36
    C  fresh cell below A (top 14, bottom -60), constant head 0.0

A and B share a horizontal connection; A drains to C through a weak vertical
connection, its only outlet.  The buoyancy term drives saline water from B
into A, so at the solution the head in A is above the head in B while the
flow is still from B into A.  When the upstream cell was selected from the
head difference, A (with zero saturation) became upstream as soon as its
head rose above the head in B, the flow dropped to zero, and the head in A
flipped between two values on every outer iteration.

The model is run as a single DISU model and as two DIS models joined by a
GWF-GWF exchange (which uses the interface model).  Both must converge to
the head in A computed independently below, and their budgets must close.
The same model with DRHODC = 0 must give the same heads as the model
without the BUY package.
"""

import re
from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["buynwt_single", "buynwt_exg"]

# geometry and properties
top = 21.0
bot_a, bot_b, bot_c = 14.0, -60.0, -60.0
delr = delc = 100.0
area = delr * delc
cl = 50.0  # distance from the cell centers of A and B to their shared face
kh = 10.0
kv_a = 1.0e-4  # weak vertical outlet from A to C
h_b, h_c = 0.36, 0.0

# density
denseref = 62.428
drhodc = 0.0446
c_a, c_b, c_c = 0.0, 35.0, 0.0


def expected_head_a():
    """Head in cell A from the Newton formulation of NPF and BUY with the
    upstream cell selected from the total flow direction.  B is upstream
    and A is dry, so the elevation of A is its bottom."""
    t_a, t_b = kh * (top - bot_a), kh * (top - bot_b)
    condsat = delc * t_a * t_b / (t_a * cl + t_b * cl)
    sat_b = (h_b - bot_b) / (top - bot_b)
    elev_a = bot_a
    elev_b = bot_b + 0.5 * sat_b * (top - bot_b)
    rho_a = 1.0 + drhodc * c_a / denseref
    rho_b = 1.0 + drhodc * c_b / denseref
    wt = 0.5
    rhoterm = wt * rho_a + (1.0 - wt) * rho_b
    elevnm = (1.0 - wt) * elev_a + wt * elev_b
    cv = area / (0.5 * (top - bot_a) / kv_a + 0.5 * (bot_a - bot_c) / kh)

    def residual(h_a):
        gbuy = (rhoterm - 1.0) * (h_b - h_a) + (rho_b - rho_a) * (
            (1.0 - wt) * h_a + wt * h_b - elevnm
        )
        return sat_b * condsat * (h_b - h_a + gbuy) - cv * (h_a - h_c)

    # the residual is linear in h_a
    r0, r1 = residual(0.0), residual(1.0)
    return -r0 / (r1 - r0)


def add_gwt(sim, name, dispkg, diskw, conc, cells):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name)
    dispkg(gwt, **diskw)
    flopy.mf6.ModflowGwtic(gwt, strt=conc)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtcnc(gwt, stress_period_data=list(zip(cells, conc)))
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtoc(gwt)


def add_gwf(sim, name, dispkg, diskw, k33, chd, rhodc, buy):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
    dispkg(gwf, **diskw)
    flopy.mf6.ModflowGwfic(gwf, strt=0.4)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kh, k33=k33)
    if buy:
        flopy.mf6.ModflowGwfbuy(
            gwf,
            denseref=denseref,
            packagedata=[(0, rhodc, 0.0, f"t{name}", "CONCENTRATION")],
        )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )


def build_sim(ws, exe, layout, rhodc=drhodc, buy=True):
    sim = flopy.mf6.MFSimulation(sim_name="buynwt", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    ims_f = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_maximum=50,
        outer_dvclose=1e-8,
        inner_maximum=50,
        inner_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
    )
    ims_t = flopy.mf6.ModflowIms(
        sim,
        filename="tran.ims",
        outer_maximum=5,
        inner_maximum=10,
        outer_dvclose=1e-6,
        inner_dvclose=1e-9,
        linear_acceleration="BICGSTAB",
    )
    if layout == "single":
        # nodes: 0 = A, 1 = C, 2 = B
        diskw = dict(
            nodes=3,
            nja=7,
            top=[top, bot_a, top],
            bot=[bot_a, bot_c, bot_b],
            area=area,
            iac=[3, 2, 2],
            ja=[0, 1, 2, 1, 0, 2, 0],
            ihc=[1, 0, 1, 1, 0, 1, 1],
            cl12=[0.0, 0.5 * (top - bot_a), cl, 0.0, 0.5 * (bot_a - bot_c), 0.0, cl],
            hwva=[0.0, area, delc, 0.0, area, 0.0, delc],
        )
        cells = [(0,), (1,), (2,)]
        add_gwf(
            sim,
            "gwf",
            flopy.mf6.ModflowGwfdisu,
            diskw,
            [kv_a, kh, kh],
            [((1,), h_c), ((2,), h_b)],
            rhodc,
            buy,
        )
        add_gwt(sim, "tgwf", flopy.mf6.ModflowGwtdisu, diskw, [c_a, c_c, c_b], cells)
        gwfs = ["gwf"]
    else:
        da = dict(
            nlay=2,
            nrow=1,
            ncol=1,
            delr=delr,
            delc=delc,
            top=top,
            botm=[bot_a, bot_c],
        )
        db = dict(nlay=1, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=bot_b)
        add_gwf(
            sim,
            "gwfa",
            flopy.mf6.ModflowGwfdis,
            da,
            [kv_a, kh],
            [((1, 0, 0), h_c)],
            rhodc,
            buy,
        )
        add_gwf(
            sim, "gwfb", flopy.mf6.ModflowGwfdis, db, kh, [((0, 0, 0), h_b)], rhodc, buy
        )
        add_gwt(
            sim,
            "tgwfa",
            flopy.mf6.ModflowGwtdis,
            da,
            [c_a, c_c],
            [(0, 0, 0), (1, 0, 0)],
        )
        add_gwt(sim, "tgwfb", flopy.mf6.ModflowGwtdis, db, [c_b], [(0, 0, 0)])
        flopy.mf6.ModflowGwfgwf(
            sim,
            exgmnamea="gwfa",
            exgmnameb="gwfb",
            newton=True,
            nexg=1,
            exchangedata=[((0, 0, 0), (0, 0, 0), 1, cl, cl, delc)],
        )
        gwfs = ["gwfa", "gwfb"]
    for name in gwfs:
        flopy.mf6.ModflowGwfgwt(
            sim, exgmnamea=name, exgmnameb=f"t{name}", filename=f"{name}.gwfgwt"
        )
    sim.register_solution_package(ims_f, gwfs)
    sim.register_solution_package(ims_t, [f"t{name}" for name in gwfs])
    return sim


def layout(idx):
    return "single" if idx == 0 else "exchange"


def gwf_names(idx):
    return ["gwf"] if idx == 0 else ["gwfa", "gwfb"]


def read_heads(ws, idx):
    return np.concatenate(
        [
            flopy.utils.HeadFile(Path(ws) / f"{name}.hds").get_data().ravel()
            for name in gwf_names(idx)
        ]
    )


def build_models(idx, test):
    return build_sim(test.workspace, "mf6", layout(idx)), None


def check_output(idx, test):
    ws = Path(test.workspace)

    # head in A agrees with the independently computed solution
    h_a = read_heads(ws, idx)[0]
    h_expected = expected_head_a()
    print(f"head in A: {h_a:.6f}, expected {h_expected:.6f}, head in B {h_b}")
    assert h_a > h_b, "head in A should be above the head in B"
    assert np.isclose(h_a, h_expected, atol=1e-5), (
        f"head in A {h_a} differs from expected {h_expected}"
    )

    # budgets close
    for name in gwf_names(idx):
        disc = [
            abs(float(v))
            for v in re.findall(
                r"PERCENT DISCREPANCY =\s+(-?\d+\.\d+)",
                (ws / f"{name}.lst").read_text(),
            )
        ]
        assert disc and max(disc) < 0.01, f"{name} budget discrepancy {disc}"

    # DRHODC = 0 gives the same heads as a model without BUY
    heads = {}
    for tag, buy in (("zero", True), ("nobuy", False)):
        ws_cmp = ws / tag
        sim = build_sim(ws_cmp, test.targets["mf6"], layout(idx), rhodc=0.0, buy=buy)
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True)
        assert success, f"{tag} comparison run failed"
        heads[tag] = read_heads(ws_cmp, idx)
    assert np.array_equal(heads["zero"], heads["nobuy"]), (
        f"DRHODC = 0 heads {heads['zero']} differ from no-BUY {heads['nobuy']}"
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
