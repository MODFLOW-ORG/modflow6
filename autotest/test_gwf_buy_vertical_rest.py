"""
Test that the BUY buoyancy term gives no flow in a fluid at rest across a
vertical connection with a partially saturated convertible cell.

In the hydraulic head formulation the buoyancy flow between cells n and m is
proportional to the pressure difference between the two cell elevations
less the weight of the fluid between them, with the average density weighted
by a factor wt. For a fluid at rest, with each cell hydrostatic in its own
density, the flow is zero only if wt places the change in density at the
shared face. The weight from CL1 and CL2 does this for full cells. The
elevation of a partially saturated convertible cell is the center of its
saturated thickness, so the weight must be recalculated for that elevation.

The model is two cells in one column: a fresh or brackish top cell, full
(confined) or partially saturated (convertible), over a saline cell. Each
case is run with the default hydraulic head formulation and with
HHFORMULATION_RHS.

    test_vertical_rest         both cells are constant heads at their exact
                               hydrostatic values, so the face flow is the
                               discrete flow at those heads; it must be zero
    test_vertical_rest_solved  the saline cell is a constant head and the
                               partially saturated top cell is solved, with
                               and without the Newton-Raphson formulation;
                               its head must be the at-rest water table
"""

import flopy
import pytest

denseref = 62.428
drhodc = 0.0446
c_salt = 35.0
area = 1.0e4
top0, bot0, bot1 = 45.0, -1.3, -500.0
kv = 1.0
eta = 0.03  # at-rest water table in a partially saturated top cell


def rho(conc):
    return 1.0 + drhodc * conc / denseref


def rest_heads(conc_top, convertible):
    """Exact hydrostatic hydraulic heads, h = p / rho + z, of the two cells."""
    if convertible:
        h0 = eta  # zero pressure at the water table
        z0 = 0.5 * (h0 + bot0)
    else:
        h0 = 10.0
        z0 = 0.5 * (top0 + bot0)
    # pressure (in reference-density head) at the shared face
    p_face = rho(conc_top) * (h0 - z0) + rho(conc_top) * (z0 - bot0)
    h1 = p_face / rho(c_salt) + bot0
    return h0, h1


def run_model(ws, exe, conc_top, convertible, newton, hhform, solve_top):
    """Return the heads and the flow into the top cell."""
    h0, h1 = rest_heads(conc_top, convertible)
    conc = [conc_top, c_salt]
    diskw = dict(
        nodes=2,
        nja=4,
        top=[top0, bot0],
        bot=[bot0, bot1],
        area=area,
        iac=[2, 2],
        ja=[0, 1, 1, 0],
        ihc=[1, 0, 1, 0],
        cl12=[0.0, 0.5 * (top0 - bot0), 0.0, 0.5 * (bot0 - bot1)],
        hwva=[0.0, area, 0.0, area],
    )
    chd = [((1,), h1)] if solve_top else [((0,), h0), ((1,), h1)]
    strt = [2.0 if solve_top else h0, h1]

    sim = flopy.mf6.MFSimulation(sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim)
    ims_f = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_maximum=100,
        outer_dvclose=1e-10,
        inner_maximum=50,
        inner_dvclose=1e-12,
        linear_acceleration="BICGSTAB",
    )
    ims_t = flopy.mf6.ModflowIms(
        sim, filename="tran.ims", linear_acceleration="BICGSTAB"
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", newtonoptions="NEWTON" if newton else None
    )
    flopy.mf6.ModflowGwfdisu(gwf, **diskw)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=[1 if convertible else 0, 0], k=10.0, k33=kv, save_flows=True
    )
    flopy.mf6.ModflowGwfbuy(
        gwf,
        denseref=denseref,
        hhformulation_rhs=hhform,
        packagedata=[(0, drhodc, 0.0, "gwt", "CONCENTRATION")],
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    flopy.mf6.ModflowGwtdisu(gwt, **diskw)
    flopy.mf6.ModflowGwtic(gwt, strt=conc)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtcnc(
        gwt, stress_period_data=[((n,), c) for n, c in enumerate(conc)]
    )
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea="gwf", exgmnameb="gwt")
    sim.register_solution_package(ims_f, ["gwf"])
    sim.register_solution_package(ims_t, ["gwt"])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success, f"model in {ws} did not run"
    heads = flopy.utils.HeadFile(ws / "gwf.hds").get_data().ravel()
    flowja = (
        flopy.utils.CellBudgetFile(ws / "gwf.cbc")
        .get_data(text="FLOW-JA-FACE")[0]
        .ravel()
    )
    return heads, flowja[1]  # row of the top cell, connection to the bottom cell


@pytest.mark.parametrize("hhform", [False, True])
@pytest.mark.parametrize("newton", [False, True])
@pytest.mark.parametrize("convertible", [False, True])
@pytest.mark.parametrize("conc_top", [0.0, 17.5])
def test_vertical_rest(function_tmpdir, targets, conc_top, convertible, newton, hhform):
    _, q = run_model(
        function_tmpdir, targets["mf6"], conc_top, convertible, newton, hhform, False
    )
    cv = area / (0.5 * (top0 - bot0) / kv + 0.5 * (bot0 - bot1) / kv)
    print(f"flow into the top cell: {q:.6g} (conductance {cv:.6g})")
    # without the correction the partially saturated cases give 9 to 19
    # ft3/d; full cells give round-off of order 1e-7 ft3/d
    assert abs(q) < 1e-6 * cv, f"flow in a fluid at rest: {q}"


@pytest.mark.parametrize("hhform", [False, True])
@pytest.mark.parametrize("newton", [False, True])
@pytest.mark.parametrize("conc_top", [0.0, 17.5])
def test_vertical_rest_solved(function_tmpdir, targets, conc_top, newton, hhform):
    heads, _ = run_model(
        function_tmpdir, targets["mf6"], conc_top, True, newton, hhform, True
    )
    print(f"solved water table {heads[0]:.8f}, at rest {eta}")
    # without the correction the head is 0.25 to 0.5 ft off
    assert abs(heads[0] - eta) < 1e-6, f"solved head {heads[0]} is not {eta}"
