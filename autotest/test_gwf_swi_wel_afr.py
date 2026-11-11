"""
WEL AUTO_FLOW_REDUCE with the SWI Package. With SWI active, the well flow
reduction is applied over the part of the cell occupied by the fluid the model
simulates rather than over the whole cell. The reduction interval for a
freshwater model extends from max(cell bottom, zeta) to the cell top, and for a
saltwater model from the cell bottom to min(cell top, zeta). The reduced rate
uses min(head, top of the reduction interval), so in an unconfined freshwater
cell the controlling thickness is (hf - zeta).

The reduction interval is updated once per outer iteration from the current
heads: directly on the first outer iteration of a time step, and moved halfway
toward the current interface on later iterations. No Newton derivative is
included for the interface-controlled part of the reduction. Because the
interval is lagged and under-relaxed, it can remain offset from the converged
interface by up to one time step when the interface moves more than the
reduction thickness in a single step. The reduction also applies to confined
cells: in standard MODFLOW a confined cell cannot dewater, but with SWI its
fluid zone can vanish as the interface moves through it.

Cases:
1. swi-welafr01 (single-fluid, unconfined): two-cell drawdown; the simulated
   rate must match the analytic flow reduction on (hf - zeta) at the converged
   heads, and the same model without SWI must show no reduction.
2. swi-welafr02 (two-fluid, freshwater well): a five-layer vertical column with
   coupled freshwater and saltwater models; at every time step the simulated
   rate must lie between the analytic rates for the current and the previous
   reduction interval as the interface rises.
3. swi-welafr03 (single-fluid, CONFINED): a deep pressurized cell
   (icelltype=0, hf far above the cell top); the rate must be reduced on
   (cell top - zeta) as the interface rises into the cell, even though the head
   never approaches the bottom of the reduction interval.
4. swi-welafr04 (two-fluid, saltwater well): the saltwater model pumps the
   bottom cell; the rate must be reduced on (zeta - cell bottom) as the
   saltwater zone drains (the top of the reduction interval equals zeta once
   the interface enters the cell).
5. swi-welafr05 (input check): an out-of-range per-well
   AUTO_FLOW_REDUCE_AUXNAME value in a CONFINED cell must be rejected by the
   WEL input check when SWI has modified the reduction interval (without SWI,
   confined cells are skipped by the check because they cannot be reduced).
"""

import flopy
import pytest
from framework import TestFramework

cases = [
    "swi-welafr01",
    "swi-welafr02",
    "swi-welafr03",
    "swi-welafr04",
    "swi-welafr05",
]
XFAIL = {4}

alphaf = 1000.0 / (1025.0 - 1000.0)  # 40.0
alphas = 1025.0 / (1025.0 - 1000.0)  # 41.0
afr = 0.1  # AUTO_FLOW_REDUCE fraction of cell thickness


def sQSaturation(top, bot, x):
    """python replica of the mf6 quadratic smoothing function."""
    b = top - bot
    s = (x - bot) / b
    if s < 0.0:
        return 0.0
    if s > 1.0:
        return 1.0
    w = x - bot
    return -2.0 / b**3 * w**3 + 3.0 / b**2 * w**2


def expected_qmult(h, cell_top, cell_bot, ict, botred, topred, frac=afr):
    """expected WEL rate multiplier for the reduction interval [botred, topred]."""
    if not (ict != 0 or botred > cell_bot or topred < cell_top):
        return 1.0
    band = frac * (cell_top - cell_bot)
    if band <= 0.0:
        return 1.0
    return sQSaturation(botred + band, botred, min(h, topred))


def dbd_ims(sim, **kwargs):
    return flopy.mf6.ModflowIms(
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
        outer_dvclose=1.0e-7,
        inner_dvclose=1.0e-8,
        linear_acceleration="bicgstab",
        **kwargs,
    )


# ---------- cases 1 and 3: single-fluid, two cells (CHD supply + well) ----------

sf = {
    # unconfined: interface at the cell bottom at the supply head
    0: dict(top=1.0, bot=-40.0, icelltype=1, hchd=1.0, q=-1.0),
    # confined: deep pressurized cell, interface enters the cell from below
    2: dict(top=-30.0, bot=-40.0, icelltype=0, hchd=1.0, q=-0.1),
    # confined + invalid per-well auxiliary reduction value: must fail the WEL
    # input check because SWI has modified the reduction interval
    4: dict(top=-30.0, bot=-40.0, icelltype=0, hchd=1.0, q=-0.1, afraux=0.0),
}
sf_k = 0.01


def build_sim_sf(ws, exe, p, with_swi):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(400.0, 40, 1.0)])
    dbd_ims(sim)
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=2,
        delr=1.0,
        delc=1.0,
        top=p["top"],
        botm=p["bot"],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=p["hchd"])
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=p["icelltype"], k=sf_k)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=0.2, transient=True)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), p["hchd"]]])
    wel_kwargs = dict(auto_flow_reduce=afr)
    spd = [[(0, 0, 1), p["q"]]]
    if "afraux" in p:
        wel_kwargs["auxiliary"] = ["afr"]
        wel_kwargs["auto_flow_reduce_auxname"] = "afr"
        spd = [[(0, 0, 1), p["q"], p["afraux"]]]
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=spd, **wel_kwargs)
    if with_swi:
        flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budget_filerecord=name + ".bud",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def check_output_sf(idx, test):
    p = sf[idx]
    gwf = test.sims[0].gwf[0]
    hf = float(gwf.output.head().get_data()[0, 0, 1])
    q_sim = float(gwf.output.budget().get_data(text="WEL")[-1]["q"][0])
    zeta = -alphaf * hf  # single-fluid, hs = 0
    botred = max(p["bot"], zeta)
    q_exp = p["q"] * expected_qmult(
        hf, p["top"], p["bot"], p["icelltype"], botred, p["top"]
    )
    print(f"hf={hf:.6f} zeta={zeta:.4f} q_sim={q_sim:.6f} q_exp={q_exp:.6f}")

    # -- the interface must be inside the cell (reduction controlled by zeta)
    assert zeta > p["bot"] + 0.5, "test should move the interface into the cell"
    if p["icelltype"] == 0:
        # -- confined: the head must stay far above the cell top, proving the
        #    reduction acted on (top - zeta) because min(h, top of the
        #    reduction interval) is the cell top
        assert hf > p["top"] + 1.0, "cell should remain pressurized"

    # -- the rate must be meaningfully reduced and match the analytic rate
    assert abs(q_sim) < 0.9 * abs(p["q"]), (
        f"rate should be reduced by the interface, got {q_sim}"
    )
    assert abs(q_sim - q_exp) < 1.0e-6 + 1.0e-4 * abs(p["q"]), (
        f"simulated rate {q_sim} != analytic rate {q_exp}"
    )

    # -- same model without SWI: no reduction (unconfined case: the head stays
    #    above the reduction thickness over the cell bottom; confined case: the
    #    flow reduction is not applied at all)
    ws = test.workspace / "noswi"
    sim2 = build_sim_sf(ws, test.targets["mf6"], p, with_swi=False)
    sim2.write_simulation(silent=True)
    success, _ = sim2.run_simulation(silent=True)
    assert success, "no-SWI companion model failed"
    q2 = float(sim2.gwf[0].output.budget().get_data(text="WEL")[-1]["q"][0])
    print(f"no-SWI: q={q2:.6f}")
    assert abs(q2 - p["q"]) < 1.0e-6, (
        f"without SWI the reduction must stay inactive, got {q2}"
    )


# ---------- cases 2 and 4: two-fluid vertical column ----------

nlay = 5
# a cell top a little above the initial water table keeps the freshwater
# well's reduction controlled by its own head (min(h, top of the reduction
# interval) = h); the case where the cell top or zeta controls is exercised by
# the confined and saltwater cases instead
tf_top = 0.05
tf_botm = [-1.0, -2.0, -3.0, -4.0, -5.0]
tf_sy = 0.2
zeta0 = -0.5

tf = {
    # fresh well in the top cell, salt supplied by a bottom CHD
    1: dict(saltwell=False, q=-0.02, nstp=20),
    # salt well in the bottom cell, fresh supplied by a top CHD. The reduction
    # thickness is half the cell (0.5 m), wider than the 0.25 m the interface
    # moves per step, so the reduced rate is resolved over several steps. With
    # the 0.1 m reduction thickness of the other cases the interface would cross
    # it within one step, and the well would pump at the full rate until the
    # zone was gone without the reduction being simulated.
    3: dict(saltwell=True, q=-0.05, nstp=30, afr=0.5),
}


def build_gwf_tf(sim, is_salt, p):
    name = "saltwater" if is_salt else "freshwater"
    has_well = is_salt == p["saltwell"]
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=1,
        ncol=1,
        delr=1.0,
        delc=1.0,
        top=tf_top,
        botm=tf_botm,
    )
    strt = 0.0 if is_salt else -zeta0 / alphaf
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    # the salt model stays confined (icelltype=0): case 4 checks that the flow
    # reduction is active when SWI narrows the interval in a confined cell
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0 if is_salt else 1, k=100.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=tf_sy, transient=True)
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    well_cell = (nlay - 1, 0, 0) if is_salt else (0, 0, 0)
    chd_cell = (nlay - 1, 0, 0) if is_salt else (0, 0, 0)
    if has_well:
        flopy.mf6.ModflowGwfwel(
            gwf,
            auto_flow_reduce=p.get("afr", afr),
            stress_period_data=[[well_cell, p["q"]]],
        )
    else:
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[chd_cell, strt]])
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budget_filerecord=name + ".bud",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def build_sim_tf(test, p):
    sim = flopy.mf6.MFSimulation(
        sim_name="mymodel", sim_ws=test.workspace, exe_name=test.targets["mf6"]
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(p["nstp"], p["nstp"], 1.0)])
    ims = dbd_ims(sim)
    build_gwf_tf(sim, False, p)
    build_gwf_tf(sim, True, p)
    flopy.mf6.ModflowSwiswi(
        sim,
        exgtype="SWI6-SWI6",
        exgmnamea="freshwater",
        exgmnameb="saltwater",
    )
    sim.register_ims_package(ims, ["freshwater", "saltwater"])
    return sim


def check_output_tf(idx, test):
    p = tf[idx]
    ws = test.workspace
    hf_obj = flopy.utils.HeadFile(ws / "freshwater.hds")
    hs_obj = flopy.utils.HeadFile(ws / "saltwater.hds")
    wellmodel = "saltwater" if p["saltwell"] else "freshwater"
    gwf = next(m for m in test.sims[0].gwf if m.name == wellmodel)
    bud = gwf.output.budget()
    times = hf_obj.times

    if p["saltwell"]:
        klay, cell_top, cell_bot, ict = nlay - 1, tf_botm[-2], tf_botm[-1], 0
    else:
        klay, cell_top, cell_bot, ict = 0, tf_top, tf_botm[0], 1

    q_sim_last = None
    prev_bounds = None
    for t in times:
        hf = float(hf_obj.get_data(totim=t)[klay, 0, 0])
        hs = float(hs_obj.get_data(totim=t)[klay, 0, 0])
        zeta = -alphaf * hf + alphas * hs
        h_well = hs if p["saltwell"] else hf
        botred = cell_bot if p["saltwell"] else max(cell_bot, zeta)
        topred = min(cell_top, zeta) if p["saltwell"] else cell_top
        if prev_bounds is None:
            prev_bounds = (botred, topred)
        q_sim = float(bud.get_data(text="WEL", totim=t)[0]["q"][0])
        q_exp = p["q"] * expected_qmult(
            h_well, cell_top, cell_bot, ict, botred, topred, p.get("afr", afr)
        )
        # the reduction interval is lagged and under-relaxed across outer
        # iterations, so when the interface moves more than the reduction
        # thickness in one step the interval used by WEL can be offset from the
        # converged interface by up to one step: accept any rate between the
        # rates computed with the current and the previous reduction interval
        q_lag = p["q"] * expected_qmult(
            h_well,
            cell_top,
            cell_bot,
            ict,
            prev_bounds[0],
            prev_bounds[1],
            p.get("afr", afr),
        )
        tol = 1.0e-6 + 1.0e-3 * abs(p["q"])
        lo, hi = min(q_exp, q_lag) - tol, max(q_exp, q_lag) + tol
        print(
            f"t={t:5.1f} h={h_well:9.6f} zeta={zeta:9.6f} "
            f"q_sim={q_sim:9.6f} q_exp={q_exp:9.6f} q_lag={q_lag:9.6f}"
        )
        assert lo <= q_sim <= hi, (
            f"t={t}: simulated rate {q_sim} outside expected range "
            f"[{lo}, {hi}] (current {q_exp}, lagged {q_lag})"
        )
        prev_bounds = (botred, topred)
        q_sim_last = q_sim

    # -- the fluid zone is finite: the flow reduction must have reduced the well rate
    assert abs(q_sim_last) < 0.5 * abs(p["q"]), (
        f"the well rate should be substantially reduced by the interface, "
        f"got {q_sim_last} of requested {p['q']}"
    )


# ---------- driver ----------


def build_models(idx, test):
    if idx in sf:
        return build_sim_sf(test.workspace, test.targets["mf6"], sf[idx], True), None
    return build_sim_tf(test, tf[idx]), None


def check_output_ck(idx, test):
    with open(test.workspace / "mfsim.lst", "r") as f:
        lines = f.readlines()
    assert any("AUTO_FLOW_REDUCE_AUXNAME" in line for line in lines), (
        "expected an AUTO_FLOW_REDUCE_AUXNAME range error in mfsim.lst"
    )


def check_output(idx, test):
    if idx in XFAIL:
        check_output_ck(idx, test)
    elif idx in sf:
        check_output_sf(idx, test)
    else:
        check_output_tf(idx, test)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=idx in XFAIL,
    )
    test.run()
