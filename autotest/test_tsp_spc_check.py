"""
Test that SPC input is checked against the flow package in every time step, not
just in periods where the SPC file has a period block.

List-based SPC's MAXBOUND must cover the flow package's bound count, or a flow
package that outgrows MAXBOUND reads past the end of SPC's array. Flows come from
a GWF-GWT/GWE exchange or an FMI budget file. These differ because a budget file's
bound counts aren't current until FMI advances, so a check made when reading the
period is too early, even if the SPC file has a block.

Array-based SPC (SPCA) can't overrun since its size is the number of cells per
layer, but the same check requires that it's paired with array-based recharge or
evapotranspiration, and that list-based SPC isn't.
"""

import flopy
import pytest

ncol = 5
qwell = 1.0
tdis_rc = [(1.0, 1, 1.0)] * 3
# the well list grows in period 2
wel_spd = {
    0: [[(0, 0, 0), qwell]],
    1: [[(0, 0, j), qwell] for j in range(3)],
}
# the error messages, which may wrap in the mf6 output
maxbound_err = "has MAXBOUND set less"
needs_array_err = "Array-based recharge must be used"
needs_list_err = "List-based packages must be used"


def value(bndno, val):
    # the variable is named for concentration in flopy for either model type
    return (bndno, "concentration", val)


well_cases = {
    # MAXBOUND covers the wells
    "enough-block": (
        3,
        {0: [value(0, 100.0)], 1: [value(0, 100.0), value(1, 200.0), value(2, 300.0)]},
    ),
    # MAXBOUND covers the wells, with no block in the period the list grows
    "enough-no-block": (3, {0: [value(0, 100.0)]}),
    # too small, with no block in the period the list grows
    "small-no-block": (1, {0: [value(0, 100.0)]}),
    # too small, with a block in that period
    "small-block": (1, {0: [value(0, 100.0)], 1: [value(0, 100.0)]}),
}


def add_sim(ws, name, exe):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=len(tdis_rc), perioddata=tdis_rc)
    return sim


def add_gwf(sim, flow_pkg):
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowIms(sim, pname="ims-gwf", filename="gwf.ims")
    sim.register_ims_package(sim.get_package("ims-gwf"), ["gwf"])
    flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=ncol)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, ncol - 1), 1.0]])
    if flow_pkg == "wel":
        flopy.mf6.ModflowGwfwel(
            gwf, maxbound=3, pname="flow", stress_period_data=wel_spd
        )
    else:
        flopy.mf6.ModflowGwfrcha(gwf, pname="flow", recharge=0.1)
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="gwf.bud",
        head_filerecord="gwf.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )


def add_tsp(sim, mtype, add_spc):
    """Add a GWT or GWE model, with add_spc adding the flow package's SPC input."""
    if mtype == "gwt":
        tsp = flopy.mf6.ModflowGwt(sim, modelname="tsp")
        pkgs = (
            flopy.mf6.ModflowGwtdis,
            flopy.mf6.ModflowGwtic,
            flopy.mf6.ModflowGwtadv,
            flopy.mf6.ModflowGwtssm,
        )
    else:
        tsp = flopy.mf6.ModflowGwe(sim, modelname="tsp")
        pkgs = (
            flopy.mf6.ModflowGwedis,
            flopy.mf6.ModflowGweic,
            flopy.mf6.ModflowGweadv,
            flopy.mf6.ModflowGwessm,
        )
    dis, ic, adv, ssm = pkgs
    flopy.mf6.ModflowIms(
        sim, pname="ims-tsp", filename="tsp.ims", linear_acceleration="BICGSTAB"
    )
    sim.register_ims_package(sim.get_package("ims-tsp"), ["tsp"])
    dis(tsp, nrow=1, ncol=ncol)
    ic(tsp, strt=0.0)
    adv(tsp)
    if mtype == "gwt":
        flopy.mf6.ModflowGwtmst(tsp, porosity=0.3)
    else:
        flopy.mf6.ModflowGweest(
            tsp,
            porosity=0.3,
            heat_capacity_water=4180.0,
            density_water=1000.0,
            heat_capacity_solid=880.0,
            density_solid=2650.0,
        )
    ssm(tsp, sources=[()], fileinput=[("flow", "flow.spc")])
    add_spc(tsp)
    return tsp


def build(ws, exe, mtype, flows, flow_pkg, add_spc):
    if flows == "exchange":
        sim = add_sim(ws, "spc", exe)
        add_gwf(sim, flow_pkg)
        add_tsp(sim, mtype, add_spc)
        exg = flopy.mf6.ModflowGwfgwt if mtype == "gwt" else flopy.mf6.ModflowGwfgwe
        exg(
            sim,
            exgtype=f"GWF6-{mtype.upper()}6",
            exgmnamea="gwf",
            exgmnameb="tsp",
            filename="flow.exg",
        )
        return sim

    flow = add_sim(ws / "flow", "flow", exe)
    add_gwf(flow, flow_pkg)
    flow.write_simulation(silent=True)
    success, buff = flow.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff)

    sim = add_sim(ws, "spc", exe)
    tsp = add_tsp(sim, mtype, add_spc)
    fmi = flopy.mf6.ModflowGwtfmi if mtype == "gwt" else flopy.mf6.ModflowGwefmi
    fmi(
        tsp,
        packagedata=[
            ("GWFHEAD", "flow/gwf.hds", None),
            ("GWFBUDGET", "flow/gwf.bud", None),
        ],
    )
    return sim


def list_spc(maxbound, spd):
    def add_spc(tsp):
        flopy.mf6.ModflowUtlspc(
            tsp, maxbound=maxbound, perioddata=spd, filename="flow.spc", pname="spc"
        )

    return add_spc


def array_spc(spd):
    def add_spc(tsp):
        flopy.mf6.ModflowUtlspca(
            tsp, concentration=spd, filename="flow.spc", pname="spc"
        )

    return add_spc


def run(sim, errmsg=None):
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    if errmsg is None:
        assert success, "\n".join(buff)
    else:
        assert not success
        assert errmsg in " ".join(" ".join(buff).split())


@pytest.mark.parametrize("mtype", ["gwt", "gwe"])
@pytest.mark.parametrize("flows", ["exchange", "budget"])
@pytest.mark.parametrize("case", well_cases)
def test_maxbound(function_tmpdir, targets, mtype, flows, case):
    maxbound, spd = well_cases[case]
    sim = build(
        function_tmpdir, targets["mf6"], mtype, flows, "wel", list_spc(maxbound, spd)
    )
    run(sim, maxbound_err if case.startswith("small") else None)


@pytest.mark.parametrize("flows", ["exchange", "budget"])
def test_array_input_valid(function_tmpdir, targets, flows):
    """SPCA with array-based recharge passes the check in every period."""
    sim = build(
        function_tmpdir, targets["mf6"], "gwt", flows, "rcha", array_spc({0: 100.0})
    )
    run(sim)


@pytest.mark.parametrize("flows", ["exchange", "budget"])
def test_array_recharge_needs_array_input(function_tmpdir, targets, flows):
    """List-based SPC with array-based recharge, with no period block at all."""
    sim = build(
        function_tmpdir, targets["mf6"], "gwt", flows, "rcha", list_spc(ncol, {})
    )
    run(sim, needs_array_err)


@pytest.mark.parametrize("flows", ["exchange", "budget"])
def test_list_package_needs_list_input(function_tmpdir, targets, flows):
    """SPCA with a list-based flow package, with no period block at all."""
    sim = build(function_tmpdir, targets["mf6"], "gwt", flows, "wel", array_spc({}))
    run(sim, needs_list_err)
