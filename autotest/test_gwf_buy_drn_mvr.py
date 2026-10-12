"""Independent routed-water check for BUY/DRN, including discharge depth.

Uniform initial concentration fixes density for the single GWF step. Two
zero-leak SFR reaches receive routed drain water. GWT/SFT/MVT run, but this
checker admits water routing only.
"""

import copy
import itertools

import flopy
import pytest
from framework import TestFramework


def factor(depth, newton):
    if depth == 0:
        return 1.0
    return 0.375 if newton else 0.5


def build(ws, cfg):
    depth, newton, multiplier, concentration, fraction = cfg
    f = factor(depth, newton)
    anchor = 1.5 + 0.5 * multiplier * f
    elevation = 2.0 if depth < 0 else 1.0
    sim = flopy.mf6.MFSimulation(sim_name="mvr", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    g = flopy.mf6.ModflowGwf(
        sim,
        modelname="flow",
        save_flows=True,
        newtonoptions="NEWTON" if newton else None,
    )
    imsg = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        outer_maximum=200,
        inner_maximum=100,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(imsg, ["flow"])
    flopy.mf6.ModflowGwfdis(g, nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfic(g, strt=anchor + 0.05)
    flopy.mf6.ModflowGwfnpf(g, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfchd(
        g,
        pname="CHD",
        auxiliary=["CONCENTRATION"],
        stress_period_data=[[(0, 0, 0), anchor, concentration]],
    )
    flopy.mf6.ModflowGwfdrn(
        g,
        pname="DRN",
        mover=True,
        auxiliary=["DEPTH", "MULT"],
        auxdepthname="DEPTH",
        auxmultname="MULT",
        stress_period_data=[[(0, 0, 1), elevation, 1.0, depth, multiplier]],
    )
    flopy.mf6.ModflowGwfbuy(
        g,
        density_filerecord="flow.density",
        packagedata=[(0, 1.3, 0.0, "trans", "CONCENTRATION")],
    )
    flopy.mf6.ModflowGwfsfr(
        g,
        pname="SFR",
        mover=True,
        nreaches=2,
        packagedata=[
            (0, (0, 0, 0), 1.0, 1.0, 0.001, 2.0, 1.0, 0.0, 0.04, 1, 1.0, 0),
            (1, (0, 0, 1), 1.0, 1.0, 0.001, 2.0, 1.0, 0.0, 0.04, 1, 1.0, 0),
        ],
        connectiondata=[(0, -1), (1, 0)],
        save_flows=True,
        budget_filerecord="sfr.cbc",
    )
    flopy.mf6.ModflowGwfmvr(
        g,
        maxmvr=1,
        maxpackages=2,
        packages=[("DRN",), ("SFR",)],
        perioddata=[("DRN", 0, "SFR", 0, "FACTOR", fraction)],
        budget_filerecord="mvr.cbc",
    )
    flopy.mf6.ModflowGwfoc(
        g,
        head_filerecord="flow.hds",
        budget_filerecord="flow.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    t = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imst = flopy.mf6.ModflowIms(
        sim, filename="trans.ims", linear_acceleration="BICGSTAB"
    )
    sim.register_ims_package(imst, ["trans"])
    flopy.mf6.ModflowGwtdis(t, nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwtic(t, strt=concentration)
    flopy.mf6.ModflowGwtmst(t, porosity=0.35)
    flopy.mf6.ModflowGwtadv(t, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(t, sources=[("CHD", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtsft(
        t,
        pname="SFT",
        flow_package_name="SFR",
        save_flows=True,
        budget_filerecord="sft.cbc",
        concentration_filerecord="sft.ucn",
        packagedata=[(0, concentration), (1, concentration)],
        reachperioddata=[(0, "STATUS", "ACTIVE"), (1, "STATUS", "ACTIVE")],
    )
    flopy.mf6.ModflowGwtmvt(t)
    flopy.mf6.ModflowGwtoc(
        t,
        budget_filerecord="trans.cbc",
        concentration_filerecord="trans.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea="flow", exgmnameb="trans")
    sim.write_simulation(silent=True)
    return sim


def read(ws):
    cb = flopy.utils.CellBudgetFile(ws / "flow.cbc", precision="double")
    sfr = flopy.utils.CellBudgetFile(ws / "sfr.cbc", precision="double")

    def rate(budget, name):
        return float(sum(x["q"].sum() for x in budget.get_data(text=name)))

    return {
        "heads": flopy.utils.HeadFile(ws / "flow.hds", precision="double")
        .get_data()
        .ravel()
        .tolist(),
        "qchd": rate(cb, "CHD"),
        "qdrn": rate(cb, "DRN"),
        "qtomvr": rate(cb, "DRN-TO-MVR"),
        "qfrommvr": rate(sfr, "FROM-MVR"),
        "qout": rate(sfr, "EXT-OUTFLOW"),
        "qsfrgwf": rate(sfr, "GWF"),
        "density": flopy.utils.HeadFile(
            ws / "flow.density", text="DENSITY", precision="double"
        )
        .get_data()
        .ravel()
        .tolist(),
    }


def errors(cfg, x):
    depth, newton, multiplier, concentration, fraction = cfg
    alpha = 1.0 + 1.3 * concentration / 1000.0
    f = factor(depth, newton)
    anchor = 1.5 + 0.5 * multiplier * f
    q = 0.5 * multiplier * f * alpha
    return {
        "head": max(abs(x["heads"][0] - anchor), abs(x["heads"][1] - 1.5)),
        "density": max(abs(v - 1000.0 * alpha) for v in x["density"]),
        "chd": abs(x["qchd"] - q),
        "drain_external": abs(x["qdrn"] + (1.0 - fraction) * q),
        "drain_mover": abs(x["qtomvr"] + fraction * q),
        "receiver_mover": abs(x["qfrommvr"] - fraction * q),
        "receiver_out": abs(x["qout"] + fraction * q),
        "receiver_leak": abs(x["qsfrgwf"]),
        "gwf_balance": abs(x["qchd"] + x["qdrn"] + x["qtomvr"]),
        "sfr_balance": abs(x["qfrommvr"] + x["qout"] + x["qsfrgwf"]),
    }


def negative(cfg, x):
    bad = copy.deepcopy(x)
    bad["heads"][1] += 0.01
    assert errors(cfg, bad)["head"] > 0.009
    # Corrupt routed water while preserving each model's water budget.
    bad = copy.deepcopy(x)
    bad["qtomvr"] -= 0.01
    bad["qdrn"] += 0.01
    bad["qfrommvr"] += 0.01
    bad["qout"] -= 0.01
    e = errors(cfg, bad)
    assert e["gwf_balance"] < 1e-8 and e["sfr_balance"] < 1e-8
    assert e["drain_mover"] > 0.009 and e["receiver_mover"] > 0.009


settings = list(
    itertools.product([-1.0, 1.0], [False, True], [1.0], [0.0, 16.0], [0.4, 1.0])
)
cases = [f"buydrnmvr{i}" for i in range(len(settings))]


def build_models(idx, test):
    return build(test.workspace, settings[idx]), None


def check_output(idx, test):
    x = read(test.workspace)
    e = errors(settings[idx], x)
    assert max(e.values()) < 1e-8, e
    negative(settings[idx], x)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        overwrite=False,
    )
    test.run()
