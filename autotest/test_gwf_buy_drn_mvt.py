"""BUY must preserve complete DRN water and solute routing through MVR/MVT."""

import copy
import itertools
from pathlib import Path

import flopy
import pytest
from framework import TestFramework


def discharge(h, cfg):
    bottom = 1.0
    u = (
        max(0.0, min(1.0, (h - bottom) / abs(cfg["depth"])))
        if cfg["depth"]
        else float(h > bottom)
    )
    factor = u * u * (2 - u) if cfg["newton"] and cfg["depth"] else u
    return (
        cfg["multiplier"]
        * (1 + 1.3 * cfg["concentration"] / 1000)
        * factor
        * max(h - bottom, 0.0)
    )


def build_case(ws, cfg):
    c, depth, mult = (cfg[k] for k in ("concentration", "depth", "multiplier"))
    sim = flopy.mf6.MFSimulation(sim_name="moveraudit", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    g = flopy.mf6.ModflowGwf(
        sim,
        modelname="flow",
        save_flows=True,
        newtonoptions="NEWTON" if cfg["newton"] else None,
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        outer_maximum=200,
        inner_maximum=100,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(ims, ["flow"])
    grid = {
        "nlay": 1,
        "nrow": 1,
        "ncol": 2,
        "top": 1.0,
        "botm": 0.0,
        "delr": 1.0,
        "delc": 1.0,
    }
    flopy.mf6.ModflowGwfdis(g, **grid)
    flopy.mf6.ModflowGwfic(g, strt=cfg.get("initial_head", 1.4))
    flopy.mf6.ModflowGwfnpf(g, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfchd(
        g,
        pname="CHD",
        auxiliary=["CONCENTRATION"],
        stress_period_data=[[(0, 0, 0), cfg.get("anchor", 1.5), c]],
    )
    flopy.mf6.ModflowGwfdrn(
        g,
        pname="DRN",
        mover=True,
        auxiliary=["CONCENTRATION", "DEPTH", "MULT"],
        auxdepthname="DEPTH",
        auxmultname="MULT",
        stress_period_data=[
            [(0, 0, 1), 1.0 if depth >= 0 else 1.0 - depth, 1.0, 0.0, depth, mult]
        ],
    )
    flopy.mf6.ModflowGwfbuy(
        g,
        packagedata=[(0, 1.3, 0.0, "trans", "CONCENTRATION")],
        density_filerecord="flow.density",
    )
    flopy.mf6.ModflowGwfsfr(
        g,
        pname="SFR",
        mover=True,
        nreaches=2,
        packagedata=[
            (0, (0, 0, 0), 1.0, 1.0, 0.01, 2.1, 0.5, 0.0, 0.03, 1, 1.0, 0),
            (1, (0, 0, 1), 1.0, 1.0, 0.01, 2.0, 0.5, 0.0, 0.03, 1, 1.0, 0),
        ],
        connectiondata=[(0, -1), (1, 0)],
        budget_filerecord="sfr.cbc",
        save_flows=True,
    )
    flopy.mf6.ModflowGwfmvr(
        g,
        maxmvr=1,
        maxpackages=2,
        packages=[("DRN",), ("SFR",)],
        perioddata=[("DRN", 0, "SFR", 0, "FACTOR", cfg["fraction"])],
        budget_filerecord="mover.cbc",
        print_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        g,
        head_filerecord="flow.hds",
        budget_filerecord="flow.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    t = flopy.mf6.ModflowGwt(sim, modelname="trans")
    it = flopy.mf6.ModflowIms(
        sim,
        filename="trans.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(it, ["trans"])
    flopy.mf6.ModflowGwtdis(t, **grid)
    flopy.mf6.ModflowGwtic(t, strt=c)
    flopy.mf6.ModflowGwtmst(t, porosity=0.35)
    flopy.mf6.ModflowGwtadv(t, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(
        t, sources=[("CHD", "AUX", "CONCENTRATION"), ("DRN", "AUX", "CONCENTRATION")]
    )
    flopy.mf6.ModflowGwtsft(
        t,
        pname="SFR",
        flow_package_name="SFR",
        packagedata=[(0, c), (1, c)],
        save_flows=True,
        concentration_filerecord="sft.ucn",
        budget_filerecord="sft.cbc",
    )
    flopy.mf6.ModflowGwtmvt(t, print_flows=True, budget_filerecord="mvt.cbc")
    flopy.mf6.ModflowGwtoc(
        t,
        concentration_filerecord="trans.ucn",
        budget_filerecord="trans.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgmnamea="flow", exgmnameb="trans", filename="test.gwfgwt"
    )
    sim.write_simulation()
    return sim


def read(ws):
    out = {}
    for key, filename, text in [
        ("heads", "flow.hds", "HEAD"),
        ("concentration", "trans.ucn", "CONCENTRATION"),
        ("density", "flow.density", "DENSITY"),
        ("stream_concentration", "sft.ucn", "CONCENTRATION"),
    ]:
        out[key] = (
            flopy.utils.HeadFile(ws / filename, text=text, precision="double")
            .get_data()
            .ravel()
            .tolist()
        )
    for key, filename, text in [
        ("qdrn", "flow.cbc", "DRN"),
        ("qdrn_mvr", "flow.cbc", "DRN-TO-MVR"),
        ("qchd", "flow.cbc", "CHD"),
        ("qmover", "mover.cbc", "MOVER-FLOW"),
        ("qstream_in", "sfr.cbc", "FROM-MVR"),
        ("qstream_out", "sfr.cbc", "EXT-OUTFLOW"),
        ("moved_mass", "mvt.cbc", "MVT-FLOW"),
        ("stream_mass_in", "sft.cbc", "FROM-MVR"),
        ("stream_mass_out", "sft.cbc", "EXT-OUTFLOW"),
        ("stream_storage", "sft.cbc", "STORAGE"),
    ]:
        cb = flopy.utils.CellBudgetFile(ws / filename, precision="double")
        out[key] = float(sum(a["q"].sum() for a in cb.get_data(text=text)))
    cb = flopy.utils.CellBudgetFile(ws / "flow.cbc", precision="double")
    out["qface"] = float(cb.get_data(text="FLOW-JA-FACE")[0].ravel()[3])
    return out


def errors(cfg, x):
    ratio = 1 + 1.3 * cfg["concentration"] / 1000
    h = x["heads"][1]
    q = discharge(h, cfg)
    c0, c1 = x["concentration"]
    qf = ratio * (cfg.get("anchor", 1.5) - h)
    fraction = cfg["fraction"]
    lo, hi = 0.5, max(2.0, cfg.get("anchor", 1.5) + 1.0)
    for _ in range(65):
        mid = (lo + hi) / 2
        if ratio * (cfg.get("anchor", 1.5) - mid) > discharge(mid, cfg):
            lo = mid
        else:
            hi = mid
    return {
        "head": abs(h - (lo + hi) / 2),
        "density": max(abs(v - 1000 * ratio) for v in x["density"]),
        "total_drain": abs(x["qdrn"] + x["qdrn_mvr"] + q),
        "diverted": abs(x["qdrn_mvr"] + fraction * q),
        "external_drain": abs(x["qdrn"] + (1 - fraction) * q),
        "mover": abs(x["qmover"] + fraction * q),
        "stream_in": abs(x["qstream_in"] - fraction * q),
        "stream_out": abs(x["qstream_out"] + fraction * q),
        "physical_cell_water": abs(qf - q),
        "internal_pressure": abs(x["qface"] - qf),
        "reported_water": abs(x["qchd"] + x["qdrn"] + x["qdrn_mvr"]),
        "moved_mass": abs(x["moved_mass"] - fraction * q * c1),
        "received_mass": abs(x["stream_mass_in"] - fraction * q * c1),
        "stream_solute": abs(
            x["stream_mass_in"] + x["stream_mass_out"] + x["stream_storage"]
        ),
        "cell0_solute": abs(
            0.35 * (c0 - cfg["concentration"])
            - x["qchd"] * cfg["concentration"]
            + x["qface"] * c0
        ),
        "cell1_solute": abs(
            0.35 * (c1 - cfg["concentration"]) - x["qface"] * c0 + q * c1
        ),
    }


def negative(cfg, x):
    bad = copy.deepcopy(x)
    bad["heads"][1] += 0.01
    assert errors(cfg, bad)["head"] > 0.009
    # This corruption leaves total drain and groundwater water closure unchanged.
    bad = copy.deepcopy(x)
    bad["qdrn"] += 0.01
    bad["qdrn_mvr"] -= 0.01
    e = errors(cfg, bad)
    assert e["reported_water"] < 1e-8 and e["diverted"] > 0.009
    bad = copy.deepcopy(x)
    bad["moved_mass"] += 0.01
    assert errors(cfg, bad)["moved_mass"] > 0.009
    return True


configs = [
    {"depth": d, "newton": n, "concentration": c, "multiplier": 1.0, "fraction": 1.0}
    for d, n, c in itertools.product([-1.0, 0.0, 1.0], [False, True], [0.0, 16.0])
]
configs += [
    {
        "depth": d,
        "newton": True,
        "concentration": 16.0,
        "multiplier": m,
        "fraction": 0.4,
    }
    for d, m in itertools.product([-1.0, 1.0], [0.0, 2.0])
]
cases = [f"buy_drn_mvr_{i:02}" for i in range(len(configs))]


def build_models(idx, test):
    return build_case(test.workspace, configs[idx]), None


def check_output(idx, test):
    cfg = configs[idx]
    x = read(Path(test.workspace))
    residuals = errors(cfg, x)
    assert max(residuals.values()) < 1e-8, residuals
    assert negative(cfg, x)


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
