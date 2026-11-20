"""TVK/TVS PERIOD values, covering: a TS-linked value continuing to track
its series when a later PERIOD block reappears only for a different,
unrelated cell (tvkcontinue, tvscontinue); a TS-linked value switching
cleanly between two different series across periods with no leftover
value from the prior one (tvkswitch); and a plain (non-TS) value
persisting across later, present-but-empty PERIOD blocks without being
applied again (tvkempty, tvsempty).
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvkcontinue", "tvkswitch", "tvscontinue", "tvkempty", "tvsempty"]

# tvkcontinue / tvscontinue: TS-linked cell (0,0,1) must keep tracking its
# series across periods 2/3, whose own PERIOD blocks reappear only for a
# different, unrelated cell (0,0,0), without repeating cell (0,0,1)'s
# own setting.
CONTINUE_TS_VALS = {
    "tvkcontinue": (10.0, 20.0, 30.0),
    "tvscontinue": (1.0e-4, 2.0e-4, 3.0e-4),
}

# tvkswitch: k_ts_a drives periods 1 and 3; k_ts_b drives period 2. Values
# are distinct at every stress period's start time so a stale link shows
# up as a value matching neither series.
TS_A_VALS = {0.0: 10.0, 1.0: 999.0, 2.0: 25.0, 3.0: 25.0}
TS_B_VALS = {0.0: 888.0, 1.0: 900.0, 2.0: 888.0, 3.0: 888.0}
SWITCH_EXPECTED = (TS_A_VALS[0.0], TS_B_VALS[1.0], TS_A_VALS[2.0])

# tvkempty: a plain (non-TS) K value set once in period 1 must persist
# across periods 2/3, whose own PERIOD blocks reappear present but empty.
# It is applied (logged) in period 1 only; heads must match a model with
# that K given in NPF.
EMPTY_K_VAL = 5.0

# tvsempty: as tvkempty, for a plain SS value, without storage change
# integration so the heads match a model with that SS given in STO.
SS0 = 1.0e-5
EMPTY_SS_VAL = 5.0e-4


def _build_base(name, ws):
    nlay, nrow, ncol = 1, 1, 3
    delr, delc = 1.0, 1.0
    top, botm = 0.0, [-1.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, print_option="NONE", linear_acceleration="BICGSTAB")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), -1.0], [(0, 0, ncol - 1), -2.0]],
        pname="CHD-1",
    )
    return sim, gwf


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim, gwf = _build_base(name, ws)
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=10.0, k33=10.0)

    if name == "tvkcontinue":
        tvk = flopy.mf6.ModflowUtltvk(
            npf,
            print_input=True,
            perioddata={
                0: [((0, 0, 1), "K", "k_ts")],
                # periods 1, 2: TVK's own PERIOD block reappears, setting a
                # different, unrelated cell's K, without repeating cell
                # (0,0,1)'s K setting
                1: [((0, 0, 0), "K", 5.0)],
                2: [((0, 0, 0), "K", 5.0)],
            },
            filename=f"{name}.tvk",
        )
        p1, p2, p3 = CONTINUE_TS_VALS[name]
        tvk.ts.initialize(
            filename="k.ts",
            timeseries=[(0.0, p1), (1.0, p2), (2.0, p3), (3.0, p3)],
            time_series_namerecord="k_ts",
            interpolation_methodrecord="stepwise",
        )
    elif name == "tvkswitch":
        tvk = flopy.mf6.ModflowUtltvk(
            npf,
            print_input=True,
            perioddata={
                0: [((0, 0, 1), "K", "k_ts_a")],
                1: [((0, 0, 1), "K", "k_ts_b")],
                2: [((0, 0, 1), "K", "k_ts_a")],
            },
            filename=f"{name}.tvk",
        )
        ts_names = ["k_ts_a", "k_ts_b"]
        ts_data = [(t, TS_A_VALS[t], TS_B_VALS[t]) for t in (0.0, 1.0, 2.0, 3.0)]
        tvk.ts.initialize(
            filename="k.ts",
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=["stepwise"] * len(ts_names),
        )
    elif name == "tvkempty":
        flopy.mf6.ModflowUtltvk(
            npf,
            print_input=True,
            perioddata={0: [((0, 0, 1), "K", EMPTY_K_VAL)], 1: [], 2: []},
            filename=f"{name}.tvk",
        )
    elif name == "tvsempty":
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=0,
            ss=SS0,
            sy=0.1,
            steady_state={0: False},
            transient={0: True},
        )
        flopy.mf6.ModflowUtltvs(
            sto,
            print_input=True,
            disable_storage_change_integration=True,
            perioddata={0: [((0, 0, 1), "SS", EMPTY_SS_VAL)], 1: [], 2: []},
            filename=f"{name}.tvs",
        )
    else:  # tvscontinue
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=0,
            ss=1e-5,
            sy=0.1,
            steady_state={0: False},
            transient={0: True},
        )
        tvs = flopy.mf6.ModflowUtltvs(
            sto,
            print_input=True,
            perioddata={
                0: [((0, 0, 1), "SS", "ss_ts")],
                # periods 1, 2: TVS's own PERIOD block reappears, setting a
                # different, unrelated cell's SS, without repeating cell
                # (0,0,1)'s SS setting
                1: [((0, 0, 0), "SS", 5.0e-4)],
                2: [((0, 0, 0), "SS", 5.0e-4)],
            },
            filename=f"{name}.tvs",
        )
        p1, p2, p3 = CONTINUE_TS_VALS[name]
        tvs.ts.initialize(
            filename="ss.ts",
            timeseries=[(0.0, p1), (1.0, p2), (2.0, p3), (3.0, p3)],
            time_series_namerecord="ss_ts",
            interpolation_methodrecord="stepwise",
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def check_output(idx, test):
    name = test.name
    lst_fname = test.workspace / f"{name}.lst"
    text = lst_fname.read_text()

    field = "SS" if name.startswith("tvs") else "K"
    pattern = re.compile(
        rf"Setting {field} value for cell\s+(\S+)\s+in stress\s+period\s+(\d+),"
        r"\s+time step\s+\d+\s*=\s*([\d.eE+-]+)"
    )
    matches = pattern.findall(text)
    assert len(matches) > 0, (
        f"No 'Setting {field} value' lines found in list file {lst_fname}"
    )

    # only track the TS-linked cell (0,0,1) -> 1-based cellstr "(1,1,2)";
    # ignore the distractor cell (0,0,0) used in the continue cases to
    # force the PERIOD block to reappear in periods 2/3
    vals_by_period = {}
    for cellstr, kper_str, val_str in matches:
        if "1,1,2" not in cellstr:
            continue
        kper = int(kper_str)
        val = float(val_str)
        vals_by_period.setdefault(kper, []).append(val)

    print(f"{name} applied {field} values by stress period: {vals_by_period}")

    if name in CONTINUE_TS_VALS:
        p1, p2, p3 = CONTINUE_TS_VALS[name]
        assert 1 in vals_by_period, f"No {field} applications logged for period 1"
        assert np.allclose(vals_by_period[1], p1), (
            f"Period 1 {field} expected {p1}, got {vals_by_period[1]}"
        )
        assert 2 in vals_by_period, (
            f"No {field} applications logged for period 2 -- TS-linked "
            f"{field} did not continue tracking after period 1 (period 2's "
            f"PERIOD block reappears for a different cell without "
            f"repeating this one). vals_by_period={vals_by_period}"
        )
        assert np.allclose(vals_by_period[2], p2), (
            f"Period 2 {field} expected {p2}, got {vals_by_period[2]}"
        )
        assert 3 in vals_by_period, f"No {field} applications logged for period 3"
        assert np.allclose(vals_by_period[3], p3), (
            f"Period 3 {field} expected {p3}, got {vals_by_period[3]}"
        )
    elif name in ("tvkempty", "tvsempty"):
        val = EMPTY_K_VAL if name == "tvkempty" else EMPTY_SS_VAL
        assert list(vals_by_period) == [1], (
            f"{field} expected to be applied in period 1 only, got {vals_by_period}"
        )
        assert np.allclose(vals_by_period[1], val), (
            f"Period 1 {field} expected {val}, got {vals_by_period[1]}"
        )

        # reference: the same value given in NPF or STO, without TVK/TVS
        ws = test.workspace / "ref"
        sim, gwf = _build_base(name, ws)
        sim.exe_name = test.targets["mf6"]
        if name == "tvkempty":
            flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=[10.0, val, 10.0], k33=10.0)
        else:
            flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0, k33=10.0)
            flopy.mf6.ModflowGwfsto(
                gwf,
                iconvert=0,
                ss=[SS0, val, SS0],
                sy=0.1,
                steady_state={0: False},
                transient={0: True},
            )
        flopy.mf6.ModflowGwfoc(
            gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
        )
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True)
        assert success, "reference simulation failed"
        h = flopy.utils.HeadFile(test.workspace / f"{name}.hds").get_alldata()
        href = flopy.utils.HeadFile(ws / f"{name}.hds").get_alldata()
        assert np.allclose(h, href, atol=1e-10), (
            f"{field} did not persist across the empty PERIOD blocks: max head "
            f"difference from the reference {np.abs(h - href).max()}"
        )
    else:  # tvkswitch
        p1, p2, p3 = SWITCH_EXPECTED
        assert 1 in vals_by_period, "No K applications logged for period 1"
        assert np.allclose(vals_by_period[1], p1), (
            f"Period 1 (k_ts_a) expected {p1}, got {vals_by_period[1]}"
        )
        assert 2 in vals_by_period, "No K applications logged for period 2"
        assert np.allclose(vals_by_period[2], p2), (
            f"Period 2 (switched to k_ts_b) expected {p2} -- if this "
            f"instead shows a k_ts_a value ({TS_A_VALS[1.0]}), a stale "
            f"link from period 1 is still driving the cell, "
            f"got {vals_by_period[2]}"
        )
        assert 3 in vals_by_period, "No K applications logged for period 3"
        assert np.allclose(vals_by_period[3], p3), (
            f"Period 3 (switched back to k_ts_a) expected {p3} -- if "
            f"this instead shows a k_ts_b value ({TS_B_VALS[2.0]}), a "
            f"stale link from period 2 is still driving the cell, "
            f"got {vals_by_period[3]}"
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
