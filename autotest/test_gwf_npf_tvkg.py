"""
Test full-grid mode grid array input with TVKG (time-varying hydraulic
conductivity given as grid arrays). Each case is compared with TVK list
input giving the equivalent static changes in each period.

Grid: 2 layers, 1 row, 5 columns; cell (2, 1, 5) is inactive. CHD at the
two ends of layer 1. One time step per period, period length 1.

Case 1 (tvkg01): static K and K33 arrays in period 2, K in period 3
    (incl. a value at the inactive cell, which is ignored).
Case 2 (tvkg02): K series (stepwise) in period 1, no block in period 2,
    K33 only in period 3; the K link keeps progressing.
Case 3 (tvkg03): K series with DNODATA at some times; a cell keeps its
    current value (NPF value if never set) where the series is DNODATA.
Case 4 (tvkg04): K series in period 1, empty block in period 2 (link
    kept), static K in period 3 (link replaced), no block in period 4.
Case 5 (tvkg05): K and K33 with the NPF K33OVERK option.
Case 6 (tvkg06): TVK6 and TVKG6 both given. Input error.
Case 7 (tvkg07): K series (linear) at one cell, DNODATA elsewhere at
    every time; period averages are applied.
Case 8 (tvkg08): K series (linear) defined at one cell at time 0 only
    (DNODATA at time 1). Input error.
Case 9 (tvkg09): K series and a static K33 array in the period 1 block,
    no block in period 2 (the K link continues), static K in period 3
    (replaces the link).
Case 10 (tvkg10): K series (stepwise) in period 1 with a value at time 0
    and DNODATA at time 0.5; the period 1 step uses the DNODATA value, so
    the cell keeps its NPF value until the series changes it in period 2.
Case 11 (tvkg11): K33 given as a ratio of K (NPF K33OVERK) in period 2, K
    only in period 3; K33 keeps the value set in period 2. The reference
    gives K33 as its value (the ratio times K = 2), without K33OVERK.
Case 12 (tvkg12): K33 ratio at cell (1, 1, 2) in period 1; in period 2,
    K at that cell and K33 from a series that is DNODATA there and
    defined at cell (1, 1, 4). The cell keeps its period 1 K33. The
    reference gives K33 as its value, without K33OVERK.
Case 13 (tvkg13): in period 2, K at cell (1, 1, 2) from a series and a K33
    ratio array. Arrays given in a PERIOD block are applied before series
    values, so the ratio uses the K in effect before period 2 (2), not the
    series value of period 1 (7) or 2 (3). The reference gives K33 as its
    value, without K33OVERK.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_arraygrid_period02 import insert_empty_period

cases = [
    "tvkg01",
    "tvkg02",
    "tvkg03",
    "tvkg04",
    "tvkg05",
    "tvkg06",
    "tvkg07",
    "tvkg08",
    "tvkg09",
    "tvkg10",
    "tvkg11",
    "tvkg12",
    "tvkg13",
]
errors = {
    "tvkg06": "TVK6 and TVKG6 cannot both be specified.",
    "tvkg08": 'Time-array series "KTAS" with LINEAR interpolation is defined '
    "at only one of times 0.0000 and 1.0000 at user node 3.",
}
xfail = set(errors)

nlay, nrow, ncol = 2, 1, 5
inactive = (1, 0, 4)
DNODATA = 3.0e30
nper = {"tvkg04": 4}

# reference TVK list changes for each case: {period: [(cellid, name, value)]}
reference = {
    "tvkg01": {
        1: [((0, 0, 2), "K", 5.0), ((1, 0, 1), "K33", 0.1)],
        2: [((0, 0, 3), "K", 0.2)],
    },
    "tvkg02": {
        0: [((0, 0, 2), "K", 2.0)],
        1: [((0, 0, 2), "K", 5.0)],
        2: [((0, 0, 2), "K", 0.5), ((1, 0, 1), "K33", 0.1)],
    },
    "tvkg03": {
        0: [((0, 0, 2), "K", 2.0)],
        1: [((0, 0, 3), "K", 4.0)],
        2: [((0, 0, 2), "K", 0.5)],
    },
    "tvkg04": {
        0: [((0, 0, 2), "K", 2.0)],
        1: [((0, 0, 2), "K", 5.0)],
        2: [((0, 0, 2), "K", 0.3)],
    },
    "tvkg05": {1: [((0, 0, 1), "K", 3.0), ((0, 0, 1), "K33", 0.2)]},
    "tvkg06": {1: [((0, 0, 2), "K", 5.0)]},
    # linear 2 -> 4 over period 1, then 4: period averages 3, 4, 4
    "tvkg07": {
        0: [((0, 0, 2), "K", 3.0)],
        1: [((0, 0, 2), "K", 4.0)],
        2: [((0, 0, 2), "K", 4.0)],
    },
    "tvkg08": {},
    "tvkg09": {
        0: [((0, 0, 2), "K", 2.0), ((1, 0, 1), "K33", 0.1)],
        1: [((0, 0, 2), "K", 5.0)],
        2: [((0, 0, 2), "K", 0.3)],
    },
    "tvkg10": {1: [((0, 0, 2), "K", 5.0)]},
    # ratio cases (base K 2): K33 values are the ratio times K
    "tvkg11": {1: [((0, 0, 1), "K33", 0.4)], 2: [((0, 0, 1), "K", 3.0)]},
    "tvkg12": {
        0: [((0, 0, 1), "K33", 0.4)],
        1: [((0, 0, 1), "K", 3.0), ((0, 0, 3), "K33", 0.8)],
    },
    "tvkg13": {1: [((0, 0, 1), "K", 3.0), ((0, 0, 1), "K33", 0.4)]},
}
ratio_cases = ("tvkg11", "tvkg12", "tvkg13")


def grid(values):
    """Full grid array, DNODATA except at the given cells."""
    a = np.full((nlay, nrow, ncol), DNODATA)
    for cell, v in values.items():
        a[cell] = v
    return a


def series(records):
    """TAS records {time: {cell: value}} as full grid arrays."""
    return {t: grid(v) for t, v in records.items()}


def build_sim(ws, name, exe, mode):
    np_ = nper.get(name, 3)
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=np_, perioddata=[(1.0, 1, 1.0)] * np_)
    flopy.mf6.ModflowIms(sim, inner_dvclose=1e-10, outer_dvclose=1e-10)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[inactive] = 0
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=10.0,
        botm=[5.0, 0.0],
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    if name == "tvkg05":
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=1.0, k33=0.5, k33overk=True)
    elif name in ratio_cases and mode == "tvkg":
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=2.0, k33=0.5, k33overk=True)
    elif name in ratio_cases:
        # reference: K33 given as its value (the ratio times K = 2)
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=2.0, k33=1.0)
    else:
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=1.0, k33=1.0)
    if mode == "tvk" or name == "tvkg06":
        npf.tvk.initialize(filename=f"{name}.tvk", perioddata=reference[name])
    if mode == "tvkg":
        add_tvkg(npf, name)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 10.0], [(0, 0, 4), 5.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation(silent=True)
    if mode == "tvkg" and name == "tvkg04":
        fpth = ws / f"{name}.tvkg"
        fpth.write_text(insert_empty_period(fpth.read_text(), 2))
    return sim


def add_tvkg(npf, name):
    tas = None
    method = "stepwise"
    if name == "tvkg01":
        k = {1: grid({(0, 0, 2): 5.0}), 2: grid({(0, 0, 3): 0.2, inactive: 9.0})}
        k33 = {1: grid({(1, 0, 1): 0.1})}
    elif name == "tvkg02":
        k = {0: "TIMEARRAYSERIES ktas"}
        k33 = {2: grid({(1, 0, 1): 0.1})}
        tas = series(
            {
                t: {(0, 0, 2): v, inactive: 99.0}
                for t, v in zip([0.0, 1.0, 2.0, 3.0], [2.0, 5.0, 0.5, 0.5])
            }
        )
    elif name == "tvkg03":
        k = {0: "TIMEARRAYSERIES ktas"}
        k33 = None
        tas = series(
            {
                0.0: {(0, 0, 2): 2.0},
                1.0: {(0, 0, 3): 4.0},
                2.0: {(0, 0, 2): 0.5},
                3.0: {(0, 0, 2): 0.5},
            }
        )
    elif name == "tvkg04":
        k = {0: "TIMEARRAYSERIES ktas", 2: grid({(0, 0, 2): 0.3})}
        k33 = None
        tas = series(
            {
                t: {(0, 0, 2): v}
                for t, v in zip([0.0, 1.0, 2.0, 3.0, 4.0], [2.0, 5.0, 0.5, 0.5, 0.5])
            }
        )
    elif name == "tvkg05":
        k = {1: grid({(0, 0, 1): 3.0})}
        k33 = {1: grid({(0, 0, 1): 0.2})}
    elif name in ("tvkg07", "tvkg08"):
        k = {0: "TIMEARRAYSERIES ktas"}
        k33 = None
        method = "linear"
        v1 = 4.0 if name == "tvkg07" else DNODATA
        tas = series(
            {
                0.0: {(0, 0, 2): 2.0},
                1.0: {(0, 0, 2): v1},
                2.0: {(0, 0, 2): 4.0},
                3.0: {(0, 0, 2): 4.0},
            }
        )
    elif name == "tvkg09":
        k = {0: "TIMEARRAYSERIES ktas", 2: grid({(0, 0, 2): 0.3})}
        k33 = {0: grid({(1, 0, 1): 0.1})}
        tas = series(
            {
                t: {(0, 0, 2): v}
                for t, v in zip([0.0, 1.0, 2.0, 3.0], [2.0, 5.0, 7.0, 7.0])
            }
        )
    elif name == "tvkg11":
        k = {2: grid({(0, 0, 1): 3.0})}
        k33 = {1: grid({(0, 0, 1): 0.2})}
    elif name == "tvkg13":
        k = {1: "TIMEARRAYSERIES ktas"}
        k33 = {1: grid({(0, 0, 1): 0.2})}
        tas = series(
            {
                t: {(0, 0, 1): v}
                for t, v in zip([0.0, 1.0, 2.0, 3.0], [7.0, 3.0, 3.0, 3.0])
            }
        )
    elif name == "tvkg12":
        k = {1: grid({(0, 0, 1): 3.0})}
        k33 = {0: grid({(0, 0, 1): 0.2}), 1: "TIMEARRAYSERIES ktas"}
        tas = series({t: {(0, 0, 3): 0.4} for t in (0.0, 1.0, 2.0, 3.0)})
    elif name == "tvkg10":
        k = {0: "TIMEARRAYSERIES ktas"}
        k33 = None
        tas = series(
            {
                0.0: {(0, 0, 2): 2.0},
                0.5: {(0, 0, 2): DNODATA},
                1.0: {(0, 0, 2): 5.0},
                3.0: {(0, 0, 2): 5.0},
            }
        )
    else:
        k = {1: grid({(0, 0, 2): 5.0})}
        k33 = None
    kwargs = {"filename": f"{name}.tvkg", "print_input": True, "k": k}
    if k33 is not None:
        kwargs["k33"] = k33
    npf.tvkg.initialize(**kwargs)
    if tas is not None:
        npf.tvkg.tas.initialize(
            filename=f"{name}.tvkg.tas",
            tas_array=tas,
            time_series_namerecord="ktas",
            interpolation_methodrecord=method,
        )


def build_models(idx, test):
    name = cases[idx]
    return build_sim(test.workspace, name, "mf6", "tvkg"), None


def check_output(idx, test):
    name = cases[idx]
    if name in xfail:
        text = " ".join((test.workspace / "mfsim.lst").read_text().split())
        assert errors[name] in text, f"expected error not reported: {errors[name]}"
        return

    # reference: TVK list input with the equivalent changes
    ws = test.workspace / "tvk"
    sim = build_sim(ws, name, test.targets["mf6"], "tvk")
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference TVK simulation failed"

    h = flopy.utils.HeadFile(test.workspace / f"{name}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{name}.hds").get_alldata()
    assert np.allclose(h, href, atol=1e-10), (
        f"TVKG heads differ from TVK: max diff {np.abs(h - href).max()}"
    )
    # the K changes must have an effect
    assert np.abs(np.diff(h, axis=0)).max() > 1e-3, "heads do not change"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=name in xfail,
        overwrite=False,
    )
    test.run()
