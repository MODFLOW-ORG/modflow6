"""
Test grid array (READARRAYGRID) period input given as full grid
time-array series (TIMEARRAYSERIES). A linked array keeps progressing
through periods that do not give it, is replaced when it is given again,
and sets boundary cells from the series array at the period start. A
series value of DNODATA at a boundary cell is an input error.

Case 1 (agtas01): WELG, Q series at cells A and B (stepwise) in period
    1; no block in period 2, so the series keeps progressing.
Case 2 (agtas02): WELG, Q series (linear), two time steps per period;
    rates are time step averages.
Case 3 (agtas03): WELG, Q series given in period 1 and again in period
    2, static Q in period 3, empty block in period 4, Q series in
    period 5.
Case 4 (agtas04): WELG, static Q with an AUX multiplier series.
Case 5 (agtas05): RIVG, STAGE series defined at all cells (more cells
    than MAXBOUND) narrowed by static COND at cells A and B.
Case 6 (agtas06): WELG, Q series at cells A and B; the later series
    array also defines cell C, which is not a boundary and is ignored.
Case 7 (agtas07): RIVG, STAGE series in period 1; static COND adds new
    cell C in period 2 while STAGE is carried. Input error.
Case 8 (agtas08): WELG, Q series DNODATA at boundary cell B at a later
    time. Input error.
Case 9 (agtas09): WELG, TIMEARRAYSERIES without a TAS6 file. Input
    error.
Case 10 (agtas10): As case 1 with a DISU grid.
Case 11 (agtas11): As case 1 with FACTOR 2.0 on the series arrays;
    DNODATA values are not scaled.
Case 12 (agtas12): As case 7, with STAGE given again with the same series
    in period 2, so new cell C gets the series value.
Case 13 (agtas13): WELG, Q series in period 1 replaced by a different
    series in period 2.
Case 14 (agtas14): WELG, Q series given in periods 1 and 4 of four 0.7
    day periods; the series record at 2.1 adds cell B. The period 4 start
    time has roundoff (2.0999999999999996), so the record must be found
    by a close time match.
Case 15 (agtas15): As case 1 with SFAC 2.0 in the series attributes;
    DNODATA values are not scaled.
"""

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework
from test_gwf_arraygrid_period import (
    bnd_rates,
    botm,
    cell_a,
    cell_b,
    cell_c,
    cell_chd,
    check_rates,
    cond_a,
    cond_b,
    delc,
    delr,
    full_values,
    grid_values,
    heads,
    ncol,
    nlay,
    node_a,
    node_b,
    node_c,
    nrow,
    riv_expected,
    stage_full,
    strt,
    top,
)
from test_gwf_arraygrid_period02 import edit_period, insert_empty_period

cases = [
    "agtas01",
    "agtas02",
    "agtas03",
    "agtas04",
    "agtas05",
    "agtas06",
    "agtas07",
    "agtas08",
    "agtas09",
    "agtas10",
    "agtas11",
    "agtas12",
    "agtas13",
    "agtas14",
    "agtas15",
]

# expected input error for each failing case
errors = {
    "agtas07": "STAGE must be specified for new boundary cell (user node 3).",
    "agtas08": 'Time-array series "QTS" is not defined (DNODATA) at user node 5',
    "agtas09": "Q TIMEARRAYSERIES requires a TAS6 file in the OPTIONS block.",
}
xfail = set(errors)

perlen = {"agtas14": 0.7}
nper = {"agtas03": 5, "agtas14": 4}
nstp = {"agtas02": 2}
tas_q = "TIMEARRAYSERIES qts"


def ab(a, b, c=None):
    """Grid values at cells A and B (and C)."""
    values = {cell_a: a, cell_b: b}
    if c is not None:
        values[cell_c] = c
    return grid_values(values)


def add_tas(pkg, name, method, records):
    """Add a full grid time-array series file to a grid package."""
    pkg.tas.initialize(
        filename=f"{pkg.filename}.tas",
        tas_array=records,
        time_series_namerecord=name,
        interpolation_methodrecord=method,
    )


def tas_tag(tag, tas_name):
    """Return an edit that gives array tag as a time-array series."""

    def edit(arrays):
        return [
            (f"  {tag} TIMEARRAYSERIES {tas_name}\n", [])
            if t.strip() == tag
            else (t, d)
            for t, d in arrays
        ]

    return edit


def check_error(test, name):
    """Check the simulation listing reports the expected input error."""
    text = " ".join((test.workspace / "mfsim.lst").read_text().split())
    assert errors[name] in text, f"expected error not reported: {errors[name]}"


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    np_ = nper.get(name, 2)
    ns = nstp.get(name, 1)

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    pl = perlen.get(name, 10.0)
    flopy.mf6.ModflowTdis(sim, nper=np_, perioddata=[(pl, ns, 1.0)] * np_)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-9, inner_dvclose=1e-9)
    if name == "agtas10":
        disu = get_disu_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
        flopy.mf6.ModflowGwfdisu(gwf, **disu)
        chd_cell = (int(np.ravel_multi_index(cell_chd, (nlay, nrow, ncol))),)
    else:
        flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        chd_cell = cell_chd
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    if name not in ("agtas05", "agtas07", "agtas12"):
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[chd_cell, strt]]})

    ext = "welg"
    if name in ("agtas01", "agtas09", "agtas10", "agtas11", "agtas15"):
        pkg = flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q={0: tas_q})
        if name != "agtas09":
            records = {0.0: ab(10.0, 20.0), 10.0: ab(30.0, 40.0), 20.0: ab(30.0, 40.0)}
            if name == "agtas10":
                records = {t: v.ravel() for t, v in records.items()}
            add_tas(pkg, "qts", "stepwise", records)
    elif name == "agtas14":
        pkg = flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q={0: tas_q, 3: tas_q})
        a_only = grid_values({cell_a: 10.0})
        records = {0.0: a_only, 2.1: ab(10.0, 40.0), 10.0: ab(10.0, 40.0)}
        add_tas(pkg, "qts", "stepwise", records)
    elif name == "agtas02":
        pkg = flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q={0: tas_q})
        add_tas(pkg, "qts", "linear", {0.0: ab(0.0, 0.0), 20.0: ab(20.0, 40.0)})
    elif name == "agtas03":
        q = {0: tas_q, 1: tas_q, 2: ab(7.0, 8.0), 4: tas_q}
        pkg = flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q=q)
        records = {10.0 * k: ab(k + 1.0, 10.0 * (k + 1)) for k in range(6)}
        add_tas(pkg, "qts", "stepwise", records)
    elif name == "agtas04":
        pkg = flopy.mf6.ModflowGwfwelg(
            gwf,
            maxbound=2,
            auxiliary="WEL_MULT",
            auxmultname="WEL_MULT",
            q={0: ab(10.0, 20.0)},
            aux={0: [ab(1.0, 1.0)]},
        )
        records = {0.0: ab(0.5, 2.0), 10.0: ab(1.0, 3.0), 20.0: ab(1.0, 3.0)}
        add_tas(pkg, "multts", "stepwise", records)
    elif name in ("agtas06", "agtas08"):
        pkg = flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q={0: tas_q})
        if name == "agtas06":
            later = ab(30.0, 40.0, 50.0)
        else:
            later = grid_values({cell_a: 30.0})
        records = {0.0: ab(10.0, 20.0), 10.0: later, 20.0: later}
        add_tas(pkg, "qts", "stepwise", records)
    elif name == "agtas13":
        pkg = flopy.mf6.ModflowGwfwelg(
            gwf, maxbound=2, q={0: tas_q, 1: "TIMEARRAYSERIES qts2"}
        )
        add_tas(
            pkg,
            "qts",
            "stepwise",
            {0.0: ab(10.0, 20.0), 10.0: ab(30.0, 40.0), 20.0: ab(30.0, 40.0)},
        )
        pkg.tas.append_package(
            filename=f"{pkg.filename}.2.tas",
            tas_array={0.0: ab(50.0, 60.0), 20.0: ab(50.0, 60.0)},
            time_series_namerecord="qts2",
            interpolation_methodrecord="stepwise",
        )
    elif name in ("agtas05", "agtas07", "agtas12"):
        ext = "rivg"
        cond = {0: ab(cond_a, cond_b)}
        stage = {0: "TIMEARRAYSERIES stagets"}
        rbot = {0: full_values(0.0)}
        if name in ("agtas07", "agtas12"):
            cond[1] = ab(cond_a, cond_b, cond_a)
        if name == "agtas12":
            stage[1] = "TIMEARRAYSERIES stagets"
            rbot[1] = full_values(0.0)
        pkg = flopy.mf6.ModflowGwfrivg(
            gwf,
            maxbound=2 if name == "agtas05" else 3,
            stage=stage,
            cond=cond,
            rbot=rbot,
        )
        records = {0.0: stage_full, 10.0: stage_full + 1.0, 20.0: stage_full + 1.0}
        add_tas(pkg, "stagets", "stepwise", records)

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # write input and edit the package file
    sim.write_simulation()
    fpth = ws / f"{name}.{ext}"
    text = fpth.read_text()
    if name == "agtas03":
        text = insert_empty_period(text, 4)
    elif name == "agtas04":
        text = edit_period(text, 1, tas_tag("WEL_MULT", "multts"))
    fpth.write_text(text)
    if name == "agtas11":
        fpth = ws / f"{name}.welg.tas"
        text = fpth.read_text()
        assert "INTERNAL  FACTOR  1.0" in text
        fpth.write_text(text.replace("INTERNAL  FACTOR  1.0", "INTERNAL  FACTOR  2.0"))
    elif name == "agtas15":
        fpth = ws / f"{name}.welg.tas"
        text = fpth.read_text()
        assert text.count("END attributes") == 1
        fpth.write_text(text.replace("END attributes", "  SFAC  2.0\nEND attributes"))

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    if name in xfail:
        check_error(test, name)
        return

    if name in ("agtas01", "agtas10"):
        expected = [{node_a: 10.0, node_b: 20.0}, {node_a: 30.0, node_b: 40.0}]
    elif name in ("agtas11", "agtas15"):
        expected = [{node_a: 20.0, node_b: 40.0}, {node_a: 60.0, node_b: 80.0}]
    elif name == "agtas14":
        expected = [{node_a: 10.0}] * 3 + [{node_a: 10.0, node_b: 40.0}]
    elif name == "agtas02":
        # linear Q = t at A and 2t at B, 5 day time steps
        expected = [
            {node_a: 2.5, node_b: 5.0},
            {node_a: 7.5, node_b: 15.0},
            {node_a: 12.5, node_b: 25.0},
            {node_a: 17.5, node_b: 35.0},
        ]
    elif name == "agtas03":
        expected = [
            {node_a: 1.0, node_b: 10.0},
            {node_a: 2.0, node_b: 20.0},
            {node_a: 7.0, node_b: 8.0},
            {},
            {node_a: 5.0, node_b: 50.0},
        ]
    elif name == "agtas04":
        expected = [{node_a: 5.0, node_b: 40.0}, {node_a: 10.0, node_b: 60.0}]
    elif name == "agtas06":
        expected = [{node_a: 10.0, node_b: 20.0}, {node_a: 30.0, node_b: 40.0}]
    elif name == "agtas13":
        expected = [{node_a: 10.0, node_b: 20.0}, {node_a: 50.0, node_b: 60.0}]

    if name in ("agtas05", "agtas12"):
        h = heads(test, name)
        cells = [
            {
                node_a: (cell_a, stage_full[cell_a] + k, cond_a),
                node_b: (cell_b, stage_full[cell_b] + k, cond_b),
            }
            for k in range(2)
        ]
        if name == "agtas12":
            cells[1][node_c] = (cell_c, stage_full[cell_c] + 1.0, cond_a)
        expected = [riv_expected(h[k], cells[k]) for k in range(2)]
        check_rates(bnd_rates(test, name, "RIV"), expected)
    else:
        check_rates(bnd_rates(test, name, "WEL"), expected)


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
