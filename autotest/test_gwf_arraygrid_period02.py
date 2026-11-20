"""
Test grid array (READARRAYGRID) period input loader paths not covered by
test_gwf_arraygrid_period.py. Most cases edit the flopy written package
file to set array order, empty blocks, or keywords flopy does not write.

Case 11 (agper11): WELG, Q and AUX at cells A and B in period 1, empty
    block in period 2, Q only in period 3. Input error (empty block
    resets; AUX must be specified again).
Case 12 (agper12): As case 11, with Q and AUX both specified in period 3.
Case 13 (agper13): WELG, AUX listed before Q, both at cells A and B.
Case 14 (agper14): WELG, AUX defined at all cells (more cells than
    MAXBOUND) listed before Q at cells A and B.
Case 15 (agper15): RIVG, every main array defined at more cells than
    MAXBOUND; overlap also exceeds MAXBOUND. Input error.
Case 16 (agper16): RIVG, every main array defined at more cells than
    MAXBOUND; overlap fits MAXBOUND.
Case 17 (agper17): RIVG, STAGE at cells A, C, B, then AUX multiplier,
    then COND at cells A and B. COND narrows the cells after AUX is
    stored.
Case 18 (agper18): WELG, Q and AUX read with the LAYERED keyword.
Case 19 (agper19): WELG, Q listed twice in a block. Input error.
Case 20 (agper20): WELG, AUX array given with the literal tag AUX instead
    of the auxiliary variable name. Input error.
Case 21 (agper21): WELG, Q at cells A and B read with FACTOR 2.0;
    DNODATA values are not scaled.
Case 22 (agper22): RIVG, STAGE at cells A and B, COND at cell C only;
    the main arrays share no defined cells. Input error.
Case 23 (agper23): WELG, Q at cells A and B, WEL_MULT at cell A only.
    Input error (an aux array must be defined at every boundary cell).
Case 24 (agper24): RIVG, STAGE and RBOT defined, COND not defined at any
    cell. Input error.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwf_arraygrid_period import (
    DNODATA,
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
    nrow,
    rbot,
    riv_expected,
    stage_a,
    stage_b,
    stage_full,
    strt,
    top,
)

cases = [
    "agper11",
    "agper12",
    "agper13",
    "agper14",
    "agper15",
    "agper16",
    "agper17",
    "agper18",
    "agper19",
    "agper20",
    "agper21",
    "agper22",
    "agper23",
    "agper24",
]

# expected input error for each failing case
errors = {
    "agper11": "WEL_MULT must be specified for new boundary cell (user node 1).",
    "agper15": "number of defined (non-DNODATA) active cells exceeds MAXBOUND=6.",
    "agper19": "Q is specified more than once in the PERIOD block.",
    "agper20": "AUX is not a valid PERIOD block array name.",
    "agper22": "main arrays have no defined (non-DNODATA) values in the same "
    "active cells.",
    "agper23": "auxiliary variable WEL_MULT is not defined (DNODATA) at boundary "
    f"cell (user node {node_b}).",
    "agper24": "COND has no defined (non-DNODATA) values in active cells.",
}
xfail = set(errors)

# all cells except one, so each array is defined at more cells than
# MAXBOUND while the overlap excludes two cells
cell_last = (0, 2, 2)


def all_but(cell, v):
    arr = full_values(v) if np.isscalar(v) else v.copy()
    arr[cell] = DNODATA
    return arr


def split_block(lines):
    """Split period block lines into [(tag line, data lines)]."""
    arrays = []
    for line in lines:
        if line.startswith("  ") and not line.startswith("   "):
            arrays.append([line, []])
        else:
            arrays[-1][1].append(line)
    return arrays


def edit_period(text, iper, edit):
    """Apply edit(arrays) to the arrays of period block iper."""
    lines = text.splitlines(keepends=True)
    begin = next(
        i for i, ln in enumerate(lines) if ln.split() == ["BEGIN", "period", str(iper)]
    )
    end = next(
        i for i in range(begin, len(lines)) if lines[i].split()[:2] == ["END", "period"]
    )
    arrays = edit(split_block(lines[begin + 1 : end]))
    body = [ln for tag, data in arrays for ln in [tag, *data]]
    return "".join(lines[: begin + 1] + body + lines[end:])


def insert_empty_period(text, iper):
    """Insert an empty period block before period block iper + 1."""
    marker = f"BEGIN period  {iper + 1}\n"
    assert marker in text
    return text.replace(marker, f"BEGIN period  {iper}\nEND period  {iper}\n\n{marker}")


def move_first(tag):
    """Return an edit that moves array tag to the front of the block."""

    def edit(arrays):
        i = next(n for n, (t, _) in enumerate(arrays) if t.strip() == tag)
        return [arrays[i]] + arrays[:i] + arrays[i + 1 :]

    return edit


def rename(tag, new):
    """Return an edit that renames array tag."""

    def edit(arrays):
        return [(f"  {new}\n" if t.strip() == tag else t, d) for t, d in arrays]

    return edit


def check_error(test, name):
    """Check the simulation listing reports the expected input error."""
    text = " ".join((test.workspace / "mfsim.lst").read_text().split())
    assert errors[name] in text, f"expected error not reported: {errors[name]}"


def add_layered(arrays):
    return [(t.rstrip("\n") + " LAYERED\n", d) for t, d in arrays]


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nper = 3 if idx in (0, 1) else 1

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-9, inner_dvclose=1e-9)
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    if idx not in (4, 5):
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[cell_chd, strt]]})

    q1 = grid_values({cell_a: 10.0, cell_b: 30.0})
    mult1 = grid_values({cell_a: 0.5, cell_b: 2.0})
    welg = dict(auxiliary="WEL_MULT", auxmultname="WEL_MULT", maxbound=2)
    if idx == 0:
        welg.update(q={0: q1, 2: q1}, aux={0: [mult1]})
    elif idx == 1:
        welg.update(q={0: q1, 2: q1}, aux={0: [mult1], 2: [mult1]})
    elif idx in (2, 7, 8, 9):
        welg.update(q={0: q1}, aux={0: [mult1]})
    elif idx == 3:
        welg.update(q={0: q1}, aux={0: [full_values(2.0)]})
    elif idx == 10:
        welg = dict(maxbound=2, q={0: q1})
    elif idx == 12:
        welg.update(q={0: q1}, aux={0: [grid_values({cell_a: 0.5})]})

    if idx in (0, 1, 2, 3, 7, 8, 9, 10, 12):
        flopy.mf6.ModflowGwfwelg(gwf, **welg)
        ext = "welg"
    elif idx in (4, 5):
        # stage undefined at cell A, cond undefined at the last cell
        flopy.mf6.ModflowGwfrivg(
            gwf,
            maxbound=6 if idx == 4 else 7,
            stage={0: all_but(cell_a, stage_full)},
            cond={0: all_but(cell_last, cond_a)},
            rbot={0: full_values(rbot)},
        )
        ext = "rivg"
    elif idx == 6:
        # cells in user node order: A (1), C (3), B (5)
        flopy.mf6.ModflowGwfrivg(
            gwf,
            maxbound=3,
            auxiliary="MULT",
            auxmultname="MULT",
            stage={0: grid_values({cell_a: stage_a, cell_c: 7.0, cell_b: stage_b})},
            cond={0: grid_values({cell_a: cond_a, cell_b: cond_b})},
            rbot={0: full_values(rbot)},
            aux={0: [grid_values({cell_a: 1.0, cell_c: 3.0, cell_b: 2.0})]},
        )
        ext = "rivg"
    elif idx in (11, 13):
        cond = grid_values({cell_c: cond_a}) if idx == 11 else grid_values({})
        flopy.mf6.ModflowGwfrivg(
            gwf,
            maxbound=3,
            stage={0: grid_values({cell_a: stage_a, cell_b: stage_b})},
            cond={0: cond},
            rbot={0: full_values(rbot)},
        )
        ext = "rivg"

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
    if idx in (0, 1):
        text = insert_empty_period(text, 2)
    elif idx in (2, 3):
        text = edit_period(text, 1, move_first("WEL_MULT"))
    elif idx == 6:
        text = edit_period(text, 1, move_first("MULT"))
        text = edit_period(text, 1, move_first("stage"))
    elif idx == 7:
        text = edit_period(text, 1, add_layered)
    elif idx == 8:
        text = edit_period(text, 1, lambda arrays: arrays[:1] + arrays)
    elif idx == 9:
        text = edit_period(text, 1, rename("WEL_MULT", "AUX"))
    elif idx == 10:
        assert "INTERNAL  FACTOR  1.0" in text
        text = text.replace("INTERNAL  FACTOR  1.0", "INTERNAL  FACTOR  2.0")
    fpth.write_text(text)

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    if name in xfail:
        check_error(test, name)
        return

    if idx == 1:
        expected = [{node_a: 5.0, node_b: 60.0}, {}, {node_a: 5.0, node_b: 60.0}]
        check_rates(bnd_rates(test, name, "WEL"), expected)
    elif idx in (2, 7):
        check_rates(bnd_rates(test, name, "WEL"), [{node_a: 5.0, node_b: 60.0}])
    elif idx == 3:
        check_rates(bnd_rates(test, name, "WEL"), [{node_a: 20.0, node_b: 60.0}])
    elif idx == 10:
        check_rates(bnd_rates(test, name, "WEL"), [{node_a: 20.0, node_b: 60.0}])
    else:
        h = heads(test, name)[0]
        if idx == 5:
            cells = {
                n + 1: (
                    np.unravel_index(n, (nlay, nrow, ncol)),
                    stage_full.flat[n],
                    cond_a,
                )
                for n in range(1, nlay * nrow * ncol - 1)
            }
        elif idx == 6:
            cells = {
                node_a: (cell_a, stage_a, cond_a * 1.0),
                node_b: (cell_b, stage_b, cond_b * 2.0),
            }
        check_rates(bnd_rates(test, name, "RIV"), [riv_expected(h, cells)])


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
