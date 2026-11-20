"""
Test grid array (READARRAYGRID) period input where arrays define
different sets of boundary cells. Values must stay associated with
their cells, boundary cells are the overlap of the main arrays given
in a period block, and undefined (DNODATA) values must not be applied.

Case 1 (agper01): WELG, Q and AUX at cells A and B in period 1. Q
    redefined at cell B only in period 2; AUX carries forward for B.
Case 2 (agper02): WELG, Q and AUX at cells A and B in period 1. AUX-only
    block in period 2 with DNODATA at boundary cell A. Input error.
Case 3 (agper03): WELG, Q at cells A and B, AUX at cell A only. Input
    error.
Case 4 (agper04): WELG, Q and AUX at cell A in period 1. Q at cells A
    and C in period 2 with no AUX for new cell C. Input error.
Case 5 (agper05): WELG, Q and AUX at cells A and B in period 1. AUX-only
    block in period 2 defined at all cells; boundary cells unchanged.
Case 6 (agper06): RIVG, STAGE and RBOT at cells A and B, COND defined at
    all cells. Boundary cells are A and B.
Case 7 (agper07): RIVG, STAGE and RBOT defined at all cells, COND at
    cells A and B. Boundary cells are A and B.
Case 8 (agper08): RIVG, all arrays at cells A and B in period 1. STAGE
    redefined at cell B only in period 2; COND and RBOT carry forward.
Case 9 (agper09): RIVG, RBOT never specified (removed from the flopy
    written file, which fills it with CONSTANT DNODATA). Input error.
Case 10 (agper10): RIVG, non-uniform STAGE defined at all cells (more
    cells than MAXBOUND, listed first), COND at cells A and B. Boundary
    cells are A and B.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "agper01",
    "agper02",
    "agper03",
    "agper04",
    "agper05",
    "agper06",
    "agper07",
    "agper08",
    "agper09",
    "agper10",
]

# expected input error for each failing case
errors = {
    "agper02": "auxiliary variable WEL_MULT is not defined (DNODATA) at "
    "boundary cell (user node 1).",
    "agper03": "auxiliary variable TAG is not defined (DNODATA) at "
    "boundary cell (user node 5).",
    "agper04": "WEL_MULT must be specified for new boundary cell (user node 3).",
    "agper09": "required array RBOT has not been specified.",
}
xfail = set(errors)

nlay, nrow, ncol = 1, 3, 3
delr = delc = 10.0
top = 10.0
botm = [0.0]
strt = 5.0
DNODATA = 3.0e30

# boundary cells (layer, row, col) and their user node numbers
cell_a = (0, 0, 0)
cell_b = (0, 1, 1)
cell_c = (0, 0, 2)
node_a = 1
node_b = 5
node_c = 3
cell_chd = (0, 2, 2)

# river properties
stage_a, stage_b = 8.0, 6.0
cond_a, cond_b = 1.0, 3.0
rbot = 0.0

# non-uniform stage defined at all cells
stage_full = np.linspace(6.0, 8.0, nlay * nrow * ncol).reshape(nlay, nrow, ncol)


def grid_values(values):
    arr = np.full((nlay, nrow, ncol), DNODATA)
    for cell, v in values.items():
        arr[cell] = v
    return arr


def full_values(v):
    return np.full((nlay, nrow, ncol), v)


def add_welg(gwf, idx):
    q1 = grid_values({cell_a: 10.0, cell_b: 30.0})
    mult1 = grid_values({cell_a: 0.5, cell_b: 2.0})
    auxiliary = "WEL_MULT"
    if idx == 0:
        q = {0: q1, 1: grid_values({cell_b: 30.0})}
        aux = {0: [mult1]}
    elif idx == 1:
        q = {0: q1}
        aux = {0: [mult1], 1: [grid_values({cell_b: 2.0})]}
    elif idx == 2:
        q = {0: q1}
        aux = {0: [grid_values({cell_a: 1.0})]}
        auxiliary = "TAG"
    elif idx == 3:
        q = {
            0: grid_values({cell_a: 10.0}),
            1: grid_values({cell_a: 10.0, cell_c: 20.0}),
        }
        aux = {0: [grid_values({cell_a: 0.5})]}
    elif idx == 4:
        q = {0: q1}
        aux = {0: [mult1], 1: [full_values(1.0)]}

    kwargs = dict(q=q, aux=aux, maxbound=2, auxiliary=auxiliary)
    if auxiliary == "WEL_MULT":
        kwargs["auxmultname"] = "WEL_MULT"
    flopy.mf6.ModflowGwfwelg(gwf, **kwargs)


def add_rivg(gwf, idx):
    stage = grid_values({cell_a: stage_a, cell_b: stage_b})
    cond = grid_values({cell_a: cond_a, cell_b: cond_b})
    rb = grid_values({cell_a: rbot, cell_b: rbot})
    kwargs = dict(maxbound=2)
    if idx == 5:
        kwargs.update(stage={0: stage}, cond={0: full_values(cond_a)}, rbot={0: rb})
    elif idx == 6:
        kwargs.update(
            stage={0: full_values(stage_a)},
            cond={0: cond},
            rbot={0: full_values(rbot)},
        )
    elif idx == 7:
        kwargs.update(
            stage={0: stage, 1: grid_values({cell_b: stage_b})},
            cond={0: cond},
            rbot={0: rb},
        )
    elif idx == 8:
        kwargs.update(stage={0: stage}, cond={0: cond})
    elif idx == 9:
        kwargs.update(stage={0: stage_full}, cond={0: cond}, rbot={0: rb})
    flopy.mf6.ModflowGwfrivg(gwf, **kwargs)


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nper = 1 if idx in (2, 5, 6, 8, 9) else 2

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
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[cell_chd, strt]]})

    if idx <= 4:
        add_welg(gwf, idx)
    else:
        add_rivg(gwf, idx)

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    if idx == 8:
        # write input and remove the rbot array flopy fills with DNODATA
        sim.write_simulation()
        fpth = ws / f"{name}.rivg"
        lines = fpth.read_text().splitlines(keepends=True)
        i = next(n for n, line in enumerate(lines) if line.strip() == "rbot")
        fpth.write_text("".join(lines[:i] + lines[i + 2 :]))

    return sim, None


def bnd_rates(test, name, text):
    """Return a {user node: rate} dict for each period."""
    cobj = flopy.utils.CellBudgetFile(
        test.workspace / f"{name}.cbc", precision="double"
    )
    return [
        {int(n): q for n, q in zip(rec["node"], rec["q"])}
        for rec in cobj.get_data(text=text)
    ]


def heads(test, name):
    """Return the head array for each period."""
    hobj = flopy.utils.HeadFile(test.workspace / f"{name}.hds")
    return [hobj.get_data(kstpkper=kk) for kk in hobj.get_kstpkper()]


def check_rates(rates, expected):
    for iper, (r, e) in enumerate(zip(rates, expected)):
        assert r.keys() == e.keys(), (
            f"period {iper + 1}: expected nodes {sorted(e)}, got {sorted(r)}"
        )
        for node, q in e.items():
            assert np.isclose(r[node], q), (
                f"period {iper + 1}: expected rate {q} at node {node}, got {r[node]}"
            )


def riv_expected(h, cells):
    """Return {user node: rate} for river cells {node: (cell, stage, cond)}."""
    return {n: cond * (stage - h[cell]) for n, (cell, stage, cond) in cells.items()}


def check_error(test, name):
    """Check the simulation listing reports the expected input error."""
    text = " ".join((test.workspace / "mfsim.lst").read_text().split())
    assert errors[name] in text, f"expected error not reported: {errors[name]}"


def check_output(idx, test):
    name = cases[idx]
    if name in xfail:
        check_error(test, name)
        return

    if idx == 0:
        expected = [{node_a: 5.0, node_b: 60.0}, {node_b: 60.0}]
        check_rates(bnd_rates(test, name, "WEL"), expected)
    elif idx == 4:
        expected = [{node_a: 5.0, node_b: 60.0}, {node_a: 10.0, node_b: 30.0}]
        check_rates(bnd_rates(test, name, "WEL"), expected)
    else:
        h = heads(test, name)
        both = {
            node_a: (cell_a, stage_a, cond_a),
            node_b: (cell_b, stage_b, cond_b),
        }
        if idx == 5:
            cells = {
                node_a: (cell_a, stage_a, cond_a),
                node_b: (cell_b, stage_b, cond_a),
            }
            expected = [riv_expected(h[0], cells)]
        elif idx == 6:
            cells = {
                node_a: (cell_a, stage_a, cond_a),
                node_b: (cell_b, stage_a, cond_b),
            }
            expected = [riv_expected(h[0], cells)]
        elif idx == 7:
            expected = [
                riv_expected(h[0], both),
                riv_expected(h[1], {node_b: both[node_b]}),
            ]
        elif idx == 9:
            cells = {
                node_a: (cell_a, stage_full[cell_a], cond_a),
                node_b: (cell_b, stage_full[cell_b], cond_b),
            }
            expected = [riv_expected(h[0], cells)]
        check_rates(bnd_rates(test, name, "RIV"), expected)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=name in xfail,
        overwrite=idx != 8,
    )
    test.run()
