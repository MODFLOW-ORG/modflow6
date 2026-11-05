"""
Test grid array (READARRAYGRID) input with values assigned to cells that
have IDOMAIN values less than one. Such values are not used and must have
no effect on boundaries defined in active cells.

Case 1 (agidom01): WELG, inactive cell between two active wells.
Case 2 (agidom02): CHDG, inactive cell between two active constant heads.
Case 3 (agidom03): WELG, inactive cell introduced in period 2 only.
Case 4 (agidom04): As case 1 with MAXBOUND 2; the inactive cell does not
    count toward MAXBOUND.
Case 5 (agidom05): WELG with AUX. In period 2, Q is defined only at the
    inactive cell and AUX at all cells; there are no wells in period 2.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["agidom01", "agidom02", "agidom03", "agidom04", "agidom05"]

nlay, nrow, ncol = 1, 3, 3
delr = delc = 10.0
top = 10.0
botm = [0.0]
strt = 5.0
DNODATA = 3.0e30

# inactive center cell
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, 1, 1] = 0

# boundary cells (layer, row, col) in user node order; center is inactive
cell_first = (0, 0, 0)
cell_inactive = (0, 1, 1)
cell_last = (0, 2, 2)


def grid_values(values):
    arr = np.full((nlay, nrow, ncol), DNODATA)
    for cell, v in values.items():
        arr[cell] = v
    return arr


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    nper = 2 if idx in (2, 4) else 1

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
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)

    if idx in (0, 3):
        q = {0: grid_values({cell_first: 10.0, cell_inactive: 20.0, cell_last: 30.0})}
        flopy.mf6.ModflowGwfwelg(gwf, q=q, maxbound=3 if idx == 0 else 2)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 2), strt]]})
    elif idx == 1:
        head = {0: grid_values({cell_first: 8.0, cell_inactive: 99.0, cell_last: 2.0})}
        flopy.mf6.ModflowGwfchdg(gwf, head=head, maxbound=3)
    elif idx == 2:
        q = {
            0: grid_values({cell_first: 10.0, cell_last: 30.0}),
            1: grid_values({cell_first: 10.0, cell_inactive: 20.0, cell_last: 30.0}),
        }
        flopy.mf6.ModflowGwfwelg(gwf, q=q, maxbound=3)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 2), strt]]})
    elif idx == 4:
        q = {
            0: grid_values({cell_first: 10.0, cell_last: 30.0}),
            1: grid_values({cell_inactive: 20.0}),
        }
        aux = {
            0: [grid_values({cell_first: 0.5, cell_last: 2.0})],
            1: [np.ones((nlay, nrow, ncol))],
        }
        flopy.mf6.ModflowGwfwelg(
            gwf, q=q, aux=aux, auxiliary="MULT", auxmultname="MULT", maxbound=2
        )
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 2), strt]]})

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def wel_rates(test, name):
    cobj = flopy.utils.CellBudgetFile(
        test.workspace / f"{name}.cbc", precision="double"
    )
    return [np.sort(rec["q"]) for rec in cobj.get_data(text="WEL")]


def check_output(idx, test):
    name = cases[idx]

    if idx in (0, 3):
        rates = wel_rates(test, name)
        assert np.allclose(rates[0], [10.0, 30.0]), (
            f"expected well rates [10, 30], got {rates[0]}"
        )
    elif idx == 1:
        hobj = flopy.utils.HeadFile(test.workspace / f"{name}.hds")
        head = hobj.get_data()
        assert np.isclose(head[cell_first], 8.0), (
            f"expected head 8.0 at {cell_first}, got {head[cell_first]}"
        )
        assert np.isclose(head[cell_last], 2.0), (
            f"expected head 2.0 at {cell_last}, got {head[cell_last]}"
        )
    elif idx == 2:
        rates = wel_rates(test, name)
        for iper, r in enumerate(rates):
            assert np.allclose(r, [10.0, 30.0]), (
                f"period {iper + 1}: expected well rates [10, 30], got {r}"
            )
    elif idx == 4:
        rates = wel_rates(test, name)
        assert np.allclose(rates[0], [5.0, 60.0]), (
            f"period 1: expected well rates [5, 60], got {rates[0]}"
        )
        assert rates[1].size == 0, f"period 2: expected no wells, got {rates[1]}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
