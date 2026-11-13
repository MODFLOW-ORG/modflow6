"""
Test that FACTOR does not scale DNODATA (3.0E+30) in grid array
(READARRAYGRID) input. WELG Q is defined at cells A and B, DNODATA
elsewhere, read with FACTOR 2.0 and MAXBOUND 2: the wells are A and B at
twice the given rates.

Case 1 (agfac01): Q read with INTERNAL FACTOR 2.0.
Case 2 (agfac02): Q read with LAYERED, INTERNAL FACTOR 2.0 for each
    layer.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["agfac01", "agfac02"]
nlay, nrow, ncol = 2, 1, 5
cell_a, cell_b, cell_chd = (0, 0, 0), (1, 0, 4), (0, 0, 2)
qa, qb, factor = 10.0, 30.0, 2.0
DNODATA = 3.0e30


def q_rows(layer):
    """Q values of one layer, DNODATA except at cells A and B."""
    q = np.full((nlay, nrow, ncol), DNODATA)
    q[cell_a] = qa
    q[cell_b] = qb
    return " ".join(f"{v:.8e}" for v in q[layer].ravel())


def welg_text(name, layered):
    lines = ["BEGIN options", "  READARRAYGRID", "END options", ""]
    lines += ["BEGIN dimensions", "  MAXBOUND  2", "END dimensions", ""]
    lines += ["BEGIN period  1"]
    if layered:
        lines += ["  q  LAYERED"]
        for k in range(nlay):
            lines += [f"    INTERNAL  FACTOR  {factor}", f"      {q_rows(k)}"]
    else:
        lines += ["  q", f"    INTERNAL  FACTOR  {factor}"]
        lines += [f"      {q_rows(k)}" for k in range(nlay)]
    lines += ["END period  1", ""]
    return "\n".join(lines)


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=10.0, botm=[5.0, 0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[cell_chd, 10.0]]})
    q = np.full((nlay, nrow, ncol), DNODATA)
    q[cell_a] = qa
    q[cell_b] = qb
    flopy.mf6.ModflowGwfwelg(gwf, maxbound=2, q={0: q})
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    sim.write_simulation()
    # replace the package file to set FACTOR (and LAYERED)
    (ws / f"{name}.welg").write_text(welg_text(name, layered=idx == 1))
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    cbc = flopy.utils.CellBudgetFile(test.workspace / f"{name}.cbc", precision="double")
    wel = cbc.get_data(text="WEL")[0]
    rates = {int(n): float(q) for n, q in zip(wel["node"], wel["q"])}
    node_a = np.ravel_multi_index(cell_a, (nlay, nrow, ncol)) + 1
    node_b = np.ravel_multi_index(cell_b, (nlay, nrow, ncol)) + 1
    expected = {node_a: qa * factor, node_b: qb * factor}
    assert rates == pytest.approx(expected), f"WEL rates {rates}, expected {expected}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        overwrite=False,
    )
    test.run()
