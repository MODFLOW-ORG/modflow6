"""Regression test for PACKAGEDATA CELLID's documented "unconnected reach"
syntax: both the literal all-zero cellid ("0 0 0" for a DIS grid, the
documented/preferred form) and the deprecated "NONE" keyword alias.

flopy's cellid API cannot write either sentinel directly -- a (0, 0, 0)
tuple is treated as a real cell reference and gets its usual 0-based-to-
1-based conversion applied, and a raw "0 0 0" string is rejected by
flopy's formatter. Both cases write a placeholder cellid via flopy, then
hand-patch it into the file afterward.

Confirms the unconnected reach is absent from the SFR GWF-exchange budget
record while its connected neighbors report nonzero exchange, matching
the documented "reach-aquifer flow is not calculated for unconnected
reaches".
"""

import re

import flopy
import pytest
from framework import TestFramework

cases = ["sfr-unconn-000", "sfr-unconn-none"]

# distinct, valid placeholder cellid for the middle reach; patched to the
# true sentinel after write_simulation()
_placeholder_cellid = (0, 0, 4)


def _patch_cellid(path, reach_num, replacement):
    """Replace reach_num's (1-based) 3-token CELLID with replacement."""
    with open(path) as f:
        lines = f.readlines()
    pattern = re.compile(rf"^(\s*{reach_num}\s+)(\S+\s+\S+\s+\S+)(\s+.*)$")
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if m:
            lines[i] = m.group(1) + replacement + m.group(3) + "\n"
            break
    else:
        raise RuntimeError(f"reach {reach_num} not found in {path}")
    with open(path, "w") as f:
        f.writelines(lines)


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 4), 4.0]]
    )

    # middle reach (IFNO 1) gets a placeholder cellid here; patched below to
    # the true unconnected sentinel, which flopy's API cannot express
    packagedata = [
        [0, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
        [1, _placeholder_cellid, 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 2, 1.0, 0],
        [2, (0, 0, 3), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
    ]
    connectiondata = [[0, -1], [1, 0, -2], [2, 1]]
    perioddata = {0: [[0, "inflow", 0.1]]}

    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
        budget_filerecord=f"{name}.sfr.bud",
        save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )

    sim.write_simulation(silent=True)
    sentinel = "0 0 0" if idx == 0 else "NONE"
    _patch_cellid(test.workspace / f"{name}.sfr", reach_num=2, replacement=sentinel)

    return sim


def check_output(idx, test):
    name = cases[idx]
    sfrbud = flopy.utils.CellBudgetFile(
        str(test.workspace / f"{name}.sfr.bud"), precision="double"
    )
    gwf_flow = sfrbud.get_data(text="GWF")[0]
    # the budget record is sparse (connected reaches only); the unconnected
    # reach (node number 2, the placeholder's real cell) must not appear
    assert 2 not in gwf_flow["node"], (
        f"unconnected reach should be absent from the sparse GWF budget "
        f"record entirely, got nodes {gwf_flow['node'].tolist()}"
    )
    assert len(gwf_flow) == 2 and all(gwf_flow["q"] != 0.0), (
        "both connected reaches should report nonzero GWF exchange, got "
        f"{gwf_flow.tolist()}"
    )


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
