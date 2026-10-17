"""
Test the seepage (SPG) boundary package.

A confined two-cell vertical column with a constant head on the top cell and an
SPG seepage boundary on the bottom cell. The seepage cell is held at its
cell-center elevation (atmospheric pressure, zero pressure head) while water
discharges into it (out of the simulated domain), and becomes a no-flow
boundary when the flow would reverse into the domain.

Two enforcement methods are exercised:

* penalty: the active (constant-head) state is approximated with a large
  penalty conductance.
* toggle:  the active state is a true constant-head cell enforced through the
  IBOUND array.

For each method two scenarios are checked:

* active:  the constant head on the top cell is above the seepage elevation, so
  water discharges through the seepage cell. The seepage cell head equals the
  seepage elevation and the discharge balances the constant-head inflow.
* no-flow: the constant head on the top cell is below the seepage elevation, so
  the seepage boundary is inactive and no water leaves through it.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["spg-pen-act", "spg-pen-nof", "spg-tog-act", "spg-tog-nof"]
toggle = [False, False, True, True]
chd_head = [5.0, 0.2, 5.0, 0.2]

# whether the seepage boundary discharges in each case
active_scenario = [True, False, True, False]

hk = 1.0
delr = 1.0
delc = 1.0
delz = 1.0
nlay, nrow, ncol = 2, 1, 1
area = delr * delc

# bottom-cell center elevation (seepage elevation): top=1, bot=0 -> 0.5
z_seep = 0.5
# vertical conductance between the two confined cells
cond = hk * area / delz


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="SECONDS", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        inner_dvclose=1e-9,
        outer_dvclose=1e-9,
        outer_maximum=100,
        inner_maximum=100,
    )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    top = nlay * delz
    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
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

    flopy.mf6.ModflowGwfic(gwf, strt=chd_head[idx])

    # confined: icelltype=0 -> saturated everywhere, constant conductance
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)

    # constant head on the top cell
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), chd_head[idx]]]},
        pname="CHD-1",
        save_flows=True,
    )

    # seepage boundary on the bottom cell
    spg_kwargs = {}
    if toggle[idx]:
        spg_kwargs["ibound_toggle"] = True
    else:
        spg_kwargs["penalty_conductance"] = 1.0e9
    flopy.mf6.ModflowGwfspg(
        gwf,
        maxbound=1,
        stress_period_data={0: [[(nlay - 1, 0, 0)]]},
        pname="SPG-1",
        save_flows=True,
        print_flows=True,
        **spg_kwargs,
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwfname = "gwf_" + name
    ws = test.workspace

    hfile = flopy.utils.HeadFile(ws / f"{gwfname}.hds")
    heads = hfile.get_data(idx=-1).flatten()

    cbb = flopy.utils.CellBudgetFile(ws / f"{gwfname}.cbc", precision="double")
    spg = cbb.get_data(text="SPG")[-1]
    q_spg = float(np.sum(spg["q"]))
    chd = cbb.get_data(text="CHD")[-1]
    q_chd = float(np.sum(chd["q"]))

    if active_scenario[idx]:
        # discharging seepage cell held at the seepage elevation
        q_expected = cond * (chd_head[idx] - z_seep)
        assert np.isclose(heads[1], z_seep, atol=1e-5), (
            f"bottom head {heads[1]} expected {z_seep}"
        )
        # seepage reported as outflow (negative), balancing the chd inflow
        assert np.isclose(q_spg, -q_expected, atol=1e-5), (
            f"SPG flow {q_spg} expected {-q_expected}"
        )
        assert np.isclose(q_chd, q_expected, atol=1e-5), (
            f"CHD flow {q_chd} expected {q_expected}"
        )
        assert np.isclose(q_chd + q_spg, 0.0, atol=1e-5)
    else:
        # seepage boundary inactive: no discharge, head equals the chd head
        assert np.isclose(heads[1], chd_head[idx], atol=1e-6), (
            f"bottom head {heads[1]} expected {chd_head[idx]}"
        )
        assert np.isclose(q_spg, 0.0, atol=1e-6), f"SPG flow {q_spg} expected 0"
        assert np.isclose(q_chd, 0.0, atol=1e-6), f"CHD flow {q_chd} expected 0"


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
