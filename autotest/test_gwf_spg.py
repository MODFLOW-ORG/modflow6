"""
Test the seepage (SPG) boundary package.

The SPG package may only be applied to Richards (UZR) cells, so the model uses
the UZR unsaturated-zone package (newton, non-convertible NPF cells).  The
seepage conductance is approximated LAK-style from the NPF hydraulic
conductivity (K11/K33) scaled by the UZR relative permeability (krel):

    cond = krel * K * HWVA / CL1

A seepage cell discharges water (out of the simulated domain) while its head is
above the cell-center elevation (zero pressure head datum), and becomes a
no-flow boundary when the head drops below it.

The model is a confined two-cell vertical column: a constant head on the top
cell and an SPG seepage boundary (vertical face, IHC=0) on the bottom cell.
Because the column stays fully saturated when discharging, krel = 1 and the
steady state reduces to an analytic two-node Robin problem:

    C_col * (chd - h_bot) = C_spg * (h_bot - z_seep)

with C_col = K33 * area / delz the cell-to-cell conductance and
C_spg = K * HWVA / CL1 the seepage conductance.

Two scenarios are checked:

* active:  the top constant head is above the seepage elevation, so water
  discharges through the seepage cell.
* no-flow: the top constant head is below the seepage elevation, so the
  seepage boundary is inactive and no water leaves through it.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework
from gwf_test_utils import get_uzr_soil_data

cases = ["spg-act", "spg-nof"]
chd_head = [5.0, 0.2]
active_scenario = [True, False]

hk = 1.0
delr = 1.0
delc = 1.0
delz = 1.0
nlay, nrow, ncol = 2, 1, 1
area = delr * delc

# seepage face geometry (bottom cell, vertical face)
spg_cl1 = 0.5  # distance from cell center (0.5) to the bottom face (0.0)
spg_hwva = area  # seepage face area
spg_ihc = 0  # vertical connection
spg_angldegx = 0.0  # reserved, unused

# bottom-cell center elevation (seepage elevation): top=1, bot=0 -> 0.5
z_seep = 0.5

# cell-to-cell and seepage conductances
cond_col = hk * area / delz
cond_spg = hk * spg_hwva / spg_cl1


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
        linear_acceleration="BICGSTAB",
        inner_dvclose=1e-9,
        outer_dvclose=1e-9,
        outer_maximum=200,
        inner_maximum=200,
    )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, newtonoptions="newton"
    )

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

    # confined: icelltype=0 (required by UZR), K drives the conductance
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)

    flopy.mf6.ModflowGwfsto(
        gwf, ss=0.0, sy=0.0, iconvert=0, steady_state=True
    )

    # unsaturated zone Richards flow: flags every cell as a Richards cell
    soil_data = get_uzr_soil_data("Celia1990-eq10-Haverkamp")
    flopy.mf6.ModflowGwfuzr(
        gwf,
        iunsat=1,
        storage_scheme="chord-slope",
        kr_averaging="harmonic",
        soil_model="Haverkamp",
        porosity=soil_data["porosity"],
        satres=soil_data["satres"],
        alphahvk=soil_data["alpha"],
        nhvk=soil_data["n"],
        betahvk=soil_data["beta"],
        khvk=soil_data["k"],
    )

    # constant head on the top cell
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), chd_head[idx]]]},
        pname="CHD-1",
        save_flows=True,
    )

    # seepage boundary on the bottom (Richards) cell
    flopy.mf6.ModflowGwfspg(
        gwf,
        maxbound=1,
        stress_period_data={
            0: [[(nlay - 1, 0, 0), spg_ihc, spg_cl1, spg_hwva, spg_angldegx]]
        },
        pname="SPG-1",
        save_flows=True,
        print_flows=True,
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
    q_spg = float(np.sum(cbb.get_data(text="SPG")[-1]["q"]))
    q_chd = float(np.sum(cbb.get_data(text="CHD")[-1]["q"]))

    if active_scenario[idx]:
        # saturated column (krel = 1): analytic two-node Robin solution
        h_bot = (cond_col * chd_head[idx] + cond_spg * z_seep) / (
            cond_col + cond_spg
        )
        q_expected = cond_spg * (h_bot - z_seep)
        assert h_bot > z_seep, "test setup should discharge seepage"

        assert np.isclose(heads[1], h_bot, atol=1e-5), (
            f"bottom head {heads[1]} expected {h_bot}"
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
        # seepage boundary inactive (head below the seepage elevation)
        assert heads[1] < z_seep, "test setup should keep seepage inactive"
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
