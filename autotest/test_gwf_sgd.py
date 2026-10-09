"""
Test the specified gradient boundary (SGD) package.

A fully saturated (confined) two-cell vertical column with a constant head
on the top cell and an SGD boundary on the bottom cell. Because the column is
confined the relative permeability is 1, so the SGD boundary removes water at
the head-independent rate

    q = |grad| * K33 * area

directed along the specified gradient (downward). The steady-state budget must
therefore balance the constant-head inflow against this drainage, and the head
drop across the vertical connection must equal |grad| * delz.

Two cases are checked: a unit downward gradient (free drainage) and a gradient
of magnitude two, to confirm the flux scales with the specified gradient.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sgd-free", "sgd-grad2"]
gradz = [-1.0, -2.0]

hk = 1.0  # cm/s
delr = 1.0
delc = 1.0
delz = 1.0
nlay, nrow, ncol = 2, 1, 1
h_top = 10.0
area = delr * delc


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

    flopy.mf6.ModflowGwfic(gwf, strt=h_top)

    # confined: icelltype=0 -> saturated everywhere, krel = 1
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)

    # constant head on the top cell
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), h_top]]},
        pname="CHD-1",
        save_flows=True,
    )

    # specified gradient (free drainage) on the bottom cell
    # columns: cellid gradx grady gradz hwva
    sgd_spd = [[(nlay - 1, 0, 0), 0.0, 0.0, gradz[idx], area]]
    flopy.mf6.ModflowGwfsgd(
        gwf,
        stress_period_data={0: sgd_spd},
        pname="SGD-1",
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

    gmag = abs(gradz[idx])
    q_expected = gmag * hk * area  # drainage (outflow) rate

    # head drop across the vertical connection equals |grad| * delz
    hfile = flopy.utils.HeadFile(ws / f"{gwfname}.hds")
    heads = hfile.get_data(idx=-1).flatten()
    h_bot_expected = h_top - gmag * delz
    assert np.isclose(heads[0], h_top, atol=1e-6), f"top head {heads[0]}"
    assert np.isclose(heads[1], h_bot_expected, atol=1e-6), (
        f"bottom head {heads[1]} expected {h_bot_expected}"
    )

    cbb = flopy.utils.CellBudgetFile(ws / f"{gwfname}.cbc", precision="double")

    # SGD outflow is reported as a negative (leaving the model) value
    sgd = cbb.get_data(text="SGD")[-1]
    q_sgd = float(np.sum(sgd["q"]))
    assert np.isclose(q_sgd, -q_expected, atol=1e-6), (
        f"SGD flow {q_sgd} expected {-q_expected}"
    )

    # constant head inflow must balance the drainage
    chd = cbb.get_data(text="CHD")[-1]
    q_chd = float(np.sum(chd["q"]))
    assert np.isclose(q_chd, q_expected, atol=1e-6), (
        f"CHD flow {q_chd} expected {q_expected}"
    )

    # global mass balance: total in == total out
    assert np.isclose(q_chd + q_sgd, 0.0, atol=1e-6)


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
