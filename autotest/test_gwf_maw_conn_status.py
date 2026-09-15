"""
Tests the MAW CONNECTION_STATUS period setting, which activates and deactivates
an individual well-to-cell connection during a simulation.

A single well is screened over all three layers of the model. Connection data is
specified for the full extent of the well, and individual connections are turned
off and on in later stress periods.

Cases:
  - maw_cs_plug     : the deepest connection is deactivated part way through the
                      simulation and reactivated later, as for a backplugged well.
  - maw_cs_deepen   : the deepest connection is inactive from the start and is
                      activated later, as for a deepened well.
  - maw_cs_override : a connection is deactivated, then the whole well is made
                      INACTIVE and ACTIVE again; the connection must still be
                      inactive afterwards.
  - maw_cs_noop     : every connection is explicitly ACTIVE, which must give the
                      same results as a model with no CONNECTION_STATUS setting.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["maw_cs_plug", "maw_cs_deepen", "maw_cs_override", "maw_cs_noop"]

nlay, nrow, ncol, nper = 3, 1, 5, 3
botm = [-10.0, -20.0, -30.0]
jcol = 2
mawrate = -100.0

# user node numbers of the well connections, which are what the budget file holds
expected_nodes = [k * nrow * ncol + jcol + 1 for k in range(nlay)]

# period data for each case, keyed by case name; None means no setting at all
settings = {
    "maw_cs_plug": {
        0: [[0, "rate", mawrate]],
        1: [[0, "connection_status", 2, "inactive"]],
        2: [[0, "connection_status", 2, "active"]],
    },
    "maw_cs_deepen": {
        0: [[0, "rate", mawrate], [0, "connection_status", 2, "inactive"]],
        1: [[0, "connection_status", 2, "active"]],
    },
    "maw_cs_override": {
        0: [[0, "rate", mawrate], [0, "connection_status", 2, "inactive"]],
        1: [[0, "status", "inactive"]],
        2: [[0, "status", "active"]],
    },
    "maw_cs_noop": {
        0: [[0, "rate", mawrate]]
        + [[0, "connection_status", k, "active"] for k in range(nlay)],
    },
}


def get_model(name, ws, perioddata):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(10.0, 5, 1.0)] * nper
    )
    flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=100.0,
        top=0.0,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=-1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=10.0, k33=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), -1.0], [(0, 0, ncol - 1), -1.0]]
    )

    conndata = [
        [0, k, (k, 0, jcol), 0.0 if k == 0 else botm[k - 1], botm[k], 1.0, 0.1]
        for k in range(nlay)
    ]
    flopy.mf6.ModflowGwfmaw(
        gwf,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        packagedata=[[0, 0.15, botm[-1], -1.0, "THIEM", nlay]],
        connectiondata=conndata,
        perioddata=perioddata,
        budget_filerecord=f"{name}.maw.bud",
        head_filerecord=f"{name}.maw.hds",
        pname="MAW-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sim = get_model(name, test.workspace, settings[name])
    mc = None
    if name == "maw_cs_noop":
        # the same model with no CONNECTION_STATUS setting at all
        mc = get_model(
            name, os.path.join(test.workspace, "mf6"), {0: [[0, "rate", mawrate]]}
        )
    return sim, mc


def gwf_terms(ws, name):
    """GWF budget entries for each stress period, as (nodes, flows)."""
    fpth = os.path.join(ws, f"{name}.maw.bud")
    assert os.path.isfile(fpth)
    cbc = flopy.utils.CellBudgetFile(fpth, precision="double")
    out = []
    for kper in range(nper):
        kk = [x for x in cbc.get_kstpkper() if x[1] == kper][-1]
        rec = cbc.get_data(kstpkper=kk, text="GWF")[0]
        out.append(
            (
                np.array([int(r[1]) for r in rec]),
                np.array([r[2] for r in rec]),
            )
        )
    return out


def check_output(idx, test):
    name = cases[idx]
    terms = gwf_terms(test.workspace, name)

    # the connection list, and its cell numbers, are the same in every period
    for kper, (nodes, _) in enumerate(terms):
        assert np.array_equal(nodes, expected_nodes), (
            f"period {kper + 1} GWF budget nodes {nodes} do not match the well "
            f"connections {expected_nodes}"
        )

    deep = nlay - 1  # the deepest connection
    if name == "maw_cs_plug":
        assert not np.isclose(terms[0][1][deep], 0.0), "connection should be active"
        assert np.isclose(terms[1][1][deep], 0.0), "deactivated connection must be dry"
        assert not np.isclose(terms[2][1][deep], 0.0), "connection should be active"
        # the remaining connections take up the full pumping rate
        assert np.allclose(terms[1][1].sum(), -mawrate, atol=1e-2), (
            f"active connections must carry the well rate, got {terms[1][1].sum()}"
        )
    elif name == "maw_cs_deepen":
        assert np.isclose(terms[0][1][deep], 0.0), "connection must start inactive"
        assert not np.isclose(terms[1][1][deep], 0.0), "connection must be activated"
        assert not np.isclose(terms[2][1][deep], 0.0), "connection must stay active"
    elif name == "maw_cs_override":
        assert np.isclose(terms[0][1][deep], 0.0), "connection must be inactive"
        assert not np.isclose(terms[0][1][0], 0.0), "other connections must flow"
        # an inactive well suspends every connection
        assert np.allclose(terms[1][1], 0.0), "inactive well must have no flow"
        # reactivating the well must not reactivate the connection
        assert np.isclose(terms[2][1][deep], 0.0), (
            "connection must still be inactive after the well is reactivated"
        )
        assert not np.isclose(terms[2][1][0], 0.0), "other connections must resume"
    elif name == "maw_cs_noop":
        base = gwf_terms(os.path.join(test.workspace, "mf6"), name)
        for kper, ((_, q), (_, qbase)) in enumerate(zip(terms, base)):
            assert np.array_equal(q, qbase), (
                f"period {kper + 1} flows differ from the model without "
                f"CONNECTION_STATUS: {q} vs {qbase}"
            )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
