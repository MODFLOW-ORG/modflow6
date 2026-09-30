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
  - maw_cs_constant : as maw_cs_override, with the well made CONSTANT instead
                      of ACTIVE.
  - maw_cs_restore  : a connection is reactivated while the well is INACTIVE;
                      it must flow once the well is made ACTIVE again.
  - maw_cs_bottom   : the well bottom is raised to the bottom of the middle
                      connection while the deepest connection is deactivated in
                      the same period; BOTTOM is listed first and the result must
                      match the same settings in the opposite order.
  - maw_cs_noop     : every connection is explicitly ACTIVE, which must give the
                      same results as a model with no CONNECTION_STATUS setting.

Separate tests check that making every connection of a well INACTIVE is an
error, since STATUS INACTIVE is used to deactivate a well, and that a BOTTOM
above the head or above the screen bottom of an active connection is an error
(xfail).
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "maw_cs_plug",
    "maw_cs_deepen",
    "maw_cs_override",
    "maw_cs_constant",
    "maw_cs_restore",
    "maw_cs_bottom",
    "maw_cs_noop",
]

nlay, nrow, ncol, nper = 3, 1, 5, 3
botm = [-10.0, -20.0, -30.0]
jcol = 2
mawrate = -100.0
dnodata = 3.0e30

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
    "maw_cs_constant": {
        0: [[0, "rate", mawrate], [0, "connection_status", 2, "inactive"]],
        1: [[0, "status", "inactive"]],
        2: [[0, "status", "constant"], [0, "well_head", -5.0]],
    },
    "maw_cs_restore": {
        0: [[0, "rate", mawrate], [0, "connection_status", 2, "inactive"]],
        1: [[0, "status", "inactive"], [0, "connection_status", 2, "active"]],
        2: [[0, "status", "active"]],
    },
    "maw_cs_bottom": {
        0: [[0, "rate", mawrate]],
        1: [[0, "bottom", botm[1]], [0, "connection_status", 2, "inactive"]],
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
    maw = flopy.mf6.ModflowGwfmaw(
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
    # flow and conductance of the shallowest and deepest connections
    maw.obs.initialize(
        filename=f"{name}.maw.obs",
        continuous={
            f"{name}.maw.obs.csv": [
                ("q1", "maw", 1, 1),
                ("q3", "maw", 1, nlay),
                ("c3", "conductance", 1, nlay),
            ]
        },
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
    if name == "maw_cs_bottom":
        # the same settings with BOTTOM after CONNECTION_STATUS
        perioddata = dict(settings[name])
        perioddata[1] = list(reversed(perioddata[1]))
        mc = get_model(name, os.path.join(test.workspace, "mf6"), perioddata)
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


def period_obs(ws, name):
    """MAW observations at the end of each stress period, keyed by name."""
    obs = flopy.utils.Mf6Obs(os.path.join(ws, f"{name}.maw.obs.csv")).get_data()
    tend = [10.0 * (kper + 1) for kper in range(nper)]
    rows = [obs[np.isclose(obs["totim"], t)][0] for t in tend]
    return {k: np.array([r[k] for r in rows]) for k in ("Q1", "Q3", "C3")}


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
        # an inactive connection reports DNODATA, as an inactive well does
        obs = period_obs(test.workspace, name)
        for key in ("Q3", "C3"):
            assert obs[key][1] == dnodata, f"{key} must be DNODATA when inactive"
            assert obs[key][0] < dnodata and obs[key][2] < dnodata, (
                f"{key} must be simulated when the connection is active"
            )
        # the remaining connections take up the full pumping rate
        assert np.allclose(terms[1][1].sum(), -mawrate, atol=1e-2), (
            f"active connections must carry the well rate, got {terms[1][1].sum()}"
        )
    elif name == "maw_cs_deepen":
        assert np.isclose(terms[0][1][deep], 0.0), "connection must start inactive"
        assert not np.isclose(terms[1][1][deep], 0.0), "connection must be activated"
        assert not np.isclose(terms[2][1][deep], 0.0), "connection must stay active"
    elif name in ("maw_cs_override", "maw_cs_constant"):
        assert np.isclose(terms[0][1][deep], 0.0), "connection must be inactive"
        assert not np.isclose(terms[0][1][0], 0.0), "other connections must flow"
        # an inactive well suspends every connection
        assert np.allclose(terms[1][1], 0.0), "inactive well must have no flow"
        # reactivating the well must not reactivate the connection
        assert np.isclose(terms[2][1][deep], 0.0), (
            "connection must still be inactive after the well is reactivated"
        )
        assert not np.isclose(terms[2][1][0], 0.0), "other connections must resume"
        obs = period_obs(test.workspace, name)
        assert np.all(obs["Q3"] == dnodata), "inactive connection must be DNODATA"
        assert obs["Q1"][1] == dnodata, "an inactive well must report DNODATA"
        assert obs["Q1"][0] < dnodata and obs["Q1"][2] < dnodata, (
            "an active connection of an active well must be simulated"
        )
    elif name == "maw_cs_restore":
        assert np.isclose(terms[0][1][deep], 0.0), "connection must start inactive"
        # the inactive well suspends the reactivated connection
        assert np.allclose(terms[1][1], 0.0), "inactive well must have no flow"
        assert not np.isclose(terms[2][1][deep], 0.0), (
            "connection reactivated while the well was inactive must flow"
        )
    elif name == "maw_cs_bottom":
        base = gwf_terms(os.path.join(test.workspace, "mf6"), name)
        for kper, ((_, q), (_, qbase)) in enumerate(zip(terms, base)):
            assert np.array_equal(q, qbase), (
                f"period {kper + 1} flows depend on the order of BOTTOM and "
                f"CONNECTION_STATUS: {q} vs {qbase}"
            )
        # the raised bottom and the inactive connection persist
        for kper in (1, 2):
            assert np.isclose(terms[kper][1][deep], 0.0), "connection must be dry"
            assert not np.isclose(terms[kper][1][0], 0.0), "others must flow"
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


def test_all_connections_inactive(function_tmpdir, targets):
    # deactivating every connection of a well is an error; the well is
    # deactivated with STATUS INACTIVE instead
    name = "maw_cs_allinact"
    perioddata = {
        0: [[0, "rate", mawrate]],
        1: [[0, "connection_status", k, "inactive"] for k in range(nlay)],
    }

    def check(test):
        text = open(os.path.join(test.workspace, "mfsim.lst")).read()
        assert "Every connection of maw well 1 is inactive" in text, (
            "expected the all-connections-inactive error"
        )
        assert "STATUS INACTIVE" in text, "expected the STATUS INACTIVE remedy"

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: (get_model(name, t.workspace, perioddata), None),
        check=check,
        xfail=True,
    )
    test.run()


@pytest.mark.parametrize(
    "name, bottom, message",
    [
        ("maw_cs_botscrn", -15.0, "is above the screen bottom of active connection"),
        ("maw_cs_bothead", 10.0, "must be >= BOTTOM_ELEVATION"),
    ],
)
def test_bottom_errors(function_tmpdir, targets, name, bottom, message):
    # a BOTTOM above the screen bottom of an active connection, or above the
    # head in the well, is an error
    perioddata = {
        0: [[0, "rate", mawrate]],
        1: [[0, "bottom", bottom]],
    }

    def check(test):
        # error messages are wrapped across lines in the listing file
        text = " ".join(open(os.path.join(test.workspace, "mfsim.lst")).read().split())
        assert message in text, f"expected BOTTOM error: {message}"

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: (get_model(name, t.workspace, perioddata), None),
        check=check,
        xfail=True,
    )
    test.run()
