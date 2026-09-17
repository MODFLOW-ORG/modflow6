"""SFR diversion (DIVERSIONRECORD) tests.

sfr_div            basic diversion fraction correctness
sfr-divaux-dedup   reissuing DIVERSION must not disturb AUXILIARY's
                    own TS tracking in the same keystring, and vice
                    versa
sfr-div-ts         a diversion (reach count > total diversions) goes
                    literal -> TS-linked -> an unreissued gap period
                    -> literal again, 2 timesteps/period throughout;
                    an unrelated reach's own diversion is untouched
                    the whole time as a cross-check
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ["sfr_div", "sfr-divaux-dedup", "sfr-div-ts"]
inflows = np.array([10, 0, 10, 0, 10])
diversion = np.array([0.5, 0.5, 0.5, 0.5, 0.0])


def _build_sfr_div(test):
    # static model data
    # temporal discretization
    nper = len(inflows)
    tdis_rc = [(1.0, 1, 1.0)] * nper

    # spatial discretization data
    nlay, nrow, ncol = 1, 1, 1
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = -10.0
    strt = 0.0

    # build MODFLOW 6 files
    name = cases[0]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    sim.simulation_data.verify_data = False

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="days",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(sim, print_option="ALL")

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=name + ".cbb", saverecord=[["BUDGET", "ALL"]]
    )

    # sfr file
    cellid = (0, 0, 0)
    nreaches = 3
    rlen = 10.0
    rwid = 10.0
    roughness = 0.001
    rbth = 1.0
    rhk = 0.0
    slope = 0.001

    sfrrch_prop = [cellid, rlen, rwid, slope, top, rbth, rhk, roughness]
    packagedata = [
        [0] + sfrrch_prop + [2, 1.0, 1],
        [1] + sfrrch_prop + [1, 0.0, 0],
        [2] + sfrrch_prop + [1, 1.0, 0],
    ]
    connectiondata = [
        [0, -1, -2],
        [1, 0],
        [2, 0],
    ]
    diversiondata = [
        [0, 0, 1, "FRACTION"],
    ]
    perioddata = {
        i: [[0, "inflow", qin], [0, "diversion", 0, qdiv]]
        for i, (qin, qdiv) in enumerate(zip(inflows, diversion))
    }

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        budgetcsv_filerecord=name + ".sfr.csv",
        budget_filerecord=name + ".sfr.cbb",
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=diversiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )

    return sim


def _check_sfr_div(test):
    # check flow for individual reach
    fname = os.path.join(test.workspace, f"{test.name}.sfr.cbb")
    with flopy.utils.CellBudgetFile(fname) as cbb:
        outflows = cbb.get_data(text="EXT-OUTFLOW")

    # check outflow for reach 2 and 3
    assert np.allclose([r.q[1] for r in outflows], -inflows * diversion), (
        "Incorrect outflow for diversion reach"
    )
    assert np.allclose([r.q[2] for r in outflows], -inflows * (1 - diversion)), (
        "Incorrect outflow for outlet reach"
    )

    # load SFR budget CSV and check overall budget
    with open(fname.replace(".cbb", ".csv")) as f:
        header = f.readline().strip().split(",")
        flux = np.loadtxt(f, delimiter=",")

    assert np.allclose(flux[:, header.index("EXT-OUTFLOW_IN")], 0), (
        "External flow IN larger than zero"
    )
    assert np.allclose(flux[:, header.index("PERCENT_DIFFERENCE")], 0), (
        "Large mass balance error in SFR"
    )


# ---------------------------------------------------------------------------
# sfr-divaux-dedup: cross-record independence (AUXILIARYRECORD/DIVERSIONRECORD
# share one keystring; reissuing one must not disturb the other)
# ---------------------------------------------------------------------------
_dedup_nper = 3
_dedup_inflow = 10.0
_dedup_div_ts_vals = [0.3] * 4
_dedup_conc_ts_vals = [50.0, 50.0, 100.0, 200.0]
_dedup_div_literal = 0.7
_dedup_conc_literal = 99.0
_dedup_div_expected = [0.3, _dedup_div_literal, _dedup_div_literal]
_dedup_conc_expected = [50.0, 100.0, _dedup_conc_literal]


def _build_divaux_dedup(ws, name):
    """AUXILIARY and DIVERSION each start TS-linked, then each is
    reissued as a literal on a different period than the other."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, nper=_dedup_nper, perioddata=[(1.0, 1, 1.0)] * _dedup_nper
    )
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=100.0, delc=100.0, top=0.0, botm=-10.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    cellid = (0, 0, 0)
    sfrrch = [cellid, 100.0, 5.0, 1e-3, 0.0, 1.0, 0.0, 0.03]
    packagedata = [
        [0] + sfrrch + [2, 1.0, 1, "conc_ts"],
        [1] + sfrrch + [1, 0.0, 0, 0.0],
        [2] + sfrrch + [1, 1.0, 0, 0.0],
    ]
    connectiondata = [[0, -1, -2], [1, 0], [2, 0]]
    diversiondata = [[0, 0, 1, "FRACTION"]]
    perioddata = {
        0: [
            [0, "inflow", _dedup_inflow],
            [0, "diversion", 0, "div_ts"],
            [0, "AUXILIARY", "conc", "conc_ts"],
        ],
        1: [[0, "inflow", _dedup_inflow], [0, "diversion", 0, _dedup_div_literal]],
        2: [
            [0, "inflow", _dedup_inflow],
            [0, "AUXILIARY", "conc", _dedup_conc_literal],
        ],
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_input=True,
        auxiliary=["conc"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=diversiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    sfr.ts.initialize(
        filename=f"{name}.sfr.ts",
        timeseries=list(
            zip([0.0, 1.0, 2.0, 3.0], _dedup_div_ts_vals, _dedup_conc_ts_vals)
        ),
        time_series_namerecord=["div_ts", "conc_ts"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim


def _check_divaux_dedup(test, name):
    fname = os.path.join(test.workspace, f"{name}.{paktest}.cbc")
    cbc = flopy.utils.CellBudgetFile(fname, precision="double")

    outflow = cbc.get_data(text="EXT-OUTFLOW")
    div_obs = np.array([-rec["q"][1] / _dedup_inflow for rec in outflow])
    assert np.allclose(div_obs, _dedup_div_expected), (
        f"DIVFLOW fraction expected {_dedup_div_expected}, got {div_obs}"
    )

    aux = cbc.get_data(text="AUXILIARY")
    auxcol = aux[0].dtype.names[-1]
    conc_obs = np.array([rec[auxcol][0] for rec in aux])
    assert np.allclose(conc_obs, _dedup_conc_expected), (
        f"AUXILIARY conc expected {_dedup_conc_expected}, got {conc_obs}"
    )


# ---------------------------------------------------------------------------
# sfr-div-ts: reach count > total diversions, exercising the offset table's
# out-of-range case; 2 timesteps/period probes for stale TS values
# ---------------------------------------------------------------------------
_divts_nper = 4
_divts_inflow = 10.0
_divts_r1_frac = 0.1
_divts_literal0 = 0.2
_divts_literal1 = 0.9
_divts_ts_times = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
_divts_ts_vals = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
# (totim, reach 1 expected fraction, reach 4 expected fraction)
_divts_checks = [
    (2.0, _divts_r1_frac, _divts_literal0),
    (3.0, _divts_r1_frac, 0.35),
    (4.0, _divts_r1_frac, 0.40),
    (5.0, _divts_r1_frac, 0.45),
    (6.0, _divts_r1_frac, 0.50),
    (8.0, _divts_r1_frac, _divts_literal1),
]


def _build_div_ts(ws, name):
    """Reach 1 (control) and reach 4 each get their own two termini, so
    each's diverted fraction is independently observable."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, nper=_divts_nper, perioddata=[(2.0, 2, 1.0)] * _divts_nper
    )
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=100.0, delc=100.0, top=0.0, botm=-10.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    cellid = (0, 0, 0)
    src = [cellid, 100.0, 5.0, 1e-3, 0.0, 1.0, 0.0, 0.03]
    term = [cellid, 100.0, 5.0, 1e-3, 0.0, 1.0, 0.0, 0.03]
    packagedata = [
        [0] + src + [2, 1.0, 1],  # reach1: source, 1 diversion
        [1] + term + [1, 1.0, 0],  # reach2: reach1's main terminus
        [2] + term + [1, 0.0, 0],  # reach3: reach1's diversion terminus
        [3] + src + [2, 1.0, 1],  # reach4: source, 1 diversion
        [4] + term + [1, 1.0, 0],  # reach5: reach4's main terminus
        [5] + term + [1, 0.0, 0],  # reach6: reach4's diversion terminus
    ]
    connectiondata = [
        [0, -1, -2],
        [1, 0],
        [2, 0],
        [3, -4, -5],
        [4, 3],
        [5, 3],
    ]
    diversiondata = [
        [0, 0, 2, "FRACTION"],
        [3, 0, 5, "FRACTION"],
    ]  # reach1->3, reach4->6
    perioddata = {
        0: [
            [0, "inflow", _divts_inflow],
            [3, "inflow", _divts_inflow],
            [0, "diversion", 0, _divts_r1_frac],
            [3, "diversion", 0, _divts_literal0],
        ],
        1: [[3, "diversion", 0, "div_ts"]],
        2: [],
        3: [[3, "diversion", 0, _divts_literal1]],
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_input=True,
        budget_filerecord=f"{name}.{paktest}.cbc",
        nreaches=6,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=diversiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    sfr.ts.initialize(
        filename=f"{name}.sfr.ts",
        timeseries=list(zip(_divts_ts_times, _divts_ts_vals)),
        time_series_namerecord=["div_ts"],
        interpolation_methodrecord=["linearend"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim


def _check_div_ts(test, name):
    # diversion destinations are reach 3 (reach 1's) and reach 6 (reach 4's)
    fname = os.path.join(test.workspace, f"{name}.{paktest}.cbc")
    cbc = flopy.utils.CellBudgetFile(fname, precision="double")
    for totim, r1_expected, r4_expected in _divts_checks:
        rec = cbc.get_data(totim=totim, text="EXT-OUTFLOW")[0]
        r1_obs = -rec["q"][rec["node"] == 3][0] / _divts_inflow
        r4_obs = -rec["q"][rec["node"] == 6][0] / _divts_inflow
        assert np.isclose(r1_obs, r1_expected), (
            f"t={totim}: reach 1 fraction expected {r1_expected}, got {r1_obs}"
        )
        assert np.isclose(r4_obs, r4_expected), (
            f"t={totim}: reach 4 fraction expected {r4_expected}, got {r4_obs}"
        )


def build_models(idx, test):
    if idx == 0:
        return _build_sfr_div(test), None
    elif idx == 1:
        return _build_divaux_dedup(test.workspace, cases[idx]), None
    return _build_div_ts(test.workspace, cases[idx]), None


def check_output(idx, test):
    if idx == 0:
        _check_sfr_div(test)
    elif idx == 1:
        _check_divaux_dedup(test, cases[idx])
    else:
        _check_div_ts(test, cases[idx])


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
