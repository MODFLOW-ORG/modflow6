"""SFR input error/validation regression tests. Each case is a minimal,
deliberately-invalid model that should fail to load, with a specific
error message expected in mfsim.lst:

  sfr-ifno-errors      duplicate + missing PACKAGEDATA reach
  sfr-cd-ifno-mm       out-of-range CONNECTIONDATA IFNO + missing reach
  sfr-conn-cycle       mutual-upstream connectivity contradiction
  sfr-obs-badreach     OBS entry referencing an out-of-range reach
  sfr-cprior-error     misspelled CPRIOR keyword in DIVERSIONS
"""

import flopy
import pytest
from framework import TestFramework

cases = [
    "sfr-ifno-errors",
    "sfr-cd-ifno-mm",
    "sfr-conn-cycle",
    "sfr-obs-badreach",
    "sfr-cprior-error",
]


def _build_ifno_errors(test, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 2), 4.0]]
    )
    # nreaches=3: reach 0 given twice (duplicate), reach 1 given once, reach
    # 2 never given (missing). All three reaches have no connections.
    packagedata = [
        [0, (0, 0, 0), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 0, 1.0, 0],
        [0, (0, 0, 0), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 0, 1.0, 0],
        [1, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 0, 1.0, 0],
    ]
    connectiondata = [[0], [1]]
    perioddata = {0: [[0, "status", "inactive"]]}
    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _build_cd_ifno_mismatch(test, name):
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
    # 3-reach linear chain; PACKAGEDATA declares valid IFNOs 0, 1, 2. The
    # middle row's IFNO is a typo (5, no matching PACKAGEDATA entry) instead
    # of 1.
    packagedata = [
        [0, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
        [1, (0, 0, 2), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 2, 1.0, 0],
        [2, (0, 0, 3), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
    ]
    connectiondata = [[0, -1], [5, 0, -2], [2, 1]]
    perioddata = {0: [[0, "inflow", 0.1]]}
    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _build_conn_cycle(test, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    # two reaches, each claiming the other as its upstream neighbor
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 1), 4.0]]
    )
    packagedata = [
        [0, (0, 0, 0), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
        [1, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
    ]
    connectiondata = [[0, 1], [1, 0]]
    perioddata = {0: [[0, "inflow", 0.1]]}
    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=2,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _build_obs_badreach(test, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    # valid 3-reach chain; OBS references a reach number far out of range
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 4), 4.0]]
    )
    packagedata = [
        [0, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
        [1, (0, 0, 2), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 2, 1.0, 0],
        [2, (0, 0, 3), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
    ]
    connectiondata = [[0, -1], [1, 0, -2], [2, 1]]
    perioddata = {0: [[0, "inflow", 0.1]]}
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    fname = f"{name}.sfr.obs"
    sfr_obs = {f"{fname}.csv": [("badreach", "stage", (99,))]}
    sfr.obs.initialize(filename=fname, continuous=sfr_obs)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


def _build_cprior_error(test, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 1), 4.0]]
    )
    # two reaches, reach 0 diverts to reach 1 with a misspelled CPRIOR
    packagedata = [
        [0, (0, 0, 0), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 1],
        [1, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 1, 1.0, 0],
    ]
    connectiondata = [[0, -1], [1, 0]]
    diversions = [[0, 0, 1, "FRACTON"]]  # typo of FRACTION
    perioddata = {0: [[0, "inflow", 0.1], [0, "diversion", 0, 0.5]]}
    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=2,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=diversions,
        perioddata=perioddata,
        pname="sfr-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{name}.cbc")
    return sim


_builders = [
    _build_ifno_errors,
    _build_cd_ifno_mismatch,
    _build_conn_cycle,
    _build_obs_badreach,
    _build_cprior_error,
]

_expected_lst_substrings = [
    [
        "PACKAGEDATA data for reach 1 specified 2 times",
        "PACKAGEDATA no data specified for reach 3",
    ],
    [
        "reach (6) must be greater than 0",
        "CONNECTIONDATA no data specified for reach 2",
    ],
    ["is not permitted."],
    ["less than or equal to", "specified value is"],
    ["Invalid cprior type"],
]


def build_models(idx, test):
    return _builders[idx](test, cases[idx])


def check_output(idx, test):
    with open(test.workspace / "mfsim.lst") as f:
        lst = f.read()
    for substr in _expected_lst_substrings[idx]:
        assert substr in lst, (
            f"expected {substr!r} in mfsim.lst for case {cases[idx]!r}; got:\n" + lst
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=True,
    )
    test.run()
