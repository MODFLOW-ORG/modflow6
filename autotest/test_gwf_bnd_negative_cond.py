"""
This test build models with DRN, RIV, or GHB packages with
negative conductance multipliers or negative conductance
values. All of these models should terminate with an error
message.

The test evaluates the error messages in the mfsim.lst file
for the correct error message.
"""

import flopy
import pytest
from framework import TestFramework

paktest = "bnd"
cases = [
    "drn-cond",
    "drn-mult",
    "riv-cond",
    "riv-mult",
    "ghb-cond",
    "ghb-mult",
]


def build_models(idx, test):
    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = []
    for _ in range(nper):
        tdis_rc.append((1.0, 1, 1.0))

    # spatial discretization data
    nlay, nrow, ncol = 1, 1, 2
    delr, delc = 100.0, 100.0
    top = 0.0
    botm = -10
    strt = 0.0

    # calculate hk
    hk = 1.0e-4

    # solver options
    hclose, rclose = 1e-9, 1e-3

    # build MODFLOW 6 files
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=hclose,
        under_relaxation="NONE",
        inner_dvclose=hclose,
        rcloserecord=rclose,
    )

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
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk)

    # chd files
    # chd data
    spd = [
        [(0, 0, 0), 1.0],
    ]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd, pname="chd-1")

    bnd_loc = (0, 0, 1)
    cond = 1.0
    mult = 1.0
    if name.endswith("mult"):
        mult *= -1.0
    else:
        cond *= -1.0
    if name.startswith("drn"):
        drn = flopy.mf6.ModflowGwfdrn(
            gwf,
            auxiliary=["mult"],
            auxmultname="mult",
            stress_period_data=[(bnd_loc, top, cond, mult)],
        )
    elif name.startswith("riv"):
        riv = flopy.mf6.ModflowGwfriv(
            gwf,
            auxiliary=["mult"],
            auxmultname="mult",
            stress_period_data=[(bnd_loc, 1.0, cond, top, mult)],
        )
    elif name.startswith("ghb"):
        ghb = flopy.mf6.ModflowGwfghb(
            gwf,
            auxiliary=["mult"],
            auxmultname="mult",
            stress_period_data=[(bnd_loc, top, cond, mult)],
        )

    return sim


def check_output(idx, test):
    print("Running error check")
    name = cases[idx]
    if name.startswith("drn"):
        pak = "DRN"
    elif name.startswith("riv"):
        pak = "RIV"
    elif name.startswith("ghb"):
        pak = "GHB"
    if name.endswith("mult"):
        tag = (
            f"1. {pak} BOUNDARY (1) CONDUCTANCE "
            + "MULTIPLIER ( -1.00    ) IS LESS THAN ZERO"
        )
    else:
        tag = f"1. {pak} BOUNDARY (1) CONDUCTANCE " + "( -1.00    ) IS LESS THAN ZERO"
    with open(test.workspace / "mfsim.lst", "r") as f:
        lines = f.readlines()
        error_count = 0
        for line in lines:
            if tag in line:
                error_count += 1

        # ensure that error msg is in mfsim.lst file
        assert error_count == 1, (
            "error count = " + str(error_count) + "but should equal 1"
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
