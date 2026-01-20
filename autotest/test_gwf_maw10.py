"""
Test to confirm that the sum of rate-actual and maw-reduction observations
is equal to the specified MAW extraction or injection pumping rate when
reported using the MAW_FLOW_REDUCE_CSV option. Injection and extraction
are tested for both MAW HEAD_LIMIT and RATE_SCALING options to reduce MAW
flow rates from the specified input (requested) values.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["maw10a", "maw10b", "maw10c", "maw10d"]
mawsetting_a = {
    0: [
        [0, "rate", -2000.0],
        [0, "head_limit", -0.4],
    ],
    1: [
        [0, "status", "inactive"],
    ],
}
mawsetting_b = {
    0: [
        [0, "rate", -2000.0],
        [0, "rate_scaling", -1.4, 1.0],
    ],
    1: [
        [0, "status", "inactive"],
    ],
}
mawsetting_c = {
    0: [
        [0, "rate", 2000.0],
        [0, "head_limit", 0.4],
    ],
    1: [
        [0, "status", "inactive"],
    ],
}
mawsetting_d = {
    0: [
        [0, "rate", 2000.0],
        [0, "rate_scaling", 0.0, 1.0],
    ],
    1: [
        [0, "status", "inactive"],
    ],
}

# simple single stress period without going inactive in 2nd stress period
# mawsetting_a = [(0, "rate", -2000.0), (0, "head_limit", -0.4)]
# mawsetting_b = [
#     [(0, "rate", -2000.0), (0, "rate_scaling", -1.4, 1.0)],
#     ["status", "inactive"],
# ]
# mawsetting_c = [[(0, "rate", 2000.0), (0, "head_limit", 0.4)], ["status", "inactive"]]
# mawsetting_d = [
#     [(0, "rate", 2000.0), (0, "rate_scaling", 0.0, 1.0)],
#     ["status", "inactive"],
# ]
mawsettings = [mawsetting_a, mawsetting_b, mawsetting_c, mawsetting_d]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 101, 101
    nper = 2
    perlen = [500.0, 500.0]
    nstp = [50, 50]
    tsmult = [1.2, 1.2]
    delr = delc = 142.0
    top = 0.0
    botm = [-1000.0]
    strt = 0.0
    hk = 10.0

    nouter, ninner = 100, 100
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=0,
        ss=1.0e-5,
        sy=0.1,
        steady_state={0: False},
        transient={0: True},
        filename=f"{name}.sto",
    )

    # MAW
    opth = f"{name}.maw.obs"
    wellbottom = -1000
    wellrecarray = [[0, 0.15, wellbottom, 0.0, "THIEM", 1]]
    wellconnectionsrecarray = [[0, 0, (0, 50, 50), 0.0, wellbottom, 0.0, 0.0]]
    wellperiodrecarray = mawsettings[idx]
    mawo_dict = {}
    mawo_dict[f"{name}.maw.obs.csv"] = [
        ("m1head", "head", (0,)),
        ("m1rate", "rate", (0,)),
    ]  # is this index one-based? Not if in a tuple
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations=mawo_dict,
        packagedata=wellrecarray,
        connectiondata=wellconnectionsrecarray,
        perioddata=wellperiodrecarray,
        mfrcsv_filerecord=f"{name}.maw-reduction.csv",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    # head observations
    obs_data0 = [("head_well_cell", "HEAD", (0, 0, 0))]
    obs_recarray = {f"{name}.obs.csv": obs_data0}
    obs = flopy.mf6.ModflowUtlobs(
        gwf,
        pname="head_obs",
        filename=f"{name}.obs",
        digits=15,
        print_input=True,
        continuous=obs_recarray,
    )

    return sim, None


# 2 tests to perform here:
# 1. within the .maw-reduction.csv file, do values of actual + reduction = requested?
# 2. do the values in .maw-reduction.csv file match with the .maw.obs.csv file at
#    each time (and all are reduction times present in the obs file)?
def check_output(idx, test):
    # MODFLOW 6 maw results
    name = cases[idx]
    fpthobs = os.path.join(test.workspace, f"{name}.maw.obs.csv")
    fpthmfr = os.path.join(test.workspace, f"{name}.maw-reduction.csv")
    try:
        tcobs = np.genfromtxt(fpthobs, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpthobs}"'
    try:
        tcmfr = np.genfromtxt(fpthmfr, names=True, delimiter=",", deletechars="")
    except:
        assert False, f'could not load data from "{fpthmfr}"'

    # test 1: does rate-requested = rate-actual + maw-reduction at each time
    # in the .maw.reduced.csv?
    for rowmfr in tcmfr:
        v1 = rowmfr["rate-requested"]
        v2 = rowmfr["rate-actual"] + rowmfr["maw-reduction"]
        errmsg = (
            "MAW flow reduction output: "
            "requested rate must equal actual rate plus reduced rate.\n"
            f"{v1} /= {v2}"
        )
        assert np.allclose(v1, v2), errmsg

    # test 2: do values of rate-actual in .maw.reduced.csv equal those
    # flow values reported in .maw.obs.csv? (and are there matching times?)
    for rowmfr in tcmfr:
        timevalmfr = rowmfr["time"]
        actvalmfr = rowmfr["rate-actual"]
        msgtime = (
            "There should be a matching time in the maw.obs.csv file for each "
            "time in the maw.reduce.csv file, but no match was found for "
            f"time = {timevalmfr} in the maw.obs.csv file"
        )
        blnFoundTimeMatch = False
        for rowobs in tcobs:
            timevalobs = rowobs["time"]
            if np.isclose(timevalmfr, timevalobs):
                blnFoundTimeMatch = True
                actvalobs = rowobs["M1RATE"]
                msgval = (
                    f"The maw.obs.csv file rate-actual value of {actvalobs} "
                    "should have matched the maw.reduce.csv file rate-actual "
                    f"value of {actvalmfr} at time {timevalobs}"
                )
                break
        assert blnFoundTimeMatch, msgtime
        assert np.isclose(actvalmfr, actvalobs), msgval


@pytest.mark.slow
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
