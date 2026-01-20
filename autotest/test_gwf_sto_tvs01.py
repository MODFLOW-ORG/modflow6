import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvs01"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 1
    perlen = [1.0, 1.0, 1.0, 1.0, 1.0]
    nper = len(perlen)
    nstp = nper * [1]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    laytyp = 1
    botm = [0.0]
    strt = [0.8]
    hk = 1.0
    ss = 1e-6
    sy = 0.01
    cellid1 = (0, 0, 0)

    nouter, ninner = 1000, 100
    hclose, rclose, relax = 1e-9, 1e-12, 0.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    transient = {}
    for i in range(nper):
        transient[i] = True

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = False
    gwf.name_file.newtonoptions = "NEWTON UNDER_RELAXATION"

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
        under_relaxation_kappa=0.07,
        under_relaxation_gamma=0.1,
        under_relaxation_momentum=0.0,
        backtracking_number=200,
        backtracking_tolerance=1.1,
        backtracking_reduction_factor=0.2,
        backtracking_residual_limit=1.0,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=laytyp, k=hk, k33=hk)

    # storage
    tvs_filename = f"{gwfname}.sto.tvs"
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        transient=transient,
    )

    # TVS
    tvsspd = {}

    # TVS SP1: No changes. Check h1 == 0.8.

    # TVS SP2: Decrease SY1. Check h1 == 3000.823.
    kper = 2
    spd = []
    spd.append([cellid1, "SY", 0.005])
    tvsspd[kper - 1] = spd

    # TVS SP3: Increase SS1. Check h1 == 300.5323.
    kper = 3
    spd = []
    spd.append([cellid1, "SS", 1e-5])
    tvsspd[kper - 1] = spd

    # TVS SP4: Increase SY1. Check h1 == 0.4 approx. (0.399976)
    kper = 4
    spd = []
    spd.append([cellid1, "SY", 0.02])
    tvsspd[kper - 1] = spd

    # TVS SP5: Revert SS1 and SY1. Check h1 == 0.8.
    kper = 5
    spd = []
    spd.append([cellid1, "SS", 1e-6])
    spd.append([cellid1, "SY", 0.01])
    tvsspd[kper - 1] = spd

    tvs = flopy.mf6.ModflowUtltvs(
        sto, print_input=True, perioddata=tvsspd, filename=tvs_filename
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    gwfname = "gwf_" + test.name

    # head
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        head = hobj.get_alldata()
    except:
        assert False, f'could not load data from "{fpth}"'

    # Check against manually calculated results
    expected_results = []
    expected_results.append(0.8)  # TVS SP1: No changes. Check initial solution.
    expected_results.append(3000.823)  # TVS SP2: Decrease SY1.
    expected_results.append(300.5323)  # TVS SP3: Increase SS1.
    expected_results.append(0.399976)  # TVS SP4: Increase SY1.
    # TVS SP5: Revert SS1 and SY1. Check that solution returns to original.
    expected_results.append(0.8)
    nper = len(expected_results)
    ex_lay = 1
    ex_row = 1
    ex_col = 1

    for kper, expected_result in enumerate(expected_results):
        h = head[kper, ex_lay - 1, ex_row - 1, ex_col - 1]
        print(kper, h, expected_result)
        assert np.isclose(h, expected_result), (
            f"Expected head {expected_result} in period {kper} but found {h}"
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
