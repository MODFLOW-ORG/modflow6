import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk03"]


def build_models(idx, test):
    nlay, nrow, ncol = 3, 1, 1
    perlen = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    nper = len(perlen)
    nstp = nper * [1]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    top = 0.0
    laytyp = 0
    botm = [-1.0, -2.0, -3.0]
    strt = [1.0, 1.0, 0.0]
    hk = 1.0
    cellid1 = (0, 0, 0)
    cellid3 = (2, 0, 0)

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

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

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
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
    tvk_filename = f"{gwfname}.npf.tvk"
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=laytyp, k=hk, k33=hk)

    # TVK
    tvkspd = {}

    # TVK SP1: No changes. Check initial solution, h2 == 0.5.

    # TVK SP2: Increase K3. Check h2 == 0.4.
    kper = 2
    spd = []
    spd.append([cellid3, "K33", 3.0])
    tvkspd[kper - 1] = spd

    # TVK SP3: Decrease K3. Check h2 == 0.75.
    kper = 3
    spd = []
    spd.append([cellid3, "K33", 0.2])
    tvkspd[kper - 1] = spd

    # TVK SP4: Revert K3 and increase K1. Check h2 == 0.6.
    kper = 4
    spd = []
    spd.append([cellid3, "K33", 1.0])
    spd.append([cellid1, "K33", 3.0])
    tvkspd[kper - 1] = spd

    # TVK SP5: Decrease K1. Check h2 == 0.25.
    kper = 5
    spd = []
    spd.append([cellid1, "K33", 0.2])
    tvkspd[kper - 1] = spd

    # TVK SP6: No changes. Check that solution remains as per SP5, h2 == 0.25.

    # TVK SP7: Revert K1. Check that solution returns to original, h2 == 0.5.
    kper = 7
    spd = []
    spd.append([cellid1, "K33", 1.0])
    tvkspd[kper - 1] = spd

    tvk = flopy.mf6.ModflowUtltvk(
        npf, print_input=True, perioddata=tvkspd, filename=tvk_filename
    )

    # chd files
    chdspd = []
    chdspd.append([cellid1, 1.0])
    chdspd.append([cellid3, 0.0])
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        save_flows=False,
        print_flows=True,
        pname="CHD-1",
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
    expected_results.append(0.5)  # TVK SP1: No changes. Check initial solution.
    expected_results.append(0.4)  # TVK SP2: Increase K3.
    expected_results.append(0.75)  # TVK SP3: Decrease K3.
    expected_results.append(0.6)  # TVK SP4: Revert K3 and increase K1.
    expected_results.append(0.25)  # TVK SP5: Decrease K1.
    # TVK SP6: No changes. Check that solution remains as per SP5.
    expected_results.append(0.25)
    # TVK SP7: Revert K1. Check that solution returns to original.
    expected_results.append(0.5)
    nper = len(expected_results)
    ex_lay = 2
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
