import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk02", "tvk02_negk", "tvk02_disv"]


def build_base(test, name, disv=False):
    """Build the standard 1x1x3 TVK model, optionally using a DISV grid."""
    nlay, nrow, ncol = 1, 1, 3
    perlen = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    nper = len(perlen)
    nstp = nper * [1]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    laytyp = 0
    botm = [-1.0]
    strt = [0.0, 1.0, 1.0]
    hk = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = [(perlen[i], nstp[i], tsmult[i]) for i in range(nper)]

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

    if disv:
        # 3-cell DISV grid equivalent to 1x1x3 DIS.
        # MODFLOW 6 requires clockwise vertex ordering viewed from above.
        # Vertices:
        #   4----5----6----7  (y=1)
        #   |    |    |    |
        #   0----1----2----3  (y=0)
        vertices = [
            (0, 0.0, 0.0),
            (1, 1.0, 0.0),
            (2, 2.0, 0.0),
            (3, 3.0, 0.0),
            (4, 0.0, 1.0),
            (5, 1.0, 1.0),
            (6, 2.0, 1.0),
            (7, 3.0, 1.0),
        ]
        # Clockwise order: top-left → top-right → bottom-right → bottom-left
        cell2d = [
            (0, 0.5, 0.5, 4, 4, 5, 1, 0),
            (1, 1.5, 0.5, 4, 5, 6, 2, 1),
            (2, 2.5, 0.5, 4, 6, 7, 3, 2),
        ]
        dis = flopy.mf6.ModflowGwfdisv(
            gwf,
            nlay=nlay,
            ncpl=3,
            top=top,
            botm=botm,
            nvert=8,
            vertices=vertices,
            cell2d=cell2d,
            filename=f"{gwfname}.disv",
        )
        cellid1 = (0, 0)
        cellid3 = (0, 2)
        strt_disv = [0.0, 1.0, 1.0]
        # initial conditions
        ic = flopy.mf6.ModflowGwfic(gwf, strt=strt_disv, filename=f"{gwfname}.ic")
    else:
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
        cellid1 = (0, 0, 0)
        cellid3 = (0, 0, 2)
        # initial conditions
        ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    tvk_filename = f"{gwfname}.npf.tvk"
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=laytyp, k=hk, k33=hk)

    # TVK period data
    tvkspd = {}

    # TVK SP1: No changes. Check initial solution, h2 == 0.5.

    # TVK SP2: Increase K1. Check h2 == 0.4.
    tvkspd[1] = [[cellid1, "K", 3.0]]

    # TVK SP3: Decrease K1. Check h2 == 0.75.
    tvkspd[2] = [[cellid1, "K", 0.2]]

    # TVK SP4: Revert K1 and increase K3. Check h2 == 0.6.
    tvkspd[3] = [[cellid1, "K", 1.0], [cellid3, "K", 3.0]]

    # TVK SP5: Decrease K3. Check h2 == 0.25.
    tvkspd[4] = [[cellid3, "K", 0.2]]

    # TVK SP6: No changes. Check that solution remains as per SP5, h2 == 0.25.

    # TVK SP7: Revert K3. Check that solution returns to original, h2 == 0.5.
    tvkspd[6] = [[cellid3, "K", 1.0]]

    tvk = flopy.mf6.ModflowUtltvk(
        npf, print_input=True, perioddata=tvkspd, filename=tvk_filename
    )

    # chd
    chdspd = [[cellid1, 0.0], [cellid3, 1.0]]
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


def build_negk(test, name):
    """Build a model where TVK sets K to a negative value.
    Expected to fail: tvk_validate_change stores an error for K <= 0."""
    nlay, nrow, ncol = 1, 1, 3
    nper = 2
    tdis_rc = [(1.0, 1, 1.0)] * nper
    hk = 1.0

    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = False

    ims = flopy.mf6.ModflowIms(
        sim,
        outer_maximum=100,
        inner_maximum=300,
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
        rcloserecord=1e-6,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=1.0,
        botm=[-1.0],
        filename=f"{gwfname}.dis",
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5, filename=f"{gwfname}.ic")

    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk)

    # TVK sets K to -1.0 for cell (0,0,0) in SP2 — validate_change should error
    tvkspd = {1: [[(0, 0, 0), "K", -1.0]]}
    flopy.mf6.ModflowUtltvk(
        npf,
        print_input=True,
        perioddata=tvkspd,
        filename=f"{gwfname}.npf.tvk",
    )

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 0.0], [(0, 0, 2), 1.0]],
        pname="CHD-1",
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST")],
    )

    return sim, None


def build_models(idx, test):
    name = cases[idx]
    if name == "tvk02_negk":
        return build_negk(test, name)
    elif name == "tvk02_disv":
        return build_base(test, name, disv=True)
    else:
        return build_base(test, name, disv=False)


def check_output(idx, test):
    name = cases[idx]

    if name == "tvk02_negk":
        # Simulation should have failed; check that the K <= 0 error is in mfsim.lst
        with open(test.workspace / "mfsim.lst", "r") as f:
            lines = f.readlines()
        error_count = sum(1 for line in lines if "is <= 0" in line)
        assert error_count >= 1, (
            f"Expected 'is <= 0' validation error but found {error_count} occurrences"
        )
        return

    gwfname = "gwf_" + name
    # head
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        head = hobj.get_alldata()
    except Exception:
        assert False, f'could not load data from "{fpth}"'

    expected_results = [
        0.5,  # SP1: No changes.
        0.4,  # SP2: Increase K1.
        0.75,  # SP3: Decrease K1.
        0.6,  # SP4: Revert K1 and increase K3.
        0.25,  # SP5: Decrease K3.
        0.25,  # SP6: No changes (retains SP5 state).
        0.5,  # SP7: Revert K3.
    ]

    for kper, expected_result in enumerate(expected_results):
        # Middle cell: layer 1, row 1, col 2 (DIS) or layer 0, cell2d 1 (DISV).
        # Both yield head[kper, 0, 0, 1] in the binary output.
        h = head[kper, 0, 0, 1]
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
        xfail=(name == "tvk02_negk"),
    )
    test.run()
