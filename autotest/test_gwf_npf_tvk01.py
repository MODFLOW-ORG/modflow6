import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk01"]
time_varying_k = [1.0, 10.0]


def build_models(idx, test):
    nlay, nrow, ncol = 3, 3, 3
    perlen = [100.0, 100.0]
    nper = len(perlen)
    nstp = nper * [1]
    tsmult = nper * [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    laytyp = 0
    botm = [0.0, -1.0, -2.0]
    strt = 1.0
    hk = 0.1

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
    gwf.name_file.save_flows = True

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
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
    )

    # tvk
    tvkspd = {}
    for kper in range(nper):
        hydraulic_conductivity = time_varying_k[kper]
        spd = []
        for k in range(nlay):
            for i in range(nrow):
                for j in range(ncol):
                    spd.append([(k, i, j), "K", hydraulic_conductivity])
        tvkspd[kper] = spd
    tvk = flopy.mf6.ModflowUtltvk(
        npf, print_input=True, perioddata=tvkspd, filename=tvk_filename
    )

    # chd files
    chdspd = []
    for k in range(nlay):
        for i in range(nrow):
            chdspd.append([(k, i, 0), top + 1])
            chdspd.append([(k, i, ncol - 1), top])
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
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    gwfname = "gwf_" + test.name

    # head
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    try:
        hobj = flopy.utils.HeadFile(fpth, precision="double")
        head = hobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # budget
    fpth = os.path.join(test.workspace, f"{gwfname}.cbc")
    try:
        bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        bud_allspd = bobj.get_data(text="CHD")
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.
    hk = time_varying_k[0]
    delc = 1.0
    delr = 1.0
    delz = 1.0
    dh = 1.0
    dl = 2 * delr

    for kper, bud in enumerate(bud_allspd):
        flow_rate_calc = time_varying_k[kper] * delc * delz * dh / dl
        print(f"Calculated q is {flow_rate_calc}")
        for node, node2, q in bud:
            print(node, node2, q, flow_rate_calc)
            assert np.isclose(flow_rate_calc, abs(q)), (
                f"Expected flow rate {flow_rate_calc} but found {q}"
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
