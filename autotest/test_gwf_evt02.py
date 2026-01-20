import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["evt02"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 3
    chdheads = list(np.linspace(1, 100))
    nper = len(chdheads)
    perlen = nper * [1.0]
    nstp = nper * [1]
    tsmult = nper * [1.0]

    delr = delc = 1.0
    strt = chdheads[0]

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
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
        top=100.0,
        botm=0.0,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=1.0)

    # chd files
    chdspd = {}
    for kper, chdval in enumerate(chdheads):
        chdspd[kper] = [[(0, 0, 0), chdval], [(0, 0, ncol - 1), chdval]]
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

    nseg = 1
    surf_rate_specified = True
    evtspd = [[(0, 0, 1), 95.0, 0.001, 90.0, 0.1]]

    # nseg = 4
    # surf_rate_specified = False
    # evtspd = [((0, 0, 1), 95., 0.001, 90., 0.25, 0.5, 0.75, 1., 0., 1.)]

    # nseg = 1
    # surf_rate_specified = False
    # evtspd = [[(0, 0, 1), 95., 0.001, 90.]]

    evt = flopy.mf6.ModflowGwfevt(
        gwf,
        print_flows=True,
        surf_rate_specified=surf_rate_specified,
        maxbound=1,
        nseg=nseg,
        stress_period_data=evtspd,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def etfunc(h, qmax, surf, exdp, petm, pxdp, petm0=1.0):
    nseg = len(petm) + 1
    d = surf - h
    if h >= surf:
        hcof = 0.0
        rhs = qmax * petm0
    elif d >= exdp:
        hcof = 0.0
        rhs = 0.0
    else:
        if nseg > 1:
            pxdp1 = 0.0
            petm1 = petm0
            for iseg in range(nseg):
                if iseg < nseg - 1:
                    pxdp2 = pxdp[iseg]
                    petm2 = petm[iseg]
                else:
                    pxdp2 = 1.0
                    petm2 = 0.0
                if d <= pxdp2 * exdp:
                    break
                pxdp1 = pxdp2
                petm1 = petm2
            hcof = -(petm1 - petm2) * qmax / ((pxdp2 - pxdp1) * exdp)
            rhs = hcof * (surf - pxdp1 * exdp) + petm1 * qmax
        else:
            hcof = -qmax / exdp
            rhs = qmax - qmax * surf / exdp
    q = h * hcof - rhs
    return q, hcof, rhs


def check_output(idx, test):
    # The nature of the bug is that the model crashes with nseg=1
    fpth = os.path.join(test.workspace, "evt02.cbc")
    assert os.path.isfile(fpth), "model did not run with nseg=1 in EVT input"


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
