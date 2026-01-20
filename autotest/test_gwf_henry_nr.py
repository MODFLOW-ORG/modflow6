"""
The Henry, Newton-Raphson problem described by Langevin et al (2020)
with a 20x40 grid instead of the 40x80 grid described in the paper.
There is freshwater inflow on the left. A sloping sea boundary on the
right moves up and down according to a simple sine function. GHBs and
DRNs alternate and move up and down along the boundary to represent
the effects of tides on the aquifer.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_henrynr01"]

# global model variables
nlay = 20
ncol = 40
fx = 0.5
fz = 0.5
lx = 2
lz = 1
sealevel = 0.85
amplitude = 0.14
frequency = 4
wellfact = 0.25


def get_idomain(nlay, nrow, ncol, lx, lz, fx, fz):
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    x1 = fx * lx
    y1 = lz
    x2 = lx
    y2 = fz * lz
    slope = (y2 - y1) / (x2 - x1)
    b = y1 - slope * x1

    delr = lx / ncol
    delv = lz / nlay
    xcenters = np.linspace(delr / 2, lx - delr / 2, ncol)
    zcenters = np.linspace(lz - delv / 2, delv / 2, nlay)

    for k in range(nlay):
        zc = zcenters[k]
        for j in range(ncol):
            xc = xcenters[j]
            zedge = slope * xc + b
            if zc > zedge:
                idomain[k, 0, j] = 0

    kidm0, iidmn0, jidmn0 = np.where(idomain == 0)
    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k, 0, j - 1] == 1:
            idomain[k, 0, j - 1] = 2

    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k + 1, 0, j] == 1:
            idomain[k + 1, 0, j] = 3

    return idomain


def sinfunc(a, b, c, d, x):
    return a * np.sin(b * (x - c)) + d


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]

    nrow = 1
    delr = lx / ncol
    delc = 1.0
    top = lz
    delz = lz / nlay
    botm = list(top - np.arange(delz, nlay * delz + delz, delz))

    perlen = [0.25] + 1000 * [0.001]
    nper = len(perlen)
    nstp = [250] + 1000 * [1]
    tsmult = 1.0
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult))

    nouter, ninner = 200, 50
    hclose, rclose, relax = 1e-9, 1e-6, 0.97

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    sim.name_file.continue_ = False

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        csv_outer_output_filerecord=name + ".ims.csv",
        outer_dvclose=hclose * 10.0,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        no_ptcrecord=True,
    )

    # --------------------  FLOW --------------------

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON", save_flows=True
    )

    idomain = get_idomain(nlay, nrow, ncol, lx, lz, fx=fx, fz=fz)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=lz)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_specific_discharge=True,
        icelltype=1,
        k=864.0,
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf, sy=0.35, iconvert=1, steady_state=[False], transient=[True]
    )

    # drn and ghb
    kidx, iidx, jidx = np.where(idomain > 1)
    xcellcenters = gwf.modelgrid.xcellcenters
    zcellcenters = gwf.modelgrid.zcellcenters
    botm = dis.botm.get_data()
    dt = 0.001
    times = np.arange(dt, 1.0 + dt, dt)
    sealevelts = [sealevel] + list(
        sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, times)
    )
    ghbspd = {}
    drnspd = {}
    for kper in range(nper):
        if kper == 0:
            sl = sealevel
        else:
            sl = sealevelts[kper]
        ghblist = []
        drnlist = []
        for k, i, j in zip(kidx, iidx, jidx):
            zcell = zcellcenters[k, i, j]
            cond = 864.0 * (delz * delc) / (0.5 * delr)
            if zcell > sl:
                drnlist.append([(k, i, j), zcell, 864.0, 0.0])
            else:
                ghblist.append([(k, i, j), sl, 864.0, 35.0, 1024.5])
        if len(ghblist) > 0:
            ghbspd[kper] = ghblist
        if len(drnlist) > 0:
            drnspd[kper] = drnlist

    # drn
    drn1 = flopy.mf6.ModflowGwfdrn(
        gwf,
        stress_period_data=drnspd,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="DRN-1",
        auxiliary="CONCENTRATION",
    )

    # ghb
    ghb1 = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="GHB-1",
        auxiliary=["CONCENTRATION", "DENSITY"],
    )

    wellist1 = []
    qwell = 5.7024 * wellfact
    qwell = qwell / nlay
    for k in range(nlay):
        wellist1.append([(k, 0, 0), qwell, 0.0])
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wellist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
        auxiliary="CONCENTRATION",
        filename=f"{name}.wel",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    return sim, None


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    name = "gwf-henry-nr"
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        compare="mf6_regression",
        verbose=False,
    )
    test.run()
