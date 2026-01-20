import os

import flopy
import pytest
from framework import TestFramework

cases = ["gwf_noptc01", "gwf_noptc02", "gwf_noptc03"]
no_ptcrecords = ["FIRST", "ALL", None]
htol = [None for _ in range(len(cases))]

# static model data
# temporal discretization
nper = 1
tdis_rc = [(1.0, 1, 1.0)]

# spatial discretization data
nlay, nrow, ncol = 1, 1, 100
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 50.0, 1.0
top = 25.0
botm = 0.0
strt = 0.0

# hydraulic properties
hk = 50.0

# all cells are active and layer 1 is convertible
ib = 1

# solver options
nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0
newtonoptions = "NEWTON"
imsla = "BICGSTAB"

# chd data
c6 = []
ccol = [0, ncol - 1]
hc = [20.0, 11.0]
for j, h in zip(ccol, hc):
    c6.append([(0, 0, j), h])
cd6 = {0: c6}
maxchd = len(cd6[0])

# recharge data
rech = {0: 0.001}


def get_model(idx, dir, no_ptcrecord):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions=newtonoptions, save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        no_ptcrecord=no_ptcrecord,
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=imsla,
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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=1, k=hk)

    # recharge
    rch = flopy.mf6.ModflowGwfrcha(gwf, readasarrays=True, recharge=rech)

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
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

    return sim


# water table recharge problem
def build_models(idx, test):
    sim = get_model(idx, test.workspace, no_ptcrecords[idx])

    # build MODFLOW-6 without no_ptc option
    pth = os.path.join(test.workspace, "mf6")
    mc = get_model(idx, pth, None)

    return sim, mc


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
    )
    test.run()
