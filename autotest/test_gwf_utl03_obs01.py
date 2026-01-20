import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["utl03_obs"]

# temporal discretization
nper = 2
tdis_rc = [(1.0, 1, 1.0), (1.0, 1, 1.0)]

# spatial discretization data
nlay, nrow, ncol = 1, 21, 21
top = 200.0
botm = 0.0

# hydraulic properties
laytyp = 1
hk = 1.0
delr = delc = 100.0
strt = botm + 20.0

# constant head data
chdlocl = [(0, i, 0) for i in range(nrow)]
chdlocr = [(0, i, ncol - 1) for i in range(nrow)]
chdl = 100.0
chdr = 75.0

c60 = []
c61 = []
for loc in chdlocl:
    c60.append([loc, chdl])
    c61.append([loc, chdl + 10])
for loc in chdlocr:
    c60.append([loc, chdr])
    c61.append([loc, chdr - 10])
cd6 = {0: c60, 1: c61}

# gwf obs
obs_data0 = [(f"h{i + 1:04d}", "HEAD", (0, 10, 10)) for i in range(1000)]
obs_data1 = [(f"h{i + 1001:04d}", "HEAD", (0, 1, 1)) for i in range(737)]

# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 0.01, 1.0


def build_mf6(idx, ws, binaryobs=True):
    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", sim_ws=ws)
    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        model_nam_file=f"{name}.nam",
        save_flows=True,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="none",
        backtracking_number=0,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        number_orthogonalizations=7,
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk)

    # gwf head observations
    if binaryobs:
        obs_recarray = {"head0.obs.bsv": obs_data0, "head1.obs.bsv": obs_data1}
    else:
        obs_recarray = {"head0.obs.csv": obs_data0, "head1.obs.csv": obs_data1}

    o = flopy.mf6.ModflowUtlobs(
        gwf,
        pname="head_obs",
        filename=f"{name}.obs",
        digits=10,
        print_input=True,
        continuous=obs_recarray,
    )

    # chd files
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=cd6)

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim


def build_model(idx, dir):
    ws = dir
    # build mf6 with ascii observation output
    sim = build_mf6(idx, ws, binaryobs=False)

    # build mf6 with binary observation output
    wsc = os.path.join(ws, "mf6")
    mc = build_mf6(idx, wsc, binaryobs=True)

    sim.write_simulation()
    mc.write_simulation()
    hack_binary_obs(idx, dir)

    return sim, mc


def build_models(idx, test):
    sim, mc = build_model(idx, test.workspace)
    sim.write_simulation()
    mc.write_simulation()
    hack_binary_obs(idx, test.workspace)
    return sim, mc


def hack_binary_obs(idx, dir):
    name = cases[idx]
    ws = dir
    wsc = os.path.join(ws, "mf6")
    fname = name + ".obs"
    fpth = os.path.join(wsc, fname)
    with open(fpth, "r") as f:
        lines = f.readlines()
    with open(fpth, "w") as f:
        for line in lines:
            line = line.rstrip()
            if "BEGIN continuous  FILEOUT" in line:
                line += "  BINARY"
            f.write(f"{line}\n")
        f.close()


def check_output(idx, test):
    # get results from the observation files
    pth = test.workspace
    files = [fn for fn in os.listdir(pth) if ".csv" in fn]
    for file in files:
        pth0 = os.path.join(pth, file)
        pth1 = os.path.join(pth, "mf6", file.replace(".csv", ".bsv"))
        d0 = flopy.utils.Mf6Obs(pth0, isBinary=False).get_data()
        d1 = flopy.utils.Mf6Obs(pth1, isBinary=True).get_data()
        names0 = d0.dtype.names
        names1 = d1.dtype.names
        msg = (
            f"The number of columns ({len(names0)}) "
            + f"in {pth0} "
            + "is not equal to "
            + f"the number of columns ({len(names1)}) "
            + f"in {pth1}."
        )
        assert len(names0) == len(names1), msg
        msg = (
            f"The number of rows ({d0.shape[0]}) "
            + f"in {pth0} "
            + "is not equal to "
            + f"the number of rows ({d1.shape[0]}) "
            + f"in {pth1}."
        )
        assert d0.shape[0] == d1.shape[0], msg
        for name in names0:
            msg = (
                f"The values for column '{name}' " + "are not within 1e-5 of each other"
            )
            assert np.allclose(d0[name], d1[name], rtol=1e-5), msg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        overwrite=False,
    )
    test.run()
