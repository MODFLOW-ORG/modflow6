import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["nur-pr2606"]
nper = 1
nlay = 2
nrow = 1
ncol = 11
xlen = 1000.0
ylen = 1000.0
top = 25.0
k11 = 1.0
H1 = 7.5
H2 = 7.5
delr = xlen / float(ncol)
delc = ylen / float(nrow)
extents = (0, xlen, 0, ylen)
shape2d = (nrow, ncol)
shape3d = (nlay, nrow, ncol)
nouter = 75
ninner = 100
hclose = 1e-9
hclose_outer = hclose * 10.0
rclose = 1e-3
botm = [0, -25]
chd_spd = [[0, i, 0, H1] for i in range(nrow)]
chd_spd += [[0, i, ncol - 1, H2] for i in range(nrow)]
wel_spd = [(1, 0, 5, -5000.0)]
icelltype = [1]
arr = np.ones(shape2d)
arr[:, 4:7] = 0
icelltype.append(arr)


def build_models(idx, test):
    name = cases[idx]
    sim_ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=str(sim_ws),
        exe_name="mf6",
    )
    flopy.mf6.ModflowTdis(sim, nper=nper)
    linear_acceleration = "bicgstab"
    newtonoptions = "newton under_relaxation"

    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        no_ptcrecord="ALL",
        linear_acceleration=linear_acceleration,
        outer_maximum=nouter,
        outer_dvclose=hclose_outer,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions=newtonoptions,
    )
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
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=icelltype,
        k=k11,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=H1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfwel(gwf, maxbound=1, stress_period_data=wel_spd)

    head_filerecord = f"{name}.hds"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    mf6sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    heads = mf6sim.get_model().output.head().get_data()
    for k in range(nlay):
        b = botm[k]
        arr = heads[k]
        for j in (4, 5, 6):
            assert arr[:, j] < b, f"head > botm in columns in ({k + 1},{1},{j + 1})"
    b = botm[-1]
    arr = heads[nlay - 1]
    for j in (0, 1, 2, 3, 7, 8, 9, 10):
        assert arr[:, j] >= b, f"head < botm in columns in ({k + 1},{1},{j + 1})"


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
