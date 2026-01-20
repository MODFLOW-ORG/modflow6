import pathlib as pl

import flopy
import numpy as np
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = ["nr_ur01", "nr_ur02"]
data_path = project_root_path / "autotest/data/ex-gwf-bump/"
nper = 1
nlay = 1
nrow = 51
ncol = 51
xlen = 100.0
ylen = 100.0
top = 25.0
k11 = 1.0
H1 = 7.5
H2 = 2.5
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
botm = np.loadtxt(data_path / "bottom.txt").reshape(shape3d)
chd_spd = [[0, i, 0, H1] for i in range(nrow)]
chd_spd += [[0, i, ncol - 1, H2] for i in range(nrow)]
base_heads = flopy.utils.HeadFile(data_path / "results.hds.cmp").get_data()
mfsplit = None


def build_models(idx, test):
    name = cases[idx]
    if idx == 1:
        sim_ws = pl.Path(f"{test.workspace}/working")
    else:
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
        icelltype=1,
        k=k11,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=H1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    head_filerecord = f"{name}.hds"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        saverecord=[("HEAD", "ALL")],
    )

    if idx == 1:
        global mfsplit
        sim.write_simulation(silent=True)
        mfsplit = flopy.mf6.utils.Mf6Splitter(sim)
        split_array = np.tri(nrow, ncol).astype(int)
        new_sim = mfsplit.split_model(split_array)
        new_sim.set_sim_path(test.workspace)
        return new_sim, None
    else:
        return sim, None


def check_output(idx, test):
    mf6sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    if idx == 1:
        head_dict = {}
        for modelname in mf6sim.model_names:
            mnum = int(modelname.split("_")[-1])
            head_dict[mnum] = mf6sim.get_model(modelname).output.head().get_data()
        heads = mfsplit.reconstruct_array(head_dict)
    else:
        heads = mf6sim.get_model().output.head().get_data()
    msg = "head comparison failed"
    assert np.allclose(base_heads, heads), msg


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
