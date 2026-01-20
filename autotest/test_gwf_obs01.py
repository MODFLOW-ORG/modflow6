import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cell_dimensions = (300,)
cases = [f"gwf_obs01{chr(ord('a') + idx)}" for idx in range(len(cell_dimensions))]
h0, h1 = 1.0, 0.0


def get_strt_array(idx):
    nrow, ncol = cell_dimensions[idx], cell_dimensions[idx]
    return np.broadcast_to(np.linspace(h0, h1, ncol), (nrow, ncol))


def get_obs(idx):
    nrow, ncol = cell_dimensions[idx], cell_dimensions[idx]
    obs_lst = []
    for i in range(nrow):
        for j in range(ncol):
            node = i * ncol + j + 1
            obs_lst.append([node, "head", (0, i, j)])
    return {f"{cases[idx]}.gwf.obs.csv": obs_lst}


def get_obs_out(idx, test):
    fpth = os.path.join(test.workspace, f"{cases[idx]}.gwf.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
        return tc.view((float, len(tc.dtype.names)))[1:]
    except:
        assert False, f'could not load data from "{fpth}"'


def get_chd(idx):
    nrow, ncol = cell_dimensions[idx], cell_dimensions[idx]
    c = [[(0, i, 0), h0] for i in range(nrow)]
    c += [[(0, i, ncol - 1), h1] for i in range(nrow)]
    return {0: c}


def build_models(idx, test):
    nlay, nrow, ncol = 1, cell_dimensions[idx], cell_dimensions[idx]
    nper = 1
    perlen = [5.0]
    nstp = [1]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    laytyp = 0
    botm = [0.0]
    hk = 1.0

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
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        no_ptcrecord="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        relaxation_factor=relax,
    )

    # create gwf model
    gwfname = name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = True

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
    )
    flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=get_obs(idx))

    # initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=get_strt_array(idx))

    # node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, icelltype=laytyp, k=hk, k33=hk
    )

    # chd files
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=get_chd(idx),
        save_flows=False,
        print_flows=True,
        pname="CHD-1",
    )

    # output control
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("BUDGET", "LAST")])

    return sim, None


def check_output(idx, test):
    hres = get_strt_array(idx).flatten()
    obs = get_obs_out(idx, test)
    assert np.allclose(hres, obs), (
        "simulated head observations do not match with known solution."
    )


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
