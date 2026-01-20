"""
Test obs package to make sure that the header in output csv files  is
correct.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cell_dimensions = (300,)
cases = ["gwf_obs02"]
h0, h1 = 1.0, 0.0
nlay, nrow, ncol = 1, 10, 10


def build_models(idx, test):
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

    # build list of obs csv files to create
    obsdict = {}
    for i in range(nrow):
        obslst = [(f"h_{i}_{j}", "head", (0, i, j)) for j in range(ncol)]
        fname = f"{name}.{i}.obs.csv"
        obsdict[fname] = obslst

    flopy.mf6.ModflowUtlobs(gwf, pname="head_obs", digits=20, continuous=obsdict)

    # initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)

    # node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, icelltype=laytyp, k=hk, k33=hk
    )

    # chd files
    chdlist = [[(0, 0, 0), h0]]
    chdlist += [[(0, nrow - 1, ncol - 1), h1]]
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist,
        save_flows=False,
        print_flows=True,
        pname="CHD-1",
    )

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        printrecord=[("BUDGET", "LAST"), ("HEAD", "LAST")],
        saverecord=[("HEAD", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    headcsv = np.empty((nlay, nrow, ncol), dtype=float)
    for i in range(nrow):
        fname = f"{test.name}.{i}.obs.csv"
        print(f"Loading and testing {fname}")
        fname = os.path.join(test.workspace, fname)
        rec = np.genfromtxt(fname, names=True, delimiter=",", deletechars="")
        for j in range(ncol):
            obsname_true = f"h_{i}_{j}".upper()
            obsname_found = rec.dtype.names[j + 1].upper()
            errmsg = 'obsname in {} is incorrect.  Looking for "{}" but found "{}"'
            errmsg = errmsg.format(fname, obsname_true, obsname_found)
            assert obsname_true == obsname_found, errmsg
        headcsv[0, i, :] = np.array(rec.tolist()[1:])

    fn = os.path.join(test.workspace, f"{test.name}.hds")
    hobj = flopy.utils.HeadFile(fn)
    headbin = hobj.get_data()

    assert np.allclose(headcsv, headbin), "headcsv not equal head from binary file"


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
