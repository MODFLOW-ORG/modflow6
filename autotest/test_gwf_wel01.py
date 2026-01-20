"""
Test to confirm that the sum wel and wel-reduced observations is equal
to the specified well pumping rate when the AUTO_FLOW_REDUCE option is
specified.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["wel01"]

# set static data
nper = 1
tdis_rc = [
    (50.0, 50, 1.05),
]

nlay, nrow, ncol = 2, 1, 1
delr = delc = 10.0
top, botm = 10.0, [9.0, 8.0]
strt = top
hk = 1.0
sy = 0.1
wellq = 0.25

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="NEWTON UNDER_RELAXATION",
    )

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
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=hk,
        k33=hk,
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=sy,
        steady_state={0: False},
        transient={0: False},
    )

    # wel files
    obs = {
        "wel.obs.csv": [
            ["q", "wel", (0, 0, 0)],
            ["qred", "wel-reduction", (0, 0, 0)],
        ],
        "wel.obs.dup.csv": [
            ["qred", "wel-reduction", (0, 0, 0)],
            ["q", "wel", (0, 0, 0)],
        ],
    }
    wel_spd = {0: [[0, 0, 0, -wellq]]}
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        auto_flow_reduce="auto_flow_reduce 0.5",
        stress_period_data=wel_spd,
        afrcsv_filerecord=f"{name}.afr.csv",
    )
    welobs = wel.obs.initialize(print_input=True, continuous=obs)

    # output control
    oc = flopy.mf6.ModflowGwfoc(gwf, printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])

    return sim, None


def check_output(idx, test):
    # MODFLOW 6 observations
    dtol = 1e-9
    for file_name in ("wel.obs.csv", "wel.obs.dup.csv"):
        fpth = os.path.join(test.workspace, file_name)
        try:
            tc = np.genfromtxt(fpth, names=True, delimiter=",")
        except:
            assert False, f'could not load data from "{fpth}"'

        qtot = tc["Q"] + tc["QRED"]

        # calculate maximum absolute error
        diff = qtot + wellq
        diffmax = np.abs(diff).max()
        msg = f"maximum absolute well rates ({diffmax}) "

        if diffmax > dtol:
            test.success = False
            msg += f"exceeds {dtol}"
            assert diffmax < dtol, msg
        else:
            test.success = True
            print("    " + msg)

    # MODFLOW 6 AFR CSV output file
    fpth = os.path.join(test.workspace, "wel01.afr.csv")
    try:
        afroutput = np.genfromtxt(fpth, names=True, delimiter=",", deletechars="")
    except:
        assert False, f'could not load data from "{fpth}"'

    a1 = afroutput["rate-requested"]
    a2 = afroutput["rate-actual"] + afroutput["wel-reduction"]
    errmsg = (
        "Auto flow reduce requested rate must equal actual rate plus reduced rate.\n"
    )
    errmsg += f"{a1} /= {a2}"
    assert np.allclose(a1, a2), errmsg


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
