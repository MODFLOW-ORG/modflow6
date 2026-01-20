import os
from types import SimpleNamespace

import flopy
import numpy as np
import pytest

cases = ["maw03a", "maw03b", "maw03c"]
budtol = 1e-2
bud_lst = ["GWF_IN", "GWF_OUT", "RATE_IN", "RATE_OUT"]


def well3(name):
    perioddata = {
        "maw03a": [(0, "rate", 2000.0)],
        "maw03b": [(0, "rate", 2000.0), (0, "head_limit", 0.4)],
        "maw03c": [(0, "rate", 2000.0), (0, "rate_scaling", 0.0, 1.0)],
    }
    wellbottom = -1000
    return SimpleNamespace(
        observations={
            f"{name}.maw.obs.csv": [
                ("m1head", "head", (0,)),
                ("m1rate", "rate", (0,)),
            ]  # is this index one-based? Not if in a tuple
        },
        packagedata=[[0, 0.15, wellbottom, 0.0, "THIEM", 1]],
        connectiondata=[[0, 0, (0, 50, 50), 0.0, wellbottom, 0.0, 0.0]],
        perioddata=perioddata[name],
    )


krylov = "CG"
nlay = 1
nrow = 101
ncol = 101
nper = 1
delr = 142
delc = 142
perlen = [1000]
nstp = [50]
tsmult = [1.2]
strt = 0
hk = 10
nouter = 100
ninner = 100
hclose = 1e-6
rclose = 1e-6
relaxation_factor = 1
compare = False


def build_model(idx, ws, mf6):
    top = 0.0
    botm = [-1000.0]

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
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
        linear_acceleration=krylov,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relaxation_factor,
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
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=0,
        ss=1.0e-5,
        sy=0.1,
        steady_state={0: False},
        transient={0: True},
        filename=f"{name}.sto",
    )

    # MAW
    well = well3(name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations=well.observations,
        packagedata=well.packagedata,
        connectiondata=well.connectiondata,
        perioddata=well.perioddata,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    # head observations
    obs_data0 = [("head_well_cell", "HEAD", (0, 0, 0))]
    obs_recarray = {f"{name}.obs.csv": obs_data0}
    obs = flopy.mf6.ModflowUtlobs(
        gwf,
        pname="head_obs",
        filename=f"{name}.obs",
        digits=15,
        print_input=True,
        continuous=obs_recarray,
    )

    return sim


def eval_results(name, workspace):
    # MODFLOW 6 maw results
    test_name = name
    fpth = os.path.join(workspace, f"{test_name}.maw.obs.csv")
    tc = np.genfromtxt(fpth, names=True, delimiter=",")

    if test_name.endswith("a"):
        # M1RATE should be 2000.
        msg = "The injection rate should be 2000. for all times"
        assert tc["M1RATE"].min() == tc["M1RATE"].max() == 2000, msg

    elif test_name.endswith("b"):
        # M1RATE should have a minimum value less than 200 and
        # M1HEAD should not exceed 0.400001
        msg = "Injection rate should fall below 200 and the head should notexceed 0.4"
        assert tc["M1RATE"].min() < 200.0, msg
        assert tc["M1HEAD"].max() < 0.400001, msg

    elif test_name.endswith("c"):
        # M1RATE should have a minimum value less than 800
        # M1HEAD should not exceed 1.0.
        msg = (
            "Min injection rate should be less than 800 and well "
            "head should not exceed 1.0"
        )
        assert tc["M1RATE"].min() < 800.0 and tc["M1HEAD"].max() < 1.0, msg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    ws = str(function_tmpdir)
    sim = build_model(idx, ws, targets["mf6"])
    sim.write_simulation()
    sim.run_simulation()
    eval_results(name, ws)
