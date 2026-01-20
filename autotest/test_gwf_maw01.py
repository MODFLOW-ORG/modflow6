import os
from types import SimpleNamespace

import flopy
import numpy as np
import pytest

cases = ["maw01", "maw01nwt", "maw01nwtur"]
budtol = 1e-2
bud_lst = ["GWF_IN", "GWF_OUT", "RATE_IN", "RATE_OUT"]
krylov = ["CG", "BICGSTAB", "BICGSTAB"]
newton = [None, "NEWTON", "NEWTON UNDER_RELAXATION"]
nlay = 1
nrow = 1
ncol = 3
nper = 3
delr = 300
delc = 300
perlen = 3 * [1]
nstp = 3 * [1]
tsmult = 3 * [1]
well = SimpleNamespace(
    observations={"maw_obs.csv": [("mh1", "head", 1)]},
    packagedata=[[0, 0.1, 50.0, 100.0, "THIEM", 1]],
    connectiondata=[[0, 0, (0, 0, 1), 100.0, 50.0, 1.0, 0.1]],
    perioddata=[[0, "rate", 0.0]],
)
strt = 100
hk = 1
nouter = 100
ninner = 300
hclose = 1e-9
rclose = 1e-3
relaxation_factor = 1
compare = False


def build_model(idx, ws, mf6):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=mf6,
        sim_ws=ws,
    )

    # create tdis package
    tdis_rc = [(perlen[i], nstp[i], tsmult[i]) for i in range(nper)]
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )
    gwf.name_file.newtonoptions = newton[idx]

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
        linear_acceleration=krylov[idx],
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
        top=100.0,
        botm=0.0,
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
        iconvert=1,
        ss=0.0,
        sy=0.1,
        steady_state={0: True},
        # transient={1: False},
        filename=f"{name}.sto",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 100.0])
    chdlist0.append([(0, 0, 2), 100.0])

    chdlist1 = []
    chdlist1.append([(0, 0, 0), 25.0])
    chdlist1.append([(0, 0, 2), 25.0])

    chdspdict = {0: chdlist0, 1: chdlist1, 2: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # wel files
    # wel = flopy.mf6.ModflowGwfwel(gwf, print_input=True, print_flows=True,
    #                              maxbound=len(ws),
    #                              periodrecarray=wd6,
    #                              save_flows=False)
    # MAW
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
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(workspace):
    # MODFLOW 6 maw results
    fpth = os.path.join(workspace, "maw_obs.csv")
    tc = np.genfromtxt(fpth, names=True, delimiter=",")

    # create known results array
    tc0 = np.array([100.0, 25.0, 100.0])

    # calculate maximum absolute error
    diff = tc["MH1"] - tc0
    diffmax = np.abs(diff).max()
    dtol = 1e-9
    msg = f"maximum absolute maw head difference {diffmax}"
    assert diffmax < dtol, msg + f" exceeds tolerance {dtol}"
    print(msg)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    ws = str(function_tmpdir)
    sim, _ = build_model(idx, ws, targets["mf6"])
    sim.write_simulation()
    sim.run_simulation()
    check_output(ws)
