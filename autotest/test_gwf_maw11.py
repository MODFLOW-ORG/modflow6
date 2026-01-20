import os
from types import SimpleNamespace as Case

import flopy
import pytest
from framework import TestFramework

well1 = Case(
    observations={"maw_obs.csv": [("mh1", "head", 1)]},
    packagedata=[
        [0, 0.1, 0.0, 100.0, "THIEM", 1],
    ],
    connectiondata=[
        [0, 0, (0, 0, 1), 100.0, 0.0, 1.0, 0.1],
    ],
    perioddata={
        0: [
            [0, "rate", -0.0],
            [0, "status", "inactive"],
            [0, "rate_scaling", 1.0, 15.0],
        ],
        1: [
            [0, "status", "constant"],
            [0, "well_head", 80.0],
        ],
    },
)

sfr_packagedata = [
    [
        0,  # irch
        (0, 0, 2),  # cellid
        10,  # rlen
        2,  # rwid
        0.001,  # slp
        99.0,  # rtp
        1.0,  # rbth,
        0.0,  # rhk,
        0.035,  # roughness
        0,  # nconn
        1.0,  # ustrf
        0,  # ndv
    ]
]

connectiondata = [[0]]

sfrbndx = [[0, "INFLOW", 0.0]]
sfr_perioddata = {0: sfrbndx}

cases = ["maw11"]
krylov = "CG"
nlay = 1
nrow = 1
ncol = 3
nper = 5
delr = 300
delc = 300
perlen = 5 * [1]
nstp = 5 * [1]
tsmult = 5 * [1]
well = well1
strt = 95.0
hk = 1
nouter = 100
ninner = 300
hclose = 1e-9
rclose = 1e-3
relaxation_factor = 1
compare = False


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", sim_ws=ws)

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
        transient={0: True},
        filename=f"{name}.sto",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 95.0])

    chdlist1 = []
    chdlist1.append([(0, 0, 0), 95.0])

    chdspdict = {0: chdlist0, 1: chdlist1, 2: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # MAW
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        budget_filerecord=f"{name}.maw.cbc",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations=well.observations,
        packagedata=well.packagedata,
        connectiondata=well.connectiondata,
        perioddata=well.perioddata,
        mover=True,
        pname="MAW-1",
    )

    budpth = f"{name}.sfr.cbc"
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        budget_filerecord=budpth,
        mover=True,
        nreaches=len(sfr_packagedata),
        packagedata=sfr_packagedata,
        connectiondata=connectiondata,
        perioddata=sfr_perioddata,
        pname="SFR-1",
        filename=f"{name}.sfr",
    )

    packages = [("MAW-1",), ("SFR-1",)]
    mvr_perioddata = [("MAW-1", 0, "SFR-1", 0, "factor", 1.0)]
    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=len(mvr_perioddata),
        budget_filerecord=f"{name}.mvr.bud",
        maxpackages=len(packages),
        print_flows=True,
        packages=packages,
        perioddata=mvr_perioddata,
        pname="MVR-1",
        filename=f"{name}.mvr",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    print("checking items in MAW binary output...")

    # get results from listing file
    name = cases[idx]
    fpth = os.path.join(test.workspace, f"{os.path.basename(name)}.maw.cbc")
    if os.path.isfile(fpth):
        mawobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        kk = mawobj.get_kstpkper()
        times = mawobj.get_times()
        bud_names = mawobj.get_unique_record_names()

    nm_lst = []
    for itm in bud_names:
        nm_lst.append(str(itm.strip(), "utf-8"))

    assert "CONSTANT-TO-MVR" in nm_lst, (
        "Expected budget term not in MAW binary output file."
    )

    rtm = mawobj.get_data(text="RATE-TO-MVR")
    ctm = mawobj.get_data(text="CONSTANT-TO-MVR")


@pytest.mark.parametrize("idx, name", list(enumerate(cases)))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
