"""
An autotest for checking that the expected behavior of RCH/ET packages
turns on or off as a lake rises or falls.  That is, when cells are exposed
because of a drying lake, then ET and RCH should occur.  Conversely, as
the lake level rises and floods a cell with ET and/or RCH specified, both
should deactivate.  Set up shown next.  Packages that are active within
cell are listed

Plan View:
---------

        +----------+----------+----------+----------+----------+
        |          |          |          |          |          |
        |          |          |          |          |          |
        |          |          |          |          |          |
        +----------+----------+----------+----------+----------+
        |          |    LAK   |    LAK   |    LAK   |          |
        |          |    ET    |    ET    |    ET    |          |
        |          |    RCH   |    RCH   |    RCH   |          |
        +----------+----------+----------+----------+----------+
        |          |    LAK   |          |    LAK   |          |
        |          |    ET    |    LAK   |    ET    |          |
        |          |    RCH   |    ONLY  |    RCH   |          |
        +----------+----------+----------+----------+----------+
        |          |    LAK   |    LAK   |    LAK   |          |
        |          |    ET    |    ET    |    ET    |          |
        |          |    RCH   |    RCH   |    RCH   |          |
        +----------+----------+----------+----------+----------+
        |          |          |          |          |          |
        |          |          |          |          |          |
        |          |          |          |          |          |
        +----------+----------+----------+----------+----------+

Profile View:
------------

        +----------+                                +----------+ 10.0 m
        |          |_____v__________________________|          |
        |          |   /     / Simulated Lake  /    |          |
        |          |  /     /     /     /     /    /|          |
        +----------+----------+  /     / +----------+----------+  9.0 m
        |          |          | /     /  |          |          |
        |          |    ET    |/     /   |    ET    |          |
        |          |    RCH   |     /    |    RCH   |          |
        +----------+----------+----------+----------+----------+  8.0 m
        |          |          |          |          |          |
        |          |          |          |          |          |
        .          .          .          .          .          .
        .   CHD    .          .          .          .   CHD    .
        .          .          .          .          .          .
        |          |          |          |          |          |
        |          |          |          |          |          |
        +----------+----------+----------+----------+----------+  0.0 m


"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-rcha", "lak-rch", "lak-et"]

# Model units
length_units = "meter"
time_units = "days"

# Model Parameters

nper = 4
nlay = 3
nrow = 5
ncol = 5
strt = 9.5
k11 = 130.0
k33 = [1179.0, 1179.0, 30.0]
ss = 3e-4
sy = 0.2
cnstnthd = [9.5, 8.75, 8.65, 9.15]  # Constant head on left side of model
recharge = 0.00116
etvdepth = 0.95
lak_strt = 9.5  # Starting lake stage
lak_etrate = 0.0  # Lake evaporation rate
lak_bedleak = 0.01  # Lakebed leakance

top = 10.0
botm = [9.0, 8.0, 0.0]

# define delr and delc
delr = 5.0
delc = 5.0

# Timing
tdis_ds = ((5000.0, 1, 1.0), (5000.0, 1, 1.0), (5000.0, 1, 1.0), (5000.0, 1, 1.0))

# Define lake location
shape3d = (nlay, nrow, ncol)
lake_map = np.zeros(shape3d, dtype=np.int32)
lake_map[0, 1 : nrow - 1, 1 : ncol - 1] = 1
lake_map[1, 2, 2] = 1

idomain = 1 - lake_map

# For the recharge specified with an array, set the layer to 2 (0-based idx is 1)
irch = np.array(
    [
        [0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0],
        [0, 1, 0, 1, 0],
        [0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0],
    ]
)

rch_lst = []
for i in range(irch.shape[0]):
    for j in range(irch.shape[1]):
        if irch[i, j] > 0:
            rch_lst.append([1, i, j, recharge])

rch_spd = {0: rch_lst}

# Prepare linearly varying evapotranspiration surface
etrchsurf = np.array(
    [
        [10.0, 10.0, 10.0, 10.0, 10.0],
        [10.0, 9.0, 9.0, 9.0, 10.0],
        [10.0, 9.0, 8.0, 9.0, 10.0],
        [10.0, 9.0, 9.0, 9.0, 10.0],
        [10.0, 10.0, 10.0, 10.0, 10.0],
    ]
)

etvrate = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0141, 0.0141, 0.0141, 0.0],
        [0.0, 0.0141, 0.0, 0.0141, 0.0],
        [0.0, 0.0141, 0.0141, 0.0141, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ]
)

# Apply constant head boundary data to the corner cells
chd_spd = {}
for tm in range(nper):
    for k in range(nlay):
        chd = []
        chd += [[k, i, j, cnstnthd[tm]] for i in [0, 4] for j in [0, 4]]
        chd_spd.update({tm: chd})

# Prepare LAK package input
lak_spd = [
    [0, "rainfall", 0.0],
    [0, "evaporation", 0.0],
]

# LAK package input
lak_tab = [
    [8.0, 0.0, 0.0],
    [8.1, 0.0167, 0.1118],
    [8.2, 0.0667, 0.4472],
    [8.3, 0.1500, 1.0062],
    [8.4, 0.2667, 1.7888],
    [8.5, 0.4167, 2.7950],
    [8.6, 0.6000, 4.0248],
    [8.7, 0.8167, 5.4782],
    [8.8, 1.0667, 7.1552],
    [8.9, 1.3500, 9.0558],
    [9.0, 1.6667, 11.180],
    [9.1, 2.2000, 13.6396],
    [9.2, 2.8000, 19.0060],
    [9.3, 3.4667, 27.7264],
    [9.4, 4.2000, 40.2480],
    [9.5, 5.0000, 57.0180],
    [9.6, 5.8667, 78.4836],
    [9.7, 6.8000, 105.092],
    [9.8, 7.8000, 137.2904],
    [9.9, 8.8667, 175.526],
    [10.0, 10.00, 220.246],
]

# Prepare LAK package input
# <ifno> <iconn> <cellid(ncelldim)> <claktype> <bedleak> <belev> <telev> ...
#        <connlen> <connwidth>
ifno = 0
connectiondata = []
# "north side"
connectiondata.append(
    [ifno, 0, (0, 0, 1), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 1, (0, 0, 2), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 2, (0, 0, 3), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 3, (1, 1, 2), "horizontal", lak_bedleak, 8.0, 9.0, delc / 2.0, delr]
)
# "east side"
connectiondata.append(
    [ifno, 4, (0, 1, 4), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 5, (0, 2, 4), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 6, (0, 3, 4), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 7, (1, 2, 3), "horizontal", lak_bedleak, 8.0, 9.0, delr / 2.0, delc]
)
# "south side"
connectiondata.append(
    [ifno, 8, (0, 4, 1), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 9, (0, 4, 2), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 10, (0, 4, 3), "horizontal", lak_bedleak, 9.0, 10.0, delc / 2.0, delr]
)
connectiondata.append(
    [ifno, 11, (1, 3, 2), "horizontal", lak_bedleak, 8.0, 9.0, delc / 2.0, delr]
)
# "west side"
connectiondata.append(
    [ifno, 12, (0, 1, 0), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 13, (0, 2, 0), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 14, (0, 3, 0), "horizontal", lak_bedleak, 9.0, 10.0, delr / 2.0, delc]
)
connectiondata.append(
    [ifno, 15, (1, 2, 1), "horizontal", lak_bedleak, 8.0, 9.0, delr / 2.0, delc]
)
# bottoms
connectiondata.append(
    [ifno, 16, (1, 1, 1), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 17, (1, 1, 2), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 18, (1, 1, 3), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 19, (1, 2, 1), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 20, (1, 2, 3), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 21, (1, 3, 1), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 22, (1, 3, 2), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 23, (1, 3, 3), "vertical", lak_bedleak, 0.0, 8.0, 0.5, delc]
)
connectiondata.append(
    [ifno, 24, (2, 2, 2), "vertical", lak_bedleak, 0.0, 8.0, 4.0, delc]
)

nconn = len(connectiondata)
# <ifno> <strt> <nlakeconn> [<aux(naux)>] [<boundname>]
lak_packagedata = [
    [0, lak_strt, nconn, "lake1"],
]

# Prepare Rch array
rech = irch * recharge

# Set solver parameters
nouter = 500
ninner = 100
hclose = 1e-9
rclose = 1e-6
relax = 0.97

#
# MODFLOW 6 flopy GWF object
#


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    print(f"Building model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating time discretization
    flopy.mf6.ModflowTdis(
        sim, nper=len(tdis_ds), perioddata=tdis_ds, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="newton",
    )

    # Instantiating solver
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        linear_acceleration="bicgstab",
        outer_maximum=nouter,
        outer_dvclose=hclose,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])

    # Instantiate discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # Instantiate node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
        k33=k33,
    )

    # Instantiate gw storage package
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=sy, ss=ss, steady_state=False)

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # Instantiate constant head boundary package
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    if idx == 0:
        # Instantiate recharge as an array package
        flopy.mf6.ModflowGwfrcha(gwf, print_flows=True, recharge=rech, irch=irch)
    elif idx == 1:
        # Instantiate list-based recharge package
        flopy.mf6.ModflowGwfrch(gwf, print_flows=True, stress_period_data=rch_spd)
    else:
        # Instantiate ET package
        flopy.mf6.ModflowGwfevta(
            gwf, ievt=irch, surface=etrchsurf, rate=etvrate, depth=etvdepth
        )

    # Instantiate LAK package
    budpth = f"{gwfname}.lak.cbc"
    tab6_filename = f"{gwfname}.laktab"
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        save_flows=True,
        print_stage=True,
        boundnames=True,
        nlakes=len(lak_packagedata),
        noutlets=0,
        ntables=1,
        tables=[0, tab6_filename],
        packagedata=lak_packagedata,
        connectiondata=connectiondata,
        perioddata=lak_spd,
        budget_filerecord=budpth,
        pname="LAK-1",
        filename=f"{gwfname}.lak",
    )
    obs_file = f"{gwfname}.lak.obs"
    csv_file = obs_file + ".csv"
    obs_dict = {
        csv_file: [
            ("stage", "stage", (0,)),
        ]
    }
    lak.obs.initialize(
        filename=obs_file, digits=10, print_input=True, continuous=obs_dict
    )

    tabinput = []
    for itm in lak_tab:
        tabinput.append([itm[0], itm[1], itm[2]])

    laktab = flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(tabinput),
        ncol=len(tabinput[0]),
        table=tabinput,
        filename=tab6_filename,
        pname="LAK_tab",
        parent_file=lak,
    )

    # Instantiate output control package
    head_filerecord = f"{gwfname}.hds"
    budget_filerecord = f"{gwfname}.cbc"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        printrecord=[("HEAD", "LAST")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    msg0 = (
        "For the current simulation, recharge should only be simulated when "
        "recharge cells that are overlain by a lake are not submerged by "
        "the lake. Currently, recharge is simulated when it should not be."
    )
    msg1 = (
        "For the current simulation, recharge should only be simulated when "
        "recharge cells that are overlain by a lake are not submerged by "
        "the lake. Currently, recharge is not being simulated when it should be."
    )
    msg2 = (
        "For the current simulation, ET should only be simulated when "
        "ET cells that are overlain by a lake are not submerged by "
        "the lake. Currently, ET is simulated when it should not be."
    )
    msg3 = (
        "For the current simulation, ET should only be simulated when "
        "ET cells that are overlain by a lake are not submerged by "
        "the lake. Currently, ET is not being simulated when it should be."
    )

    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name
    ws = test.workspace

    # read flow results from model
    # sim1 = flopy.mf6.MFSimulation.load(sim_ws=test.workspace, load_only=["dis"])
    # gwf = sim1.get_model(gwfname)

    # Get heads
    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    hdobj = flopy.utils.binaryfile.HeadFile(fname, precision="double")
    hds = hdobj.get_alldata()

    # Get final lake stage
    lk_pth0 = os.path.join(ws, f"{gwfname}.lak.obs.csv")
    lkstg = np.genfromtxt(lk_pth0, names=True, delimiter=",")
    lkstg_val = lkstg["STAGE"]

    # Lift out budget items
    lst_pth = os.path.join(ws, f"{gwfname}.lst")
    budget_key = "VOLUME BUDGET FOR ENTIRE MODEL AT END OF TIME STEP"
    lst = flopy.utils.MfListBudget(lst_pth, budgetkey=budget_key)

    if idx == 0:
        rcha_in = lst.get_dataframes()[0]["RCHA_IN"].values.tolist()

        # recharge should only be simulated when cells overlain by a lake
        # and with rechrate > 0 is not submerged by the lake
        for tm in range(len(tdis_ds)):
            if lkstg_val[tm] > botm[0]:
                # rech cells inundated by overlying lake, ensure no rech
                assert rcha_in[tm] == 0.0, msg0
            elif lkstg_val[tm] < botm[0]:
                # rech cells exposed, ok for specified recharge to enter model
                assert rcha_in[tm] > 0.0, msg1

    elif idx == 1:
        rch_in = lst.get_dataframes()[0]["RCH_IN"].values.tolist()

        # recharge should only be simulated when cells overlain by a lake
        # and with rechrate > 0 is not submerged by the lake
        for tm in range(len(tdis_ds)):
            if lkstg_val[tm] > botm[0]:
                # rech cells inundated by overlying lake, ensure no rech
                assert rch_in[tm] == 0.0, msg0
            elif lkstg_val[tm] < botm[0]:
                # rech cells exposed, ok for specified recharge to enter model
                assert rch_in[tm] > 0.0, msg1

    elif idx == 2:
        et_out = lst.get_dataframes()[0]["EVTA_OUT"].values.tolist()

        # ET should only be simulated when a cell(s) with etrate > 0
        # is not submerged by the lake
        for tm in range(len(tdis_ds)):
            if lkstg_val[tm] > botm[0]:
                # rech cells inundated by overlying lake, ensure no rch
                assert et_out[tm] == 0.0, msg2
            elif lkstg_val[tm] < botm[0]:
                # rech cells exposed, ok for specified recharge to enter model
                assert et_out[tm] > 0.0, msg3


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
