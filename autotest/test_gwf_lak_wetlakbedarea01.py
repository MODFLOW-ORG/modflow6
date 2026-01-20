"""
A simple 2 layer by 1 row by 2 column model.  Upper-right cell is the only
active LAK cell.  Lake starts out initially dry and then is wetted by a
rising water table.  A constant head boundary in the lower left corner cell
is used to raise water table. This autotest checks to ensure that the wetted
areas between the lake and the 2 connected cells (1 vertical, 1 horizontal)
is correct.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-1cellkbd"]

# Model units
length_units = "feet"
time_units = "days"

# Model Parameters

nlay = 2
nrow = 1
ncol = 2
strt = 9.96
k11 = 130.0
k33 = [1179.0, 1179.0]
ss = 3e-4
sy = 0.2
chd_lr = 9.0
lak_strt = 9.0  # Starting lake stage
lak_bedleak = 10.0  # Lakebed leakance

idomain = np.full((nlay, nrow, ncol), 1)
idomain[0, 0, 1] = 0  # deactivate upper-right corner of 2x1x2 model

top = 20.0
botm = [10.0, 0.0]

# define delr and delc
delr = 10.0
delc = 10.0

# Timing
perlen = [1] * 10
nper = len(perlen)
nstp = 1
tsmult = 1.0

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp, tsmult))

# set dt0, dtmin, dtmax, dtadj, dtfailadj
dt0 = 5
dtmin = 1.001e-5
dtmax = 10.0
dtadj = 2.0
dtfailadj = 5.0

# Prepare constant head boundary data information
chd_spd = {}
chd_inc = [
    9.999999,
    10.0,
    10.000001,
    10.00001,
    10.0001,
    10.001,
    10.01,
    10.1,
    10.11,
    10.12,
]
for i, t in enumerate(range(len(perlen))):
    chd_spd.update({i: [nlay - 1, nrow - 1, 0, chd_inc[i]]})

# Prepare LAK package input
use_embedded_lak = False
if use_embedded_lak:
    idomain[0, 0, 1] = 1

lak_spd = [
    [0, "rainfall", 0.0],
    [0, "evaporation", 0.0],
]
lak_tab = [
    [10.1, 0.0, 0.0, 0.0],
    [10.11, 0.0002272, 0.001, 0.001],
    [10.12, 0.0006799, 0.01, 0.01],
    [10.13, 0.0014466, 0.1, 0.1],
    [10.15, 0.0036611, 0.2, 0.2],
    [10.2, 0.011438787, 0.3, 0.3],
    [10.25, 0.0280716, 0.4, 0.4],
    [10.3, 0.068889924, 0.5, 0.5],
    [10.35, 0.1690610, 0.6, 0.6],
    [10.4, 0.4148885490, 0.7, 0.7],
]

# Set solver parameters
nouter = 500
ninner = 100
hclose = 1e-6
rclose = 1e-6
relax = 0.97


def resolve_lvl(stg, hd, toplay):
    ss = min(stg, toplay)
    hh = min(hd, toplay)
    thk = max(ss, hh)
    return thk


def calc_qSat(top, bot, thk):
    teps = 1e-6
    tbmin = 0.0
    b = top - bot
    if b > 0.0:
        if thk < bot:
            br = 0.0
        elif thk > top:
            br = 1.0
        else:
            br = (thk - bot) / b

        av = 1.0 / (1.0 - teps)
        bri = 1.0 - br
        if br < tbmin:
            br = tbmin

        if br < teps:
            y = av * 0.5 * (br * br) / teps
        elif br < (1.0 - teps):
            y = av * br + 0.5 * (1.0 - av)
        elif br < 1.0:
            y = 1.0 - ((av * 0.5 * (bri * bri)) / teps)
        else:
            y = 1.0

    return y


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
    ats_filerecord = None
    tdis = flopy.mf6.ModflowTdis(
        sim,
        ats_filerecord=ats_filerecord,
        nper=nper,
        perioddata=tdis_rc,
        time_units=time_units,
    )

    if True:
        ats_filerecord = gwfname + ".ats"
        atsperiod = [
            (1, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (2, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (3, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (4, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (5, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (6, dt0, dtmin, dtmax, dtadj, dtfailadj),
            (7, dt0, dtmin, dtmax, dtadj, dtfailadj),
        ]
        tdis.ats.initialize(
            maxats=len(atsperiod),
            perioddata=atsperiod,
            filename=ats_filerecord,
        )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION",
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
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=sy, ss=ss, transient={0: True})

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # Instantiate constant head boundary package
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # Instantiate LAK package
    lak_conn = []
    if use_embedded_lak:
        lak_conn.append([0, 0, (0, 0, 1), "embeddedv", lak_bedleak, 0.0, 0.0, 1.0, 0.0])
    else:
        lak_conn.append(
            [0, 0, (0, 0, 0), "horizontal", lak_bedleak, 10.0, 20.0, 10.0, 10.0]
        )
        lak_conn.append([0, 1, (1, 0, 1), "vertical", lak_bedleak, 0.0, 0.0, 0.0, 0.0])

    lak_packagedata = [0, lak_strt, len(lak_conn)]
    budpth = f"{gwfname}.lak.cbc"
    tab6_filename = f"{gwfname}.laktab"
    if use_embedded_lak:
        # LAK package input requires tables option when using embedded lakes.
        lak = flopy.mf6.ModflowGwflak(
            gwf,
            save_flows=True,
            print_stage=True,
            nlakes=1,
            noutlets=0,
            packagedata=lak_packagedata,
            connectiondata=lak_conn,
            perioddata=lak_spd,
            ntables=1,
            tables=[0, tab6_filename],
            budget_filerecord=budpth,
            time_conversion=86400,
            length_conversion=3.28081,
            surfdep=0.05,
            pname="LAK-1",
            filename=f"{gwfname}.lak",
        )
    else:
        # Don't need to use the "TABLES" option for non-embedded lakes
        lak = flopy.mf6.ModflowGwflak(
            gwf,
            save_flows=True,
            print_stage=True,
            nlakes=1,
            noutlets=0,
            packagedata=lak_packagedata,
            connectiondata=lak_conn,
            perioddata=lak_spd,
            budget_filerecord=budpth,
            time_conversion=86400,
            length_conversion=3.28081,
            # surfdep=0.05,
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
    if use_embedded_lak:
        tabinput = []
        for itm in lak_tab:
            tabinput.append([itm[0], itm[1], itm[2], itm[3]])

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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name

    # read flow results from model
    sim1 = flopy.mf6.MFSimulation.load(sim_ws=test.workspace, load_only=["dis"])
    gwf = sim1.get_model(gwfname)

    # get final lake stage
    lk_pth0 = os.path.join(test.workspace, f"{gwfname}.lak.obs.csv")
    lkstg = np.genfromtxt(lk_pth0, names=True, delimiter=",")
    lkstg_time = lkstg["time"].tolist()
    lkstg_val = lkstg["STAGE"].tolist()

    # Store only the values at the end of the time step
    indices = [i for i, val in enumerate(lkstg_time) if not val.is_integer()]
    for i in indices[::-1]:
        lkstg_time.pop(i)
        lkstg_val.pop(i)

    # Get heads
    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    hdobj = flopy.utils.binaryfile.HeadFile(fname, precision="double")
    hds = hdobj.get_alldata()

    # Get lake/gwf exchange information
    fname = gwfname + ".lak.cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    lakobj = flopy.utils.binaryfile.CellBudgetFile(fname, precision="double")
    wetted_lakebed_area = lakobj.get_data(text="gwf")

    wetted_out = []
    for i in np.arange(len(wetted_lakebed_area)):
        vals = []
        for j in np.arange(len(wetted_lakebed_area[i])):
            val = wetted_lakebed_area[i][j][-1]
            vals.append(val)
        wetted_out.append(vals)
    wetted_out = np.array(wetted_out)

    # Compare MF6 output to answer calculated here
    msg = (
        "Compare value written by MF6 to a value calculated here based on "
        "either lake stage or gw head"
    )
    for tm in np.arange(wetted_out.shape[0]):
        for conn in np.arange(wetted_out.shape[1]):
            stg = lkstg_val[tm]
            # horizontal connections are stored first
            if conn == 0:
                gwh = hds[tm, 0, 0, 0]
                thk = resolve_lvl(stg, gwh, top)
                sat = calc_qSat(top, botm[0], thk)
                wa = sat * delc * (top - botm[0])
                mf6_wa = wetted_out[tm, conn]
                assert np.isclose(mf6_wa, wa, atol=1e-6), msg

            # vertical connection analysis
            elif conn == 1:
                gwh = hds[tm, 1, 0, 1]
                if gwh >= botm[0] or stg >= botm[0]:
                    # For a wetted vertical connection, it doesn't matter
                    # which direction the gradient is in, the wetted area
                    # is always delr * delc
                    wa = delc * delr
                    assert wetted_out[tm, conn] == wa, msg

    # Lake is dry in the first two stress periods. Both horiz. and vert.
    # connections
    msg = (
        "Lake starts out dry and wets-up after stress period 2. Horizontal "
        "wetted lakebed area should be equal to 0.0 in the first two stress "
        "periods"
    )
    assert np.all(wetted_out[0:2, 0] == 0.0), msg

    msg = (
        "With stage rising in the lake continuously, so too should the wetted area "
        "of the horizontal connection"
    )
    monotonicIncrease = np.diff(wetted_out[2:, 0])
    assert np.all(monotonicIncrease > 0), msg


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
