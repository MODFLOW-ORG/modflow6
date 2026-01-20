"""
An adaptation of the LAK package problem 1 supplemented with an additional
layer that has variable thinkness to help test that the shared wetted area
between a lakebed and groundwater cells in contact with the lake are written
to the LAK cbc output file correctly.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-wetlkbd"]

# Model units
length_units = "feet"
time_units = "days"

# Model Parameters

nper = 1
nlay = 6
nrow = 17
ncol = 17
strt = 115.0
k11 = 130.0
k33 = [1179.0, 1179.0, 30.0, 30.0, 30.0, 30.0]
ss = 3e-4
sy = 0.2
H1 = 160.0  # Constant head on left side of model
H2 = 140.0  # Constant head on right side of model
recharge = 0.0116
etvrate = 0.0141
etvdepth = 15.0
lak_strt = 110.0  # Starting lake stage
lak_etrate = 0.0103  # Lake evaporation rate
lak_bedleak = 0.1  # Lakebed leakance

left_new_lay_elv = 140
right_new_lay_elv = 120
new_lay_bot_elv = np.linspace(left_new_lay_elv, right_new_lay_elv, ncol)

top = np.ones((nrow, ncol)) * 500

botm = []
botm1 = []
for i in np.arange(nrow):
    botm1.append(new_lay_bot_elv)

botm1 = np.array(botm1)
botm.append(botm1)
# Account for remaining bottom layers
botm_elevs = [107.0, 97.0, 87.0, 77.0, 67.0]
for i in np.arange(1, nlay):
    botm_lay = np.ones((nrow, ncol)) * botm_elevs[i - 1]
    botm.append(botm_lay)

botm = np.array(botm)

# define delr and delc
delr = np.array(
    [
        250.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.00,
        500.00,
        500.00,
        500.0,
        500.00,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        250.0,
    ]
)
delc = np.array(
    [
        250.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        500.00,
        500.00,
        500.00,
        500.0,
        500.00,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        250.0,
    ]
)

# Timing
tdis_ds = ((5000.0, 1, 1.0),)

# Define dimensions
extents = (0.0, delr.sum(), 0.0, delc.sum())
shape3d = (nlay, nrow, ncol)

# Create the array defining the lake location
lake_map = np.ones(shape3d, dtype=np.int32) * -1
lake_map[0, 6:11, 6:11] = 0
lake_map[1, 6:11, 6:11] = 0
lake_map[2, 7:10, 7:10] = 0
lake_map = np.ma.masked_where(lake_map < 0, lake_map)

# Prepare linearly varying evapotranspiration surface
xlen = delr.sum() - 0.5 * (delr[0] + delr[-1])
x = 0.0
s1d = H1 * np.ones(ncol, dtype=float)
for idx in range(1, ncol):
    x += 0.5 * (delr[idx - 1] + delr[idx])
    frac = x / xlen
    s1d[idx] = H1 + (H2 - H1) * frac

surf = np.tile(s1d, (nrow, 1))
surf[lake_map[0] == 0] = 0
surf[lake_map[1] == 0] = 0
surf[lake_map[2] == 0] = 0

# Prepare constant head boundary data information
chd_spd = []
for k in range(nlay):
    chd_spd += [[k, i, 0, H1] for i in range(nrow)]
    chd_spd += [[k, i, ncol - 1, H2] for i in range(nrow)]

# Prepare LAK package input
lak_spd = [
    [0, "rainfall", recharge],
    [0, "evaporation", lak_etrate],
]

# Prepare Rch array
rech = np.ones((nrow, ncol)) * recharge
rech[lake_map[0] == 0] = 0

# Set solver parameters
nouter = 500
ninner = 100
hclose = 1e-9
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

    else:
        if x < bot:
            y = 0.0
        else:
            y = 1.0

    return y


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
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=sy, ss=ss, steady_state=True)

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # Instantiate constant head boundary package
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # Instantiate recharge package
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    # Instantiate ET package
    flopy.mf6.ModflowGwfevta(gwf, surface=surf, rate=etvrate, depth=etvdepth)

    # Instantiate LAK package
    (idomain_wlakes, pakdata_dict, lak_conn) = flopy.mf6.utils.get_lak_connections(
        gwf.modelgrid, lake_map, bedleak=lak_bedleak
    )
    global lak_con
    lak_con = lak_conn
    lak_packagedata = [[0, lak_strt, pakdata_dict[0]]]
    budpth = f"{gwfname}.lak.cbc"
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
    gwf.dis.idomain = idomain_wlakes

    # Instantiate output control package
    head_filerecord = f"{gwfname}.hds"
    budget_filerecord = f"{gwfname}.cbc"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
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
    lkstg_val = lkstg["STAGE"]

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
    lak_wetted_interface_area = lakobj.get_data(text="gwf")

    checks_out = []
    for i in np.arange(len(lak_wetted_interface_area[0])):
        dat = lak_wetted_interface_area[0][i]
        checks_out.append(dat[3])

    # Assumes that the order in which connection data is passed to the LAK instantiation
    # is the same order that values are written to the binary output file (*.lak.cbc)
    msg = (
        "The wetted interfacial areas saved in the binary output file "
        "(.cbc) do not match the values calculated in the autotest script"
    )
    for ii, itm in enumerate(lak_con):
        k, i, j = itm[2]
        ctype = itm[3]
        if ctype[0] == "h":
            botelv = botm[k, i, j]
            if k == 0:
                topelv = top[i, j]
            else:
                topelv = botm[k - 1, i, j]

            gwhd = hds[0, k, i, j]
            thk = resolve_lvl(lkstg_val, gwhd, topelv)
            sat = calc_qSat(topelv, botelv, thk)
            width = itm[-1]
            warea = sat * ((topelv - botelv) * width)
        elif ctype[0] == "v":
            length = delr[j]
            width = delc[i]
            warea = length * width

        assert np.isclose(warea, checks_out[ii], atol=1e-5), msg


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
