import math
import os

import flopy
import numpy as np
import pytest
from flopy.utils.lgrutil import Lgr
from framework import TestFramework

cases = ["mltmvr", "mltmvr5050", "mltmvr7525"]
sim_workspaces = []
gwf_names = []

# ----------------
# Universal input
# ----------------
numdays = 1
perlen = [1] * numdays
nper = len(perlen)
nstp = [1] * numdays
tsmult = [1.0] * numdays

icelltype = [1, 0, 0]

# Aquifer properties
hk = 1
k33 = 1

# Solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 0.97

# ------------------------------------------
# Static input associated with parent model
# ------------------------------------------
nlayp = 3
nrowp = 15
ncolp = 15
delrp = 1544.1 / ncolp
delcp = 1029.4 / nrowp
x = [round(x, 3) for x in np.linspace(50.0, 45.0, ncolp)]
topp = np.repeat(x, nrowp).reshape((15, 15)).T
z = [round(z, 3) for z in np.linspace(50.0, 0.0, nlayp + 1)]
botmp = [topp - z[len(z) - 2], topp - z[len(z) - 3], topp - z[0]]
idomainp = np.ones((nlayp, nrowp, ncolp), dtype=np.int32)
# Zero out where the child grid will reside
idomainp[0:2, 6:11, 2:8] = 0

xorigin = 2 * delrp
yorigin = 4 * delcp

# ------------------------------------------
# Common SFR data for all parent models
# ------------------------------------------

# Package_data information
sfrcells = [
    (0, 0, 1),
    (0, 1, 1),
    (0, 2, 1),
    (0, 2, 2),
    (0, 3, 2),
    (0, 4, 2),
    (0, 4, 3),
    (0, 5, 3),
    (0, 8, 8),
    (0, 8, 9),
    (0, 8, 10),
    (0, 8, 11),
    (0, 7, 11),
    (0, 7, 12),
    (0, 6, 12),
    (0, 6, 13),
    (0, 6, 14),
    (0, 5, 14),
]
rlen = [
    65.613029,
    72.488609,
    81.424789,
    35.850410,
    75.027390,
    90.887520,
    77.565651,
    74.860397,
    120.44695,
    112.31332,
    109.00368,
    91.234566,
    67.486000,
    24.603355,
    97.547943,
    104.97595,
    8.9454498,
    92.638367,
]
rwid = 5
rgrd1 = 0.12869035e-02
rgrd2 = 0.12780087e-02
rbtp = [
    49.409676,
    49.320812,
    49.221775,
    49.146317,
    49.074970,
    48.968212,
    48.859821,
    48.761742,
    45.550678,
    45.401943,
    45.260521,
    45.132568,
    45.031143,
    44.972298,
    44.894241,
    44.764832,
    44.692032,
    44.627121,
]
rbth = 1.5
rbhk = 0.1
man = 0.04
ustrf = 1.0
ndv = 0

# -----------------------------------------------
# Child model SFR data (common to all scenarios)
# -----------------------------------------------
connsc = []
for i in np.arange(89):
    if i == 0:
        connsc.append((i, -1 * (i + 1)))
    elif i == 88:
        connsc.append((i, i - 1))
    else:
        connsc.append((i, i - 1, -1 * (i + 1)))

# Package_data information
sfrcellsc = [
    (0, 0, 3),
    (0, 1, 3),
    (0, 1, 2),
    (0, 2, 2),
    (0, 2, 1),
    (0, 3, 1),
    (0, 4, 1),
    (0, 5, 1),
    (0, 6, 1),
    (0, 7, 1),
    (0, 7, 2),
    (0, 7, 3),
    (0, 7, 4),
    (0, 6, 4),
    (0, 5, 4),
    (0, 4, 4),
    (0, 3, 4),
    (0, 3, 5),
    (0, 3, 6),
    (0, 4, 6),
    (0, 4, 7),
    (0, 5, 7),
    (0, 5, 8),
    (0, 6, 8),
    (0, 7, 8),
    (0, 7, 7),
    (0, 8, 7),
    (0, 8, 6),
    (0, 8, 5),
    (0, 8, 4),
    (0, 9, 4),
    (0, 9, 3),
    (0, 10, 3),
    (0, 11, 3),
    (0, 12, 3),
    (0, 13, 3),
    (0, 13, 4),
    (0, 14, 4),
    (0, 14, 5),
    (0, 14, 6),
    (0, 13, 6),
    (0, 13, 7),
    (0, 12, 7),
    (0, 11, 7),
    (0, 11, 8),
    (0, 10, 8),
    (0, 9, 8),
    (0, 8, 8),
    (0, 7, 8),
    (0, 7, 9),
    (0, 6, 9),
    (0, 5, 9),
    (0, 4, 9),
    (0, 3, 9),
    (0, 2, 9),
    (0, 2, 10),
    (0, 1, 10),
    (0, 0, 10),
    (0, 0, 11),
    (0, 0, 12),
    (0, 0, 13),
    (0, 1, 13),
    (0, 2, 13),
    (0, 3, 13),
    (0, 4, 13),
    (0, 5, 13),
    (0, 6, 13),
    (0, 6, 12),
    (0, 7, 12),
    (0, 8, 12),
    (0, 9, 12),
    (0, 10, 12),
    (0, 11, 12),
    (0, 12, 12),
    (0, 12, 13),
    (0, 13, 13),
    (0, 13, 14),
    (0, 13, 15),
    (0, 12, 15),
    (0, 11, 15),
    (0, 10, 15),
    (0, 10, 16),
    (0, 9, 16),
    (0, 9, 15),
    (0, 8, 15),
    (0, 7, 15),
    (0, 6, 15),
    (0, 6, 16),
    (0, 6, 17),
]

rlenc = [
    24.637711,
    31.966246,
    26.376442,
    11.773884,
    22.921772,
    24.949730,
    23.878050,
    23.190311,
    24.762365,
    24.908625,
    34.366299,
    37.834534,
    6.7398176,
    25.150850,
    22.888292,
    24.630053,
    24.104542,
    35.873375,
    20.101446,
    35.636936,
    39.273537,
    7.8477302,
    15.480835,
    22.883194,
    6.6126003,
    31.995899,
    9.4387379,
    35.385513,
    35.470993,
    23.500074,
    18.414469,
    12.016913,
    24.691732,
    23.105467,
    23.700483,
    19.596104,
    5.7555680,
    34.423119,
    36.131992,
    7.4424477,
    35.565659,
    1.6159637,
    32.316132,
    20.131876,
    6.5242062,
    25.575630,
    25.575630,
    24.303566,
    1.9158504,
    21.931326,
    23.847176,
    23.432203,
    23.248718,
    23.455051,
    15.171843,
    11.196334,
    34.931976,
    4.4492774,
    36.034172,
    38.365566,
    0.8766859,
    30.059759,
    25.351671,
    23.554117,
    24.691738,
    26.074226,
    13.542957,
    13.303432,
    28.145079,
    24.373089,
    23.213642,
    23.298107,
    24.627758,
    27.715137,
    1.7645065,
    39.549232,
    37.144009,
    14.943290,
    24.851254,
    23.737432,
    15.967736,
    10.632832,
    11.425938,
    20.009295,
    24.641207,
    27.960585,
    4.6452723,
    36.717735,
    34.469074,
]
rwidc = 5
rgrdc = 0.14448310e-02
rbtpc = [
    48.622822,
    48.581932,
    48.539783,
    48.512222,
    48.487160,
    48.452576,
    48.417301,
    48.383297,
    48.348656,
    48.312775,
    48.269951,
    48.217793,
    48.185593,
    48.162552,
    48.127850,
    48.093521,
    48.058315,
    48.014984,
    47.974548,
    47.934284,
    47.880165,
    47.846127,
    47.829273,
    47.801556,
    47.780251,
    47.752357,
    47.722424,
    47.690044,
    47.638855,
    47.596252,
    47.565975,
    47.543991,
    47.517471,
    47.482941,
    47.449127,
    47.417850,
    47.399536,
    47.370510,
    47.319538,
    47.288059,
    47.256992,
    47.230129,
    47.205616,
    47.167728,
    47.148472,
    47.125282,
    47.088329,
    47.052296,
    47.033356,
    47.016129,
    46.983055,
    46.948902,
    46.915176,
    46.881439,
    46.853535,
    46.834484,
    46.801159,
    46.772713,
    46.743465,
    46.689716,
    46.661369,
    46.639019,
    46.598988,
    46.563660,
    46.528805,
    46.492130,
    46.463512,
    46.444118,
    46.414173,
    46.376232,
    46.341858,
    46.308254,
    46.273632,
    46.235821,
    46.214523,
    46.184677,
    46.129272,
    46.091644,
    46.062897,
    46.027794,
    45.999111,
    45.979897,
    45.963959,
    45.941250,
    45.908993,
    45.870995,
    45.847439,
    45.817558,
    45.766132,
]
rbthc = 1.5
rbhkc = 0.1
manc = 0.04
ustrfc = 1.0
ndvc = 0


# ---------------------------------------------------
# Scenario specific parent model SFR connection data
# ---------------------------------------------------
connsp_base = [
    (0, -1),
    (1, 0, -2),
    (2, 1, -3),
    (3, 2, -4),
    (4, 3, -5),
    (5, 4, -6),
    (6, 5, -7),
    (7, 6),
    (8, -9),
    (9, 8, -10),
    (10, 9, -11),
    (11, 10, -12),
    (12, 11, -13),
    (13, 12, -14),
    (14, 13, -15),
    (15, 14, -16),
    (16, 15, -17),
    (17, 16),
]

connsp_mvr = [
    (0, -1),
    (1, 0, -2),
    (2, 1, -3),
    (3, 2, -4),
    (4, 3, -5),
    (5, 4, -6),
    (6, 5, -7),
    (7, 6),
    (8, -9),
    (9, 8, -10),
    (10, 9, -11),
    (11, 10, -12),
    (12, 11, -13),
    (13, 12, -14),
    (14, 13, -15),
    (15, 14),
    (16, -17),
    (17, 16),
]

scen_conns = [connsp_base, connsp_mvr, connsp_mvr]

# ---------------------------------------------------
# Scenario specific MVR connection data
# (for simulation- and gwf-level MVRs)
# ---------------------------------------------------
# parent model gwf mvr
# static data
mvrpack = [["WEL-1"], ["SFR-parent"]]
maxpackages = len(mvrpack)
maxmvr = 10

# scenario specific data
parent_mvr_frac = [None, 0.50, 0.75]


def get_parent_mvr_info(frac):
    # return the appropriate mvr info for the current scenario
    mvrperioddata = [("WEL-1", 0, "SFR-parent", 10, "FACTOR", 1.0)]
    if frac is not None:
        mvrperioddata.append(("SFR-parent", 15, "SFR-parent", 16, "FACTOR", frac))

    mvrspd = {0: mvrperioddata}

    return mvrspd


# child model gwf mvr_scen
mvrpackc = [["WEL-2"], ["SFR-child"]]
maxpackagesc = len(mvrpackc)
mvrperioddatac = [("WEL-2", 0, "SFR-child", 53, "FACTOR", 1.0)]
mvrspdc = {0: mvrperioddatac}


# simulation mvr
def generate_parentmod_sfr_input(conns):
    pkdat = []
    for i in np.arange(len(rlen)):
        if i < 8:
            rgrd = rgrd1
        else:
            rgrd = rgrd2

        cln_list = len(
            [itm for itm in conns[i] if itm is not None and itm is not np.nan]
        )
        ncon = cln_list - 1
        pkdat.append(
            (
                i,
                sfrcells[i],
                rlen[i],
                rwid,
                rgrd,
                rbtp[i],
                rbth,
                rbhk,
                man,
                ncon,
                ustrf,
                ndv,
            )
        )

    return pkdat


def generate_childmod_sfr_input():
    pkdatc = []
    for i in np.arange(len(rlenc)):
        cln_list = len(
            [itm for itm in connsc[i] if itm is not None and itm is not np.nan]
        )
        nconc = cln_list - 1
        pkdatc.append(
            (
                i,
                sfrcellsc[i],
                rlenc[i],
                rwidc,
                rgrdc,
                rbtpc[i],
                rbthc,
                rbhkc,
                manc,
                nconc,
                ustrfc,
                ndvc,
            )
        )

    return pkdatc


def instantiate_base_simulation(sim_ws, gwfname, gwfnamec):
    # All pckgs between 3 test models the same except for parent model SFR input
    # static model data
    sim_workspaces.append(sim_ws)
    gwf_names.append(gwfname)
    sim = flopy.mf6.MFSimulation(
        sim_name="gwf",
        version="mf6",
        exe_name="mf6",
        sim_ws=sim_ws,
        continue_=False,
    )

    # Instantiate time discretization package
    tdis_rc = []
    for i in range(len(perlen)):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # Instantiate the gwf model (parent model)
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="NEWTON",
        model_nam_file=f"{gwfname}.nam",
    )

    # Create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    # Instantiate the discretization package
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlayp,
        nrow=nrowp,
        ncol=ncolp,
        delr=delrp,
        delc=delcp,
        top=topp,
        botm=botmp,
        idomain=idomainp,
        filename=f"{gwfname}.dis",
    )

    # Instantiate initial conditions package
    strt = [topp - 0.25, topp - 0.25, topp - 0.25]
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # Instantiate node property flow package
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        alternative_cell_averaging="AMT-LMK",
        icelltype=icelltype,
        k=hk,
        k33=k33,
        save_specific_discharge=False,
        filename=f"{gwfname}.npf",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # Instantiate constant head package
    rowList = np.arange(0, nrowp).tolist()
    layList = np.arange(0, nlayp).tolist()
    chdspd_left = []
    chdspd_right = []

    # Loop through rows, the left & right sides will appear in separate,
    # dedicated packages
    hd_left = 49.75
    hd_right = 44.75
    for l in layList:
        for r in rowList:
            # first, do left side of model
            chdspd_left.append([(l, r, 0), hd_left])
            # finally, do right side of model
            chdspd_right.append([(l, r, ncolp - 1), hd_right])

    chdspd = {0: chdspd_left}
    chd1 = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
        filename=f"{gwfname}.chd1.chd",
    )
    chdspd = {0: chdspd_right}
    chd2 = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-2",
        filename=f"{gwfname}.chd2.chd",
    )

    welspd_mf6 = []
    #                 [(layer,   row, column),  flow]
    welspd_mf6.append([(3 - 1, 8 - 1, 10 - 1), -5.0])
    wel_mf6_spd = {0: welspd_mf6}
    maxbound = len(welspd_mf6)
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=False,
        print_flows=True,
        maxbound=maxbound,
        mover=True,
        auto_flow_reduce=0.1,
        stress_period_data=wel_mf6_spd,  # wel_spd established in the MVR setup
        boundnames=False,
        save_flows=True,
        pname="WEL-1",
        filename=f"{gwfname}.wel",
    )

    # ---------------------------
    # Now work on the child grid
    # ---------------------------
    ncpp = 3
    ncppl = [3, 3, 0]

    lgr = Lgr(
        nlayp,
        nrowp,
        ncolp,
        delrp,
        delcp,
        topp,
        botmp,
        idomainp,
        ncpp=ncpp,
        ncppl=ncppl,
        xllp=0.0,
        yllp=0.0,
    )

    # Get child grid info:
    delrc, delcc = lgr.get_delr_delc()
    idomainc = lgr.get_idomain()  # child idomain
    topc, botmc = lgr.get_top_botm()  # top/bottom of child grid

    # Instantiate the gwf model (child model)
    gwfc = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfnamec,
        save_flows=True,
        newtonoptions="NEWTON",
        model_nam_file=f"{gwfnamec}.nam",
    )

    # Instantiate the discretization package
    child_dis_shp = lgr.get_shape()
    nlayc = child_dis_shp[0]
    nrowc = child_dis_shp[1]
    ncolc = child_dis_shp[2]
    disc = flopy.mf6.ModflowGwfdis(
        gwfc,
        nlay=nlayc,
        nrow=nrowc,
        ncol=ncolc,
        delr=delrc,
        delc=delcc,
        top=topc,
        botm=botmc,
        idomain=idomainc,
        xorigin=xorigin,
        yorigin=yorigin,
        filename=f"{gwfnamec}.dis",
    )

    # Instantiate initial conditions package
    strtc = [
        topc - 0.25,
        topc - 0.25,
        topc - 0.25,
        topc - 0.25,
        topc - 0.25,
        topc - 0.25,
    ]
    icc = flopy.mf6.ModflowGwfic(gwfc, strt=strtc, filename=f"{gwfnamec}.ic")

    # Instantiate node property flow package
    icelltypec = [1, 0, 0, 0, 0, 0]
    npfc = flopy.mf6.ModflowGwfnpf(
        gwfc,
        save_flows=False,
        alternative_cell_averaging="AMT-LMK",
        icelltype=icelltypec,
        k=hk,
        k33=k33,
        save_specific_discharge=False,
        filename=f"{gwfnamec}.npf",
    )

    # output control
    occ = flopy.mf6.ModflowGwfoc(
        gwfc,
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    welspd_mf6c = []
    #                 [(layer,   row, column),  flow]
    welspd_mf6c.append([(6 - 1, 4 - 1, 9 - 1), -10.0])
    wel_mf6_spdc = {0: welspd_mf6c}
    maxboundc = len(welspd_mf6c)
    welc = flopy.mf6.ModflowGwfwel(
        gwfc,
        print_input=False,
        print_flows=True,
        maxbound=maxboundc,
        mover=True,
        auto_flow_reduce=0.1,
        stress_period_data=wel_mf6_spdc,  # wel_spd established in the MVR setup
        boundnames=False,
        save_flows=True,
        pname="WEL-2",
        filename=f"{gwfnamec}.wel",
    )

    # exchange data
    exchange_data = lgr.get_exchange_data()

    # Establish GWF-GWF exchange
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        print_flows=True,
        print_input=True,
        exgmnamea=gwfname,
        exgmnameb=gwfnamec,
        nexg=len(exchange_data),
        exchangedata=exchange_data,
        pname="EXG-1",
        filename="gwf.exg",
    )

    return sim, gwf, gwfc


def add_parent_sfr(gwf, gwfname, conns):
    # Instantiate a scenario-specific sfr package
    pkdat = generate_parentmod_sfr_input(conns)
    sfrspd = {0: [[0, "INFLOW", 40.0]]}
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=False,
        print_flows=True,
        mover=True,
        pname="SFR-parent",
        length_conversion=1.0,
        time_conversion=86400.0,
        boundnames=False,
        nreaches=len(conns),
        packagedata=pkdat,
        connectiondata=conns,
        perioddata=sfrspd,
        filename=f"{gwfname}.sfr",
    )


def add_child_sfr(gwfc, gwfnamec):
    # Instantiate child model sfr package (same for all scenarios)
    pkdatc = generate_childmod_sfr_input()
    sfrspd = {0: [[0, "INFLOW", 0.0]]}
    sfrc = flopy.mf6.ModflowGwfsfr(
        gwfc,
        print_stage=False,
        print_flows=True,
        mover=True,
        pname="SFR-child",
        length_conversion=1.0,
        time_conversion=86400.0,
        boundnames=False,
        nreaches=len(connsc),
        packagedata=pkdatc,
        connectiondata=connsc,
        perioddata=sfrspd,
        filename=f"{gwfnamec}.sfr",
    )


def add_parent_mvr(gwf, gwfname, frac):
    # get scenario specific mvr data
    mvrspd = get_parent_mvr_info(frac)
    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackages,
        packages=mvrpack,
        perioddata=mvrspd,
        filename=f"{gwfname}.mvr",
    )


def add_child_mvr(gwfc, gwfnamec):
    mvrc = flopy.mf6.ModflowGwfmvr(
        gwfc,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackagesc,
        packages=mvrpackc,
        perioddata=mvrspdc,
        filename=f"{gwfnamec}.mvr",
    )


def add_sim_mvr(sim, gwfname, gwfnamec, remaining_frac=None):
    # simulation-level mvr data
    mvrpack_sim = [[gwfname, "SFR-parent"], [gwfnamec, "SFR-child"]]
    maxpackages_sim = len(mvrpack_sim)

    # Set up static SFR-to-SFR connections that remain fixed for entire simulation
    if remaining_frac is not None:
        sim_mvr_perioddata = [  # don't forget to use 0-based values
            [
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                7,
                mvrpack_sim[1][0],
                mvrpack_sim[1][1],
                0,
                "FACTOR",
                1.00,
            ],
            [
                mvrpack_sim[1][0],
                mvrpack_sim[1][1],
                88,
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                8,
                "FACTOR",
                1.00,
            ],
            [
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                15,
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                16,
                "FACTOR",
                remaining_frac,
            ],
        ]
    else:
        sim_mvr_perioddata = [  # don't forget to use 0-based values
            [
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                7,
                mvrpack_sim[1][0],
                mvrpack_sim[1][1],
                0,
                "FACTOR",
                1.00,
            ],
            [
                mvrpack_sim[1][0],
                mvrpack_sim[1][1],
                88,
                mvrpack_sim[0][0],
                mvrpack_sim[0][1],
                8,
                "FACTOR",
                1.00,
            ],
        ]

    mvrspd = {0: sim_mvr_perioddata}
    maxmvr = 3
    gwfgwf = sim.get_exchange_file("gwf.exg")
    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackages,
        packages=mvrpack_sim,
        perioddata=mvrspd,
        filename="gwf.mvr",
    )


def build_models(idx, test):
    scen_nm, conns, frac = (cases[idx], scen_conns[idx], parent_mvr_frac[idx])
    scen_nm_parent = "gwf_" + scen_nm + "_p"
    scen_nm_child = "gwf_" + scen_nm + "_c"
    sim, gwf, gwfc = instantiate_base_simulation(
        test.workspace, scen_nm_parent, scen_nm_child
    )
    # add the sfr packages
    add_parent_sfr(gwf, scen_nm_parent, conns)
    add_child_sfr(gwfc, scen_nm_child)
    # add the mover packages (simulation level and gwf level)
    add_parent_mvr(gwf, scen_nm_parent, frac)
    add_child_mvr(gwfc, scen_nm_child)
    if frac is not None:
        add_sim_mvr(sim, scen_nm_parent, scen_nm_child, 1.0 - frac)
    else:
        add_sim_mvr(sim, scen_nm_parent, scen_nm_child)

    return sim, None


def check_output(idx, test):
    gwf_srch_str1 = " SFR-PARENT PACKAGE - SUMMARY OF FLOWS FOR EACH CONTROL VOLUME"
    gwf_srch_str2 = " WATER MOVER PACKAGE (MVR) FLOW RATES   "
    sim_srch_str = " WATER MOVER PACKAGE (MVR) FLOW RATES "

    # cur_ws, gwfparent = ex[idx], gwf_names[idx]
    cur_ws = test.workspace
    gwfparent = "gwf_" + cases[idx] + "_p"
    with (
        open(os.path.join(cur_ws, gwfparent + ".lst"), "r") as gwf_lst,
        open(os.path.join(cur_ws, "mfsim.lst"), "r") as sim_lst,
    ):
        gwf_lst_lines = gwf_lst.readlines()
        sim_lst_lines = sim_lst.readlines()

    # Convert lists of lines to iterable objects
    gwf_lst = iter(gwf_lst_lines)
    sim_lst = iter(sim_lst_lines)

    # Peel mvr values from gwf lst file to be compared between scenarios
    done = False
    for line in gwf_lst:
        # adv file pointer to search line
        if gwf_srch_str1 in line:
            # once at the identified line, continue searching until sfr obj 17 queued
            while True:
                line = next(gwf_lst)
                m_arr = line.strip().split()
                if m_arr[0] == "18":
                    # store the 3rd value on the line
                    # (it should be the same across all scenarios)
                    parent_sfr_last_reach_flow = float(m_arr[2])
                    break

        # the second search string will only appear in the 50/50 and 75/25 scenarios
        if gwf_srch_str2 in line and idx > 0:
            # once at srch_str2, continue searching until 2nd mvr connection queued
            while True:
                line = next(gwf_lst)
                m_arr = line.strip().split()
                if m_arr[0] == "2":
                    parent_sfr_mvr_amount = float(m_arr[4])
                    done = True
                    break

        if done:
            break

    # now cycle through the simulation lst file
    # only the 50/50 and 75/25 scenarios will have the desired line
    if idx > 0:
        for line in sim_lst:
            if sim_srch_str in line:
                while True:
                    line = next(sim_lst)
                    m_arr = line.strip().split()
                    if m_arr[0] == "3":
                        sim_mvr_amount = float(m_arr[4])
                        done = True
                        break

    # perform the comparisons:
    #  o check flow entering last reach of parent model
    #    - should be the same across all 3 simulation
    #    - for current version of model, this amount was roughly 214.25
    #  o check the relative proportion of flows between gwf- and
    #    simulation-level mvrs
    #    - 50/50: ~107 units of flow in each
    #    - 75/25: 75% goes through the gwf mvr, 25% through the simulation mvr
    q_target = 214.25
    assert math.isclose(
        parent_sfr_last_reach_flow,
        q_target,
        rel_tol=0.1,
    ), (
        "Flow in the last reach of scenario "
        + cases[idx]
        + " = "
        + str(parent_sfr_last_reach_flow)
        + ", whereas the target flow "
        + " = "
        + str(q_target)
        + ".  Something changed, quitting."
    )

    # 50/50
    if idx == 1:
        gwf_transferred_50 = parent_sfr_mvr_amount / (
            parent_sfr_mvr_amount + sim_mvr_amount
        )
        sim_transferred_50 = sim_mvr_amount / (parent_sfr_mvr_amount + sim_mvr_amount)
        assert np.allclose(
            np.array([gwf_transferred_50, sim_transferred_50]),
            np.array([0.5, 0.5]),
            rtol=0.1,
        ), (
            "There should be a 50/50 split in the amount "
            "of water transferred by the GWF- and simulation-level MVRs."
        )

    # 75/25
    elif idx == 2:
        gwf_transferred_75 = parent_sfr_mvr_amount / (
            parent_sfr_mvr_amount + sim_mvr_amount
        )
        sim_transferred_75 = sim_mvr_amount / (parent_sfr_mvr_amount + sim_mvr_amount)
        assert np.allclose(
            np.array([gwf_transferred_75, sim_transferred_75]),
            np.array([0.75, 0.25]),
            rtol=0.1,
        ), (
            "There should be a 75/25 split in the amount of water "
            "transferred by the GWF- and simulation-level MVRs."
        )


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
