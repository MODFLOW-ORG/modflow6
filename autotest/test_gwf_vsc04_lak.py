"""
Simple single lake model.  Lake cut into top two layers of a 5 layer
model.  Model is loosely based on the first example problem in
Merritt and Konikow (2000) which also is one of the MT3D-USGS test
problems.  This test developed to isolate lake-aquifer interaction;
no SFR or other advanced packages.  Problem set up to have groundwater
pass through the lake: gw inflow on the left side, gw outflow on the
right side of the lake.  Uses constant stage boundary in the lake to
ensure desired flow conditions for testing budget changes with and
without VSC active.

starting groundwater temperature: 30.0
left chd boundary inflow temperature: 30.0
starting lake temperature: 4.0
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["no-vsc04-lak", "vsc04-lak"]
viscosity_on = [False, True]

# Model units
length_units = "m"
time_units = "days"

# model domain and grid definition
delr = [
    76.2,
    304.8,
    304.8,
    304.8,
    304.8,
    304.8,
    152.4,
    152.4,
    152.4,
    152.4,
    152.4,
    304.8,
    304.8,
    304.8,
    304.8,
    304.8,
    76.2,
]

delc = [
    76.2,
    304.8,
    304.8,
    304.8,
    304.8,
    304.8,
    152.4,
    152.4,
    152.4,
    152.4,
    152.4,
    304.8,
    304.8,
    304.8,
    304.8,
    304.8,
    76.2,
]

fixedstrthds = [
    35.052,
    34.9267,
    34.7216,
    34.5062,
    34.2755,
    34.0237,
    33.8143,
    33.6657,
    33.5077,
    33.3394,
    33.1599,
    32.8728,
    32.4431,
    31.9632,
    31.4353,
    30.8627,
    30.48,
]

nrow = len(delc)
ncol = len(delr)
top = np.ones((nrow, ncol)) * 35.6616
bot1 = np.ones_like(top) * 32.6136
bot2 = np.ones_like(top) * 29.5656
bot3 = np.ones_like(top) * 26.5176
bot4 = np.ones_like(top) * 23.4696
bot5 = np.ones_like(top) * 20.4216
botm = np.array([bot1, bot2, bot3, bot4, bot5])
nlay = botm.shape[0]
ibound = np.ones_like(botm)

# deactivate gw cells where lake cells are active
ibound[0, 6:11, 6:11] = 0  # layer 1
ibound[1, 7:10, 7:10] = 0  # layer 2

strthd = np.zeros_like(ibound)
for j in np.arange(ncol):
    strthd[:, :, j] = fixedstrthds[j]

# setup lake array
lakibnd = np.zeros_like(ibound)
lakibnd[0] = 1 - ibound[0]  # layer 1
lakibnd[1] = 1 - ibound[1]  # layer 2

# NPF parameters
k11 = 9.144  # = 30 ft/day
k33 = 0.9144  # = 30 ft/day
ss = 3e-4
sy = 0.20
hani = 1
laytyp = 1

# Package boundary conditions
chdl = 35.052
chdr = 30.48
viscref = 8.904e-4

# time params
transient = {0: True}
nstp = [100]
tsmult = [1.02]
perlen = [5000]

# solver params
nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-3, 1e-4, 0.97

# Transport related parameters
al = 1  # longitudinal dispersivity ($m$)
ath1 = al  # horizontal transverse dispersivity
atv = al  # vertical transverse dispersivity
mixelm = 0  # Upstream vs TVD (Upstream selected)
initial_temperature = 35.0  # Initial temperature (unitless)
porosity = 0.20  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rho_water = 1000  # Density of water ($kg/m^3$)
rho_solids = 2650  # Density of the aquifer material ($kg/m^3$)
C_p_w = 4180  # Heat Capacity of water ($J/kg/C$)
C_s = 880  # Heat capacity of the solids ($J/kg/C$)
D_m = K_therm / (porosity * rho_water * C_p_w)
rhob = (1 - porosity) * rho_solids  # Bulk density ($kg/m^3$)
K_d = C_s / (rho_water * C_p_w)  # Partitioning coefficient ($m^3/kg$)
leftTemp = 30.0  # Temperature of inflow from left constant head ($C$)

# Viscosity related parameters
tviscref = 20.0


def build_models(idx, test):
    global lak_lkup_dict

    # Base simulation and model name and workspace
    name = cases[idx]

    print(f"Building model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=test.workspace, exe_name="mf6", version="mf6"
    )

    tdis_rc = []
    for i in range(len(nstp)):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(
        sim, nper=len(nstp), perioddata=tdis_rc, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, newtonoptions="newton"
    )

    # Instantiating solver
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="cooley",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
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
        idomain=ibound,
        filename=f"{gwfname}.dis",
    )

    # Instantiate node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
        k33=k33,
    )

    # Instantiate storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        transient=transient,
    )

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strthd)

    # Instantiate viscosity package
    if viscosity_on[idx]:
        vsc_filerecord = f"{gwfname}.vsc.bin"
        vsc_pd = [(0, 0.0, tviscref, gwtname, "TEMPERATURE")]
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=viscref,
            viscosity_filerecord=vsc_filerecord,
            thermal_formulation="nonlinear",
            thermal_a2=10.0,
            thermal_a3=248.37,
            thermal_a4=133.16,
            nviscspecies=len(vsc_pd),
            packagedata=vsc_pd,
            pname="vsc",
            filename=f"{gwfname}.vsc",
        )

    # Instantiate output control package
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 17, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "LAST")],
    )

    # Instantiate constant head package
    # (for driving gw flow from left to right)
    chdlistl = []
    chdlistr = []
    for k in np.arange(nlay):
        for i in np.arange(nrow):
            # left side
            if botm[k, i, 0] <= chdl:
                chdlistl.append([(k, i, 0), chdl, leftTemp])
            # right side
            if botm[k, i, -1] <= chdr:
                chdlistr.append([(k, i, ncol - 1), chdr, 10.0])

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlistl,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-L",
        auxiliary="TEMPERATURE",
        filename=f"{gwfname}.left.chd",
    )

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlistr,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-R",
        auxiliary="TEMPERATURE",
        filename=f"{gwfname}.right.chd",
    )

    # Instantiate lake package
    lakeconnectiondata = []
    nlakecon = [0]  # Expand this to [0, 0, ...] for each additional lake
    ilakconn = -1
    lak_leakance = 0.1
    lak_lkup_dict = {}
    for k in [0, 1]:
        for i in range(nrow):
            for j in range(ncol):
                if lakibnd[k, i, j] == 0:
                    continue
                else:
                    ilak = int(lakibnd[k, i, j] - 1)
                    # back
                    if i > 0:
                        if lakibnd[k, i - 1, j] == 0 and ibound[k, i - 1, j] == 1:
                            ilakconn += 1
                            # by setting belev==telev, MF6 will automatically
                            # re-assign elevations based on cell dimensions
                            h = [
                                ilak,  # <ifno>
                                ilakconn,  # <iconn>
                                (k, i - 1, j),  # <cellid(ncelldim)>
                                "horizontal",  # <claktype>
                                lak_leakance,  # <bedleak>
                                0.0,  # <belev>
                                0.0,  # <telev>
                                delc[i] / 2.0,  # <connlen>
                                delr[j],  # <connwidth>
                            ]
                            lakeconnectiondata.append(h)
                            lak_lkup_dict.update({ilakconn: (k, i, j)})

                    # left
                    if j > 0:
                        if lakibnd[k, i, j - 1] == 0 and ibound[k, i, j - 1] == 1:
                            ilakconn += 1
                            h = [
                                ilak,
                                ilakconn,
                                (k, i, j - 1),
                                "horizontal",
                                lak_leakance,
                                0.0,
                                0.0,
                                delr[j] / 2.0,
                                delc[i],
                            ]
                            lakeconnectiondata.append(h)
                            lak_lkup_dict.update({ilakconn: (k, i, j)})

                    # right
                    if j < ncol - 1:
                        if lakibnd[k, i, j + 1] == 0 and ibound[k, i, j + 1] == 1:
                            ilakconn += 1
                            h = [
                                ilak,
                                ilakconn,
                                (k, i, j + 1),
                                "horizontal",
                                lak_leakance,
                                0.0,
                                0.0,
                                delr[j] / 2.0,
                                delc[i],
                            ]
                            lakeconnectiondata.append(h)
                            lak_lkup_dict.update({ilakconn: (k, i, j)})

                    # front
                    if i < nrow - 1:
                        if lakibnd[k, i + 1, j] == 0 and ibound[k, i + 1, j] == 1:
                            ilakconn += 1
                            h = [
                                ilak,
                                ilakconn,
                                (k, i + 1, j),
                                "horizontal",
                                lak_leakance,
                                0.0,
                                0.0,
                                delc[i] / 2.0,
                                delr[j],
                            ]
                            lakeconnectiondata.append(h)
                            lak_lkup_dict.update({ilakconn: (k, i, j)})

                # vertical
                if lakibnd[k, i, j] == 1 and ibound[k + 1, i, j] == 1:
                    ilakconn += 1
                    v = [
                        ilak,
                        ilakconn,
                        (k + 1, i, j),
                        "vertical",
                        lak_leakance,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                    ]
                    lakeconnectiondata.append(v)
                    lak_lkup_dict.update({ilakconn: (k, i, j)})

    strtStg = 33.75
    lakpackagedata = [[0, strtStg, len(lakeconnectiondata), 4.0, "lake1"]]
    lak_pkdat_dict = {"filename": "lak_pakdata.in", "data": lakpackagedata}

    lakeperioddata = {
        0: [
            (0, "STATUS", "CONSTANT"),  # RAINFALL 0.005 & 0.00504739035
            (0, "STAGE", 33.5),
        ]
    }

    lak_obs = {
        f"{gwfname}.lakeobs": [
            ("lakestage", "stage", "lake1"),
            ("gwexchng", "lak", "lake1"),
        ]
    }
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        auxiliary="TEMPERATURE",
        time_conversion=86400.0,
        print_stage=True,
        print_flows=True,
        budget_filerecord=gwfname + ".lak.bud",
        length_conversion=1.0,
        mover=False,
        pname="LAK-1",
        boundnames=True,
        nlakes=len(lakpackagedata),
        noutlets=0,
        packagedata=lak_pkdat_dict,
        connectiondata=lakeconnectiondata,
        perioddata=lakeperioddata,
        observations=lak_obs,
        filename=f"{gwfname}.lak",
    )

    # pull in the tabfile defining the lake stage, vol, & surface area
    fname = os.path.join("data", "vsc04-laktab", "stg-vol-surfarea.dat")
    tabinput = []
    with open(fname, "r") as f:
        # peel off the hdr line
        hdr = next(f)
        for line in f:
            m_arr = line.strip().split(",")
            #                        <stage>, <volume>,  <sarea>,
            tabinput.append([float(m_arr[0]), m_arr[1], m_arr[2]])

    tab6_filename = f"{gwfname}.laktab"
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(tabinput),
        ncol=3,
        table=tabinput,
        filename=tab6_filename,
        pname="LAK_tab",
        parent_file=lak,
    )

    # create gwt model
    # ----------------
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, model_nam_file=f"{gwtname}.nam")
    gwt.name_file.save_flows = True

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=ibound,
        filename=f"{gwtname}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    strtconc = leftTemp
    flopy.mf6.ModflowGwtic(gwt, strt=strtconc, filename=f"{gwtname}.ic")

    # Instantiate mobile storage and transfer package
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")

    # Instantiating MODFLOW 6 transport advection package
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()

    # Instantiate advection package
    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{gwtname}.adv")

    # Instantiate dispersion package
    flopy.mf6.ModflowGwtdsp(gwt, alh=al, ath1=ath1, atv=atv, filename=f"{gwtname}.dsp")

    # Instantiate source/sink mixing package
    sourcerecarray = [
        ("CHD-L", "AUX", "TEMPERATURE"),
        ("CHD-R", "AUX", "TEMPERATURE"),
    ]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm")

    # Instantiating MODFLOW 6 transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 17, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwtname}.oc",
    )

    # Instantiating MODFLOW 6 lake transport (lkt) package
    lktpackagedata = [(0, 4.0, "lake1")]

    lktperioddata = {0: [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 4.0)]}

    # note: for specifying lake number, use fortran indexing!
    lkt_obs = {
        f"{gwtname}.lakobs": [
            ("resTemp", "concentration", 1),
            ("resGwMassExchng", "lkt", "lake1"),
        ]
    }

    flopy.mf6.ModflowGwtlkt(
        gwt,  # Set time_conversion for use with Manning's eqn.
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="TEMPERATURE",
        budget_filerecord=gwtname + ".lkt.bud",
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=False,
        print_concentration=True,
        packagedata=lktpackagedata,
        lakeperioddata=lktperioddata,
        observations=lkt_obs,
        pname="LKT-1",
        filename=f"{gwtname}.lkt",
    )

    # GWF-GWT exchange
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name

    fname = gwfname + ".lak.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    budobj = flopy.utils.CellBudgetFile(fname, precision="double")
    outbud = budobj.get_data(text="             GWF")

    # Establish known answer:
    stored_ans = np.array(
        [
            [1.0, 9.20e1, 7.34424130e-01, 1.15181189e06],
            [1.0, 1.08e2, 2.44117249e00, 1.15181189e06],
            [1.0, 3.98e2, 2.19216490e01, 1.15181189e06],
            [1.0, 9.30e1, -6.78268488e-02, 1.15181189e06],
            [1.0, 3.99e2, 3.17042406e-01, 1.15181189e06],
            [1.0, 9.40e1, -7.47709994e-01, 1.15181189e06],
            [1.0, 4.00e2, -6.88938281e00, 1.15181189e06],
            [1.0, 9.50e1, -1.51336530e00, 1.15181189e06],
            [1.0, 4.01e2, -1.57614834e01, 1.15181189e06],
            [1.0, 9.60e1, -2.54095715e00, 1.15181189e06],
            [1.0, 1.14e2, -3.95961161e00, 1.15181189e06],
            [1.0, 4.02e2, -5.60853100e01, 1.15181189e06],
            [1.0, 1.25e2, 2.35138538e00, 1.15181189e06],
            [1.0, 4.15e2, 1.69275311e01, 1.15181189e06],
            [1.0, 1.31e2, -3.66648779e00, 1.15181189e06],
            [1.0, 4.19e2, -3.45225854e01, 1.15181189e06],
            [1.0, 1.42e2, 2.32550672e00, 1.15181189e06],
            [1.0, 4.32e2, 1.65405908e01, 1.15181189e06],
            [1.0, 1.48e2, -3.58087615e00, 1.15181189e06],
            [1.0, 4.36e2, -3.27154545e01, 1.15181189e06],
            [1.0, 1.59e2, 2.35138505e00, 1.15181189e06],
            [1.0, 4.49e2, 1.69275277e01, 1.15181189e06],
            [1.0, 1.65e2, -3.66648777e00, 1.15181189e06],
            [1.0, 4.53e2, -3.45225840e01, 1.15181189e06],
            [1.0, 1.76e2, 2.44117266e00, 1.15181189e06],
            [1.0, 1.94e2, 7.34424819e-01, 1.15181189e06],
            [1.0, 4.66e2, 2.19216459e01, 1.15181189e06],
            [1.0, 1.95e2, -6.78264346e-02, 1.15181189e06],
            [1.0, 4.67e2, 3.17038149e-01, 1.15181189e06],
            [1.0, 1.96e2, -7.47709250e-01, 1.15181189e06],
            [1.0, 4.68e2, -6.88938656e00, 1.15181189e06],
            [1.0, 1.97e2, -1.51336458e00, 1.15181189e06],
            [1.0, 4.69e2, -1.57614826e01, 1.15181189e06],
            [1.0, 1.82e2, -3.95961151e00, 1.15181189e06],
            [1.0, 1.98e2, -2.54095654e00, 1.15181189e06],
            [1.0, 4.70e2, -5.60853022e01, 1.15181189e06],
            [1.0, 3.99e2, 4.03508517e-03, 1.15181189e06],
            [1.0, 4.15e2, 2.15441304e-01, 1.15181189e06],
            [1.0, 7.05e2, 8.25215117e-01, 1.15181189e06],
            [1.0, 4.00e2, -8.76830539e-02, 1.15181189e06],
            [1.0, 7.06e2, -3.90793309e-01, 1.15181189e06],
            [1.0, 4.01e2, -2.00600698e-01, 1.15181189e06],
            [1.0, 4.19e2, -4.39378360e-01, 1.15181189e06],
            [1.0, 7.07e2, -2.43955302e00, 1.15181189e06],
            [1.0, 4.32e2, 2.10516610e-01, 1.15181189e06],
            [1.0, 7.22e2, 8.37390920e-01, 1.15181189e06],
            [1.0, 7.23e2, -6.85716153e-02, 1.15181189e06],
            [1.0, 4.36e2, -4.16378511e-01, 1.15181189e06],
            [1.0, 7.24e2, -1.73090130e00, 1.15181189e06],
            [1.0, 4.49e2, 2.15441262e-01, 1.15181189e06],
            [1.0, 4.67e2, 4.03503099e-03, 1.15181189e06],
            [1.0, 7.39e2, 8.25212943e-01, 1.15181189e06],
            [1.0, 4.68e2, -8.76831016e-02, 1.15181189e06],
            [1.0, 7.40e2, -3.90795105e-01, 1.15181189e06],
            [1.0, 4.53e2, -4.39378341e-01, 1.15181189e06],
            [1.0, 4.69e2, -2.00600688e-01, 1.15181189e06],
            [1.0, 7.41e2, -2.43955388e00, 1.15181189e06],
        ]
    )

    # talley some flows on the left and right sides of the lake for comparison
    # test
    left_chk_ans = []
    right_chk_ans = []
    left_chk_no_vsc = []
    right_chk_no_vsc = []
    left_chk_with_vsc = []
    right_chk_with_vsc = []

    if idx == 0:
        no_vsc_bud_last = np.array(outbud[-1].tolist())
        no_vsc_bud_np = np.array(no_vsc_bud_last.tolist())

        for ii in np.arange(stored_ans.shape[0]):
            k, i, j = lak_lkup_dict[ii]

            # left side of lake
            if j < 7:
                if no_vsc_bud_np[ii, 2] > 0 and stored_ans[ii, 2] > 0:
                    left_chk_no_vsc.append(no_vsc_bud_np[ii, 2])
                    left_chk_ans.append(stored_ans[ii, 2])

            # right side of lake
            if j > 9:
                if no_vsc_bud_np[ii, 2] < 0 and stored_ans[ii, 2] < 0:
                    right_chk_no_vsc.append(no_vsc_bud_np[ii, 2])
                    right_chk_ans.append(stored_ans[ii, 2])

        # Check that all the flows entering the lak in the 'with vsc' model are greater
        # than their 'no vsc' counterpart
        assert np.allclose(
            np.array(left_chk_ans), np.array(left_chk_no_vsc), atol=1e-3
        ), "Lake inflow in no-VSC LAK simulation do not match established solution."

        # Check that all the flows leaving the lak in the 'with vsc' model are less
        # than their 'no vsc' counterpart (keep in mind values are negative, which
        # affects how the comparison is made)
        assert np.allclose(
            np.array(right_chk_ans), np.array(right_chk_no_vsc), atol=1e-3
        ), "Lake outflow in no-VSC LAK simulation do not match established solution."

    elif idx == 1:
        with_vsc_bud_last = np.array(outbud[-1].tolist())
        with_vsc_bud_np = np.array(with_vsc_bud_last.tolist())

        for ii in np.arange(stored_ans.shape[0]):
            k, i, j = lak_lkup_dict[ii]

            # left side of lake
            if j < 7:
                if stored_ans[ii, 2] > 0 and with_vsc_bud_np[ii, 2] > 0:
                    left_chk_no_vsc.append(stored_ans[ii, 2])
                    left_chk_with_vsc.append(with_vsc_bud_np[ii, 2])

            # right side of lake
            if j > 9:
                if stored_ans[ii, 2] < 0 and with_vsc_bud_np[ii, 2] < 0:
                    right_chk_no_vsc.append(stored_ans[ii, 2])
                    right_chk_with_vsc.append(with_vsc_bud_np[ii, 2])

        # Check that all the flows entering the lak in the 'with vsc' model are greater
        # than their 'no vsc' counterpart
        assert np.greater(
            np.array(left_chk_with_vsc), np.array(left_chk_no_vsc)
        ).all(), "Lake inflow did no increase with VSC turned on and should have."

        # Check that all the flows leaving the lak in the 'with vsc' model are less
        # than their 'no vsc' counterpart (keep in mind values are negative, which
        # affects how the comparison is made)
        assert np.greater(
            np.array(right_chk_with_vsc), np.array(right_chk_no_vsc)
        ).all(), "Lake outflow did no decrease with VSC turned on and should have."


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
