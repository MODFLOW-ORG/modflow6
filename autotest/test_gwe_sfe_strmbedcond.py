# Test conduction between an advanced package feature, in this case stream
# reaches with varying channel geometries and the host GWE gw cells.
# This test should include:
#  - with gw-sw interaction
#  - hot gw cell warming an upstream reach
#  - thermally hot stream water warming host gw cells
#

import math
import os

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["sfe-conductn", "sfe-conducti", "sfe-conducto", "sfe-conductm", "sfe_conductk"]
#
# The last letter in the names above indicates the following
# n = "no gw/sw exchange"
# i = "gwf into strm"
# o = "strm to gw"
# m = "mixed" (i.e., convection one direction, conductive gradient the other direction?)
# k = "known answer"

k11 = 500.0
rhk = [0.0, k11, k11, k11, 0.0]
strt_gw_temp = [4.0, 4.0, 4.0, 20.0, 1.0]
strm_temp = [18.0, 18.0, 20.0, 4.0, 10.0]
chd_condition = ["n", "i", "o", "m", "k"]
surf_Q_in = [
    [8.64, 0.0],
    [8640.0, 0.0],
    [8640.0, 0.0],
    [8640.0, 0.0],
    [8640.0, 0.0],
]  # 86400 m^3/d = 1 m^3/s = 35.315 cfs


def get_x_frac(x_coord1, rwid):
    x_xsec1 = [val / rwid for val in x_coord1]
    return x_xsec1


def get_xy_pts(x, y, rwid):
    x_xsec1 = get_x_frac(x, rwid)
    x_sec_tab = [[xx, hh] for xx, hh in zip(x_xsec1, y)]
    return x_sec_tab


# Model units
length_units = "m"
time_units = "days"

# model domain and grid definition
Lx = 90.0
Ly = 90.0
nrow = 3
ncol = 3
nlay = 1
delr = Lx / ncol
delc = Ly / nrow
xmax = ncol * delr
ymax = nrow * delc
X, Y = np.meshgrid(
    np.linspace(delr / 2, xmax - delr / 2, ncol),
    np.linspace(ymax - delc / 2, 0 + delc / 2, nrow),
)
ibound = np.ones((nlay, nrow, ncol))
# Because eqn uses negative values in the Y direction, need to do a little manipulation
Y_m = -1 * np.flipud(Y)
top = np.array(
    [
        [101.50, 101.25, 101.00],
        [101.25, 101.00, 100.75],
        [101.50, 101.25, 101.00],
    ]
)

botm = np.array(
    [
        [98.5, 98.25, 98.0],
        [98.25, 98.0, 97.75],
        [98.5, 98.25, 98.0],
    ]
)
strthd = 98.75
chd_on = True

# NPF parameters
ss = 0.00001
sy = 0.20
hani = 1
laytyp = 1

# Package boundary conditions
sfr_evaprate = 0.1
rwid = [9.0, 10.0, 20]
# Channel geometry: trapezoidal
x_sec_tab1 = get_xy_pts(
    [0.0, 2.0, 4.0, 5.0, 7.0, 9.0],
    [0.66666667, 0.33333333, 0.0, 0.0, 0.33333333, 0.66666667],
    rwid[0],
)

x_sec_tab2 = get_xy_pts(
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    [0.5, 0.25, 0.0, 0.0, 0.25, 0.5],
    rwid[1],
)

x_sec_tab3 = get_xy_pts(
    [0.0, 4.0, 8.0, 12.0, 16.0, 20.0],
    [0.33333333, 0.16666667, 0.0, 0.0, 0.16666667, 0.33333333],
    rwid[2],
)
x_sec_tab = [x_sec_tab1, x_sec_tab2, x_sec_tab3]


def calc_wp(j, stg):
    if j < 1:
        rise = 1 / 3
        run = 2
        bot_wid = 1.0
    elif j < 2:
        rise = 1 / 4
        run = 2
        bot_wid = 2.0
    else:
        rise = 1 / 6
        run = 4
        bot_wid = 4.0

    ang = math.atan2(rise, run)
    hyp_len = stg / math.sin(ang)
    wp = hyp_len * 2 + bot_wid

    return wp


def process_line(line):
    m_arr = line.strip().split()
    if any("=" in itm and len(itm) > 1 for itm in m_arr):
        m_arr = [
            float(itm.split("=")[-1]) if len(itm.split("=")) > 1 else itm
            for itm in m_arr
        ]
        nm = m_arr[-2]
    else:
        nm = m_arr[-3]
    val = m_arr[-1]
    return {nm: float(val)}


def get_bud(fname, srchStr):
    in_bud_lst = {}
    out_bud_lst = {}
    with open(fname, "r") as f:
        for line in f:
            if srchStr in line:
                # Read the package budget
                line = next(f)
                while "TOTAL IN =" not in line:
                    if "=" in line:
                        in_bud_lst.update(process_line(line))

                    line = next(f)

                # Get "total in"
                dct = process_line(line)
                T_in = dct["IN"]

                line = next(f)
                while "TOTAL OUT =" not in line:
                    if "=" in line:
                        out_bud_lst.update(process_line(line))

                    line = next(f)

                # Get "total out"
                dct = process_line(line)
                T_out = dct["OUT"]

                break

    return T_in, T_out, in_bud_lst, out_bud_lst


def trenddetector(list_of_index, array_of_data, order=1):
    result = np.polyfit(list_of_index, list(array_of_data), order)
    slope = result[-2]
    return float(slope)


# Transport related parameters
porosity = sy  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rhow = 1000  # Density of water ($kg/m^3$)
rhos = 2650  # Density of the aquifer material ($kg/m^3$)
Cpw = 4180  # Heat capacity of water ($J/kg/C$)
Cps = 880  # Heat capacity of the solids ($J/kg/C$)
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)
# Thermal conductivity of the streambed material ($W/m/C$)
K_therm_strmbed = [1.5, 1.75, 2.0]
rbthcnd = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001]

# time params
steady = {0: False, 1: False}
transient = {0: True, 1: True}
nstp = [1, 1]
tsmult = [1, 1]
perlen = [1, 1]

nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-3, 1e-4, 0.97

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
    gwename = "gwe-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating time discretization
    tdis_rc = []
    for i in range(len(nstp)):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(
        sim, nper=len(nstp), perioddata=tdis_rc, time_units=time_units
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
    )

    # Instantiate node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
    )

    # Instantiate storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state=steady,
        transient=transient,
    )

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=strthd)

    # Instantiate output control package
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate constant head boundary package
    if chd_condition[idx] == "n":
        chdelev1 = top[0, 0] - 3.0
        chdelev2 = top[0, -1] - 3.0
    elif chd_condition[idx] == "i":
        chdelev1 = top[0, 0] - 0.05
        chdelev2 = top[0, -1] - 0.05
    elif chd_condition[idx] == "o":
        chdelev1 = top[0, 0] - 3.0
        chdelev2 = top[0, -1] - 3.0
    elif chd_condition[idx] == "m":
        chdelev1 = top[0, 0] - 3.0  # convection from stream to gw,
        chdelev2 = top[0, -1] - 3.0  # conduction from gw to strm
    elif chd_condition[idx] == "k":
        chdelev1 = top[0, 0] - 3.0
        chdelev2 = top[0, -1] - 3.0

    # Instantiate constant head boundary package
    if chd_on:
        chdlist1 = [
            [(0, 0, 0), chdelev1, strt_gw_temp[idx]],
            [(0, 0, ncol - 1), chdelev2, strt_gw_temp[idx]],
            [(0, nrow - 1, 0), chdelev1, strt_gw_temp[idx]],
            [(0, nrow - 1, ncol - 1), chdelev2, strt_gw_temp[idx]],
        ]
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist1,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="CHD-1",
            auxiliary="TEMPERATURE",
            filename=f"{gwfname}.chd",
        )

    # Instantiate streamflow routing package
    # Determine the middle row and store in rMid (account for 0-base)
    rMid = 1
    # sfr data
    nreaches = ncol
    rlen = delr
    roughness = 0.035
    rbth = 1.0
    strmbd_hk = rhk[idx]
    strm_up = 100.25
    strm_dn = 99
    # divide by 10 to further reduce slop
    slope = (strm_up - strm_dn) / ((ncol - 1) * delr) / 10
    ustrf = 1.0
    ndv = 0
    strm_incision = 1.0

    # use trapezoidal cross-section for channel geometry
    sfr_xsec_tab_nm1 = f"{gwfname}.xsec.tab1"
    sfr_xsec_tab_nm2 = f"{gwfname}.xsec.tab2"
    sfr_xsec_tab_nm3 = f"{gwfname}.xsec.tab3"
    sfr_xsec_tab_nm = [sfr_xsec_tab_nm1, sfr_xsec_tab_nm2, sfr_xsec_tab_nm3]
    crosssections = []
    for n in range(nreaches):
        # 3 reaches, 3 cross section types
        crosssections.append([n, sfr_xsec_tab_nm[n]])

    # Setup the tables
    for n in range(len(x_sec_tab)):
        flopy.mf6.ModflowUtlsfrtab(
            gwf,
            nrow=len(x_sec_tab[n]),
            ncol=2,
            table=x_sec_tab[n],
            filename=sfr_xsec_tab_nm[n],
            pname="sfrxsectable" + str(n + 1),
        )

    packagedata = []
    for irch in range(nreaches):
        nconn = 1
        if 0 < irch < nreaches - 1:
            nconn += 1
        rp = [
            irch,
            (0, rMid, irch),
            rlen,
            rwid[irch],
            slope,
            top[rMid, irch] - strm_incision,
            rbth,
            strmbd_hk,
            roughness,
            nconn,
            ustrf,
            ndv,
        ]
        packagedata.append(rp)

    connectiondata = []
    for irch in range(nreaches):
        rc = [irch]
        if irch > 0:
            rc.append(irch - 1)
        if irch < nreaches - 1:
            rc.append(-(irch + 1))
        connectiondata.append(rc)

    sfr_perioddata = {}
    for t in np.arange(len(surf_Q_in[idx])):
        sfrbndx = []
        for i in np.arange(nreaches):
            if i == 0:
                sfrbndx.append([i, "INFLOW", surf_Q_in[idx][t]])
            # sfrbndx.append([i, "EVAPORATION", sfr_evaprate])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname}.sfr.obs.csv": [
            ("rch1_depth", "depth", 1),
            ("rch2_depth", "depth", 2),
            ("rch3_depth", "depth", 3),
        ],
        "digits": 8,
        "print_input": True,
        "filename": name + ".sfr.obs",
    }

    budpth = f"{gwfname}.sfr.cbc"
    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400,
        budget_filerecord=budpth,
        mover=False,
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        crosssections=crosssections,
        perioddata=sfr_perioddata,
        observations=sfr_obs,
        pname="SFR-1",
        filename=f"{gwfname}.sfr",
    )

    # --------------------------------------------------
    # Setup the GWE model for simulating heat transport
    # --------------------------------------------------
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename)

    # Instantiating solver for GWT
    imsgwe = flopy.mf6.ModflowIms(
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
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwename])

    # Instantiating DIS for GWE
    flopy.mf6.ModflowGwedis(
        gwe,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        pname="DIS",
        filename=f"{gwename}.dis",
    )

    # Instantiate Mobile Storage and Transfer package
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=porosity,
        heat_capacity_water=Cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=Cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # Instantiate Energy Transport Initial Conditions package
    flopy.mf6.ModflowGweic(gwe, strt=strt_gw_temp[idx])

    # Instantiate Advection package
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")

    # Instantiate Dispersion package (also handles conduction)
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        ktw=0.5918,
        kts=0.2700,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    # [b/c at least one boundary back is active (SFR), ssm must be on]
    sourcerecarray = [("CHD-1", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{gwename}.ssm")

    # Instantiate Streamflow Energy Transport package
    sfepackagedata = []
    for irno in range(ncol):
        t = (irno, strm_temp[idx], K_therm_strmbed[irno], rbthcnd[idx])
        sfepackagedata.append(t)

    sfeperioddata = []
    for irno in range(ncol):
        if irno == 0:
            sfeperioddata.append((irno, "INFLOW", strm_temp[idx]))
        # sfeperioddata.append((irno, sfr_applied_bnd[idx], sfe_applied_temp[idx]))

    sfe_obs = {
        (gwename + ".sfe.obs.csv",): [
            (f"sfe-{i + 1}-temp", "TEMPERATURE", i + 1) for i in range(3)
        ]
        + [
            ("sfe-extin", "EXT-INFLOW", 1),
            ("sfe-rain", "RAINFALL", 1),
            ("sfe-roff", "RUNOFF", 1),
            ("sfe-evap", "EVAPORATION", 1),
            ("sfe-extout", "EXT-OUTFLOW", 3),
            ("sfe-sfe", "SFE", 2),
            ("sfe-strmbd1", "STRMBD-COND", 1),
            ("sfe-strmbd2", "STRMBD-COND", 2),
            ("sfe-strmbd3", "STRMBD-COND", 3),
        ],
    }

    flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=False,
        save_flows=True,
        print_input=False,
        print_flows=False,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=sfepackagedata,
        reachperioddata=sfeperioddata,
        flow_package_name="SFR-1",
        observations=sfe_obs,
        pname="SFE-1",
        filename=f"{gwename}.sfe",
    )

    # Instantiate Output Control package for transport
    flopy.mf6.ModflowGweoc(
        gwe,
        temperature_filerecord=f"{gwename}.ucn",
        saverecord=[("TEMPERATURE", "ALL")],
        temperatureprintrecord=[("COLUMNS", 3, "WIDTH", 20, "DIGITS", 8, "GENERAL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwename}.oc",
    )

    # Instantiate Gwf-Gwe Exchange package
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{gwename}.gwfgwe",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    fname = gwfname + ".sfr.cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    sfrobj = flopy.utils.binaryfile.CellBudgetFile(fname, precision="double")
    sfr_wetted_interface_area = sfrobj.get_data(text="gwf")

    # Retrieve simulated stage of each reach
    sfr_pth0 = os.path.join(test.workspace, f"{gwfname}.sfr.obs.csv")
    sfrstg = np.genfromtxt(sfr_pth0, names=True, delimiter=",")

    # Extract shared wetted interfacial areas
    shared_area = []
    for t in range(len(sfr_wetted_interface_area)):
        sp_area = []
        for i in range(ncol):
            sp_area.append(sfr_wetted_interface_area[t][i][3])

        shared_area.append(sp_area)

    shared_area = np.array(shared_area)

    # Calculate wetted streambed area for comparison
    wa_lst = []
    for j, stg in enumerate(list(sfrstg[0])[1:]):
        wp = calc_wp(j, stg)
        wa = wp * delr
        msg = (
            "Wetted streambed area for reach "
            + str(j)
            + "in stress period 1 does not match explicitly-calculated answer"
        )

        assert np.isclose(wa, shared_area[0, j], atol=1e-4), msg
        wa_lst.append(wa)

    # Sub-scenario checks
    # initialize search term
    srchStr = (
        "SFE-1 BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   1"
    )
    fname = "gwe-" + name + ".lst"
    fname = os.path.join(test.workspace, fname)

    # pull in SFE CSV output for comparison with lst file budget term
    fpth = os.path.join(test.workspace, gwename + ".sfe.obs.csv")
    df = pd.read_csv(fpth)

    # gw exchng (item 'GWF') should be zero in heat transport budget
    T_in, T_out, in_bud_lst, out_bud_lst = get_bud(fname, srchStr)
    assert np.isclose(T_in, T_out, atol=0.1), "There is a heat budget discrepancy"

    # compare individual streambed conduction obs with total
    df["sum_strmbd_cond"] = df.iloc[:, -3:].sum(axis=1)
    if name[-1] != "m":
        assert np.isclose(
            out_bud_lst["STRMBD-COND"],
            abs(df.loc[0, "sum_strmbd_cond"]),
            atol=0.0001,
        ), "There is a streambed conductance discrepancy " + str(
            out_bud_lst["STRMBD-COND"] - abs(df.loc[0, "sum_strmbd_cond"])
        )
    else:
        assert np.isclose(
            in_bud_lst["STRMBD-COND"], abs(df.loc[0, "sum_strmbd_cond"]), atol=0.0001
        ), "There is a streambed conductance discrepancy " + str(
            out_bud_lst["STRMBD-COND"] - abs(df.loc[0, "sum_strmbd_cond"])
        )

    # Get temperature of streamwater
    fname1 = "gwe-" + name + ".sfe.bin"
    fname1 = os.path.join(test.workspace, fname1)
    sfeobj = flopy.utils.HeadFile(fname1, precision="double", text="TEMPERATURE")
    sfe_temps = sfeobj.get_alldata()

    # Get temperature of gw
    fname2 = "gwe-" + name + ".ucn"
    fname2 = os.path.join(test.workspace, fname2)
    gwobj = flopy.utils.HeadFile(fname2, precision="double", text="TEMPERATURE")
    gw_temps = gwobj.get_alldata()

    msg1 = "Budget item 'GWF' should be 0.0 for this scenario"
    msg2 = "Thermal conduction is occurring in the wrong direction"
    msg3 = (
        "There should be a decreasing temperatures trend in "
        "downstream direction owing to conductive losses"
    )
    msg4 = (
        "There should be an increasing temperature trend in "
        "the row of cells hosting the stream owing to increasing "
        "conductive losses from the stream to the aquifer "
        "(i.e., greater shared wetted areas)"
    )
    msg5 = (
        "The conductive exchange of energy calculated by GWE in the 5th "
        "sub-test does not match an externally calculated solution"
    )
    if name[-1] == "n":
        # no gw/sw convective exchange, simulates conductive exchange only
        assert in_bud_lst["GWF"] == 0.0, msg1
        assert out_bud_lst["GWF"] == 0.0, msg1

        # Determine gw/sfe temperature gradient direction
        if sfe_temps[0, 0, 0, 0] > gw_temps[0, 0, 0, 0]:
            # conduction will be from stream to gw
            assert in_bud_lst["STRMBD-COND"] == 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] > 0.0, msg2

            slp = trenddetector(
                np.arange(0, sfe_temps.shape[-1]), sfe_temps[0, 0, 0, :]
            )
            assert slp < 0.0, msg3

            slp = trenddetector(np.arange(0, gw_temps.shape[-2]), gw_temps[0, 0, 1, :])
            assert slp > 0.0, msg4

        else:
            assert in_bud_lst["STRMBD-COND"] > 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] == 0.0, msg2

    # streamflow gain from aquifer ("into stream")
    if name[-1] == "i":
        msg = "Budget item 'GWF' should reflect heat entering stream"
        assert in_bud_lst["GWF"] > 0.0, msg
        assert out_bud_lst["GWF"] == 0.0, msg

        # Determine gw/sfe temperature gradient direction
        if sfe_temps[0, 0, 0, 0] > gw_temps[0, 0, 0, 0]:
            # conduction will be from stream to gw
            assert in_bud_lst["STRMBD-COND"] == 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] > 0.0, msg2

            slp = trenddetector(
                np.arange(0, sfe_temps.shape[-1]), sfe_temps[0, 0, 0, :]
            )
            assert slp < 0.0, msg3

            slp = trenddetector(np.arange(0, gw_temps.shape[-2]), gw_temps[0, 0, 1, :])
            assert slp > 0.0, msg4

        else:
            assert in_bud_lst["STRMBD-COND"] > 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] == 0.0, msg2

    # streamflow loss to aquifer ("out of stream")
    if name[-1] == "o":
        msg = "Budget item 'GWF' should reflect heat exiting stream"
        assert in_bud_lst["GWF"] == 0.0, msg
        assert out_bud_lst["GWF"] > 0.0, msg

        # Determine gw/sfe temperature gradient direction
        if sfe_temps[0, 0, 0, 0] > gw_temps[0, 0, 0, 0]:
            # conduction will be from stream to gw
            assert in_bud_lst["STRMBD-COND"] == 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] > 0.0, msg2

            slp = trenddetector(
                np.arange(0, sfe_temps.shape[-1]), sfe_temps[0, 0, 0, :]
            )
            assert slp < 0.0, msg3

            slp = trenddetector(np.arange(0, gw_temps.shape[-2]), gw_temps[0, 0, 1, :])
            assert slp < 0.0, msg4

        else:
            assert in_bud_lst["STRMBD-COND"] > 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] == 0.0, msg2

    # Reverse temperature gradient  (cold stream, warm aquifer)
    # Loss of streamwater to aquifer
    # Thus, convection from strm to gw, conduction from gw to strm
    if name[-1] == "m":  # 'm' for mixed
        msg = "Budget item 'GWF' should reflect heat exiting stream"
        assert in_bud_lst["GWF"] == 0.0, msg
        assert out_bud_lst["GWF"] > 0.0, msg

        # Determine gw/sfe temperature gradient direction
        if sfe_temps[0, 0, 0, 0] > gw_temps[0, 0, 0, 0]:
            # conduction will be from stream to gw
            assert in_bud_lst["STRMBD-COND"] == 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] > 0.0, msg2

        else:
            assert in_bud_lst["STRMBD-COND"] > 0.0, msg2
            assert out_bud_lst["STRMBD-COND"] == 0.0, msg2

            slp = trenddetector(
                np.arange(0, sfe_temps.shape[-1]), sfe_temps[0, 0, 0, :]
            )
            assert slp > 0.0, msg3

            slp = trenddetector(np.arange(0, gw_temps.shape[-2]), gw_temps[0, 0, 1, :])
            assert slp > 0.0, msg4

    if name[-1] == "k":  # 'k' for known
        wa = shared_area[0][0]
        K_t_sb = K_therm_strmbed[0]
        sbthermthk = rbthcnd[idx]

        # final stream temperature, reach 1
        strm_temp = df.at[0, "SFE-1-TEMP"]
        gw_temp = gw_temps[0, 0, 1, 0]

        thermcond = wa * K_t_sb / sbthermthk * (gw_temp - strm_temp)

        assert np.isclose(thermcond, df.at[0, "SFE-STRMBD1"], atol=1e-4), (
            msg5
            + ". Values are thermcond: "
            + str(thermcond)
            + "   GWE: "
            + str(df.at[0, "SFE-STRMBD1"])
            + "   wa: "
            + str(wa)
            + "   ktsb: "
            + str(K_t_sb)
            + "   sbthk: "
            + str(sbthermthk)
        )


# - No need to change any code below
@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
