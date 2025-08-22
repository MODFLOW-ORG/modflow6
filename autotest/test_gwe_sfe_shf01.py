# Test the use of the sensible heat flux utility used in conjunction with the
# SFE advanced package.  This test checks to make sure that simulated sensible
# heat flux amounts are inline with what would be expected given changes in
# specific input parameters (while holding the others constant).  Relative
# changes in stress periods 2 through 4 are compared back to stress period 1
# or among the 3 reaches within each stress period

# This test to include:
#  - typical ways of specifying input
#  - hot gw cell warming an upstream reach
#  - thermally hot stream water warming host gw cells
#


import os

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["sfe-shf", "sfe-shf-ts"]
#
# The last letter in the names above indicates the following
# n = "no gw/sw exchange"
# i = "gwf into strm"
# o = "strm to gw"
# m = "mixed" (i.e., convection one direction, conductive gradient the other direction?)


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

# Boundary conditions
strt_gw_temp = [4.0, 4.0]
chd_condition = ["n", "i", "o", "m"]

# NPF parameters
ss = 0.00001
sy = 0.20
hani = 1
laytyp = 1
k11 = 500.0
# SFR/SFE
rhk = [0.0, k11]
strm_temp = [18.0, 18.0]
rlen = delr
surf_Q_in = [8.64, 86.4, 8.64, 8.64]  # 86400 m^3/d = 1 m^3/s = 35.315 cfs
# SHF
# Remember, from the first to the second stress period, the flow rate increases
# Stress periods 3 and 4 return to the same flow rate as stress period 1
# For stress period 3, increase wpd (everything else remains as is)
# For stress period 4, increase tatm (everything else remains as is)
wspd = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [4.0, 5.0, 6.0], [1.0, 1.0, 1.0]]
tatm = [[10.0, 10.0, 10.0], [10.0, 10.0, 10.0], [10.0, 10.0, 10.0], [15.0, 20.0, 30.0]]


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
rbthcnd = 0.0001

# time params
steady = {0: False, 1: False}
transient = {0: True, 1: True}
nstp = [1, 1, 1, 1]
tsmult = [1, 1, 1, 1]
perlen = [1, 1, 1, 1]

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
    chdelev1 = top[0, 0] - 0.05
    chdelev2 = top[0, -1] - 0.05
    gw_temp = strt_gw_temp[idx]
    if chd_on:
        chdlist1 = [
            [(0, 0, 0), chdelev1, gw_temp],
            [(0, nrow - 1, 0), chdelev1, gw_temp],
            [(0, 0, ncol - 1), chdelev2, gw_temp],
            [(0, nrow - 1, ncol - 1), chdelev2, gw_temp],
        ]
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist1,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="CHD",
            auxiliary="TEMPERATURE",
            filename=f"{gwfname}.chd",
        )

    # Instantiate streamflow routing package
    # Determine the middle row and store in rMid (account for 0-base)
    rMid = 1
    # sfr data
    nreaches = ncol
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
    for t in np.arange(len(surf_Q_in)):
        sfrbndx = []
        for i in np.arange(nreaches):
            if i == 0:
                sfrbndx.append([i, "INFLOW", surf_Q_in[t]])
            # sfrbndx.append([i, "EVAPORATION", sfr_evaprate])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname}.sfr.obs.csv": [
            ("rch1_width", "wet-width", 1),
            ("rch2_width", "wet-width", 2),
            ("rch3_width", "wet-width", 3),
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
        pname="SFR",
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
    sourcerecarray = [("CHD", "AUX", "TEMPERATURE")]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{gwename}.ssm")

    # Instantiate Streamflow Energy Transport package
    sfepackagedata = []
    for irno in range(ncol):
        t = (irno, strm_temp[idx], K_therm_strmbed[irno], rbthcnd)
        sfepackagedata.append(t)

    sfeperioddata = []
    for irno in range(ncol):
        if irno == 0:
            sfeperioddata.append((irno, "INFLOW", strm_temp[idx]))

    # Instantiate SFE observation points
    sfe_obs = {
        f"{gwename}.sfe.obs.csv": [
            ("rch1_outftemp", "temperature", 1),
            ("rch1_outfener", "ext-outflow", 1),
            ("rch1_shf", "shf", 1),
            ("rch2_outftemp", "temperature", 2),
            ("rch2_outfener", "ext-outflow", 2),
            ("rch2_shf", "shf", 2),
            ("rch3_outftemp", "temperature", 3),
            ("rch3_outfener", "ext-outflow", 3),
            ("rch3_shf", "shf", 3),
        ],
        "digits": 8,
        "print_input": True,
        "filename": gwename + ".sfe.obs",
    }

    shf_filename = f"{gwename}.sfe.shf"
    sfe = flopy.mf6.modflow.ModflowGwesfe(
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
        observations=sfe_obs,
        flow_package_name="SFR",
        pname="SFE",
        filename=f"{gwename}.sfe",
    )

    # shf
    shf_spd = {}
    for kper in range(len(nstp)):
        spd = []
        for irno in range(ncol):
            spd.append([irno, "WSPD", wspd[kper][irno]])
            spd.append([irno, "TATM", tatm[kper][irno]])
        shf_spd[kper] = spd

    shf = flopy.mf6.ModflowUtlshf(
        sfe,
        print_input=True,
        density_air=1.225,
        heat_capacity_air=717.0,
        drag_coefficient=0.002,
        reachperioddata=shf_spd,
        filename=shf_filename,
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

    # Retrieve simulated top width of each reach
    sfr_pth0 = os.path.join(test.workspace, f"{gwfname}.sfr.obs.csv")
    assert os.path.isfile(sfr_pth0)
    sfroutdf = pd.read_csv(sfr_pth0)

    # assert that each successive reach grows in surface area
    msg0 = "surface area not increasing from reach 1 to reach 2 as expected."
    msg1 = "surface area not increasing from reach 2 to reach 3 as expected."
    for i in np.arange(sfroutdf.shape[0]):
        assert sfroutdf.loc[i, "RCH1_WIDTH"] < sfroutdf.loc[i, "RCH2_WIDTH"], msg0
        assert sfroutdf.loc[i, "RCH2_WIDTH"] < sfroutdf.loc[i, "RCH3_WIDTH"], msg1

    # Retrieve select SFE output
    sfe_pth0 = os.path.join(test.workspace, f"{gwename}.sfe.obs.csv")
    assert os.path.isfile(sfe_pth0), "missing sfe observation output file"
    sfeoutdf = pd.read_csv(sfe_pth0)

    # assert that each successive reach grows in sensible heat flux
    # based on the problem setup, each stress period should reflect an increase
    # in sensible heat flux relative to the first stress period for different
    # reasons.
    # From SP1 to SP2: inflow increase, therefore surface area increases
    #                  (trapezoidal cross section)
    # From SP1 to SP3: wind speed increases in each successive reach
    # From SP1 to SP4: temperature of the atmosphere increases each successive
    #                  reach
    for i in np.arange(sfeoutdf.shape[0]):
        assert abs(sfeoutdf.loc[i, "RCH1_SHF"]) < abs(sfeoutdf.loc[i, "RCH2_SHF"]), (
            "magnitude of shf not increasing from reach 1 to reach 2 as expected."
        )
        assert abs(sfeoutdf.loc[i, "RCH2_SHF"]) < abs(sfeoutdf.loc[i, "RCH3_SHF"]), (
            "magnitude of shf not increasing from reach 2 to reach 3 as expected."
        )

    # as a result of the amount of energy leaving from shf, ensure that
    # temperatures are changing accordingly in each successive reach
    msg2 = (
        "temperatures should be decreasing in the downstream direction as a "
        "result of sensible heat flux losses"
    )
    for i in np.arange(sfeoutdf.shape[0]):
        assert sfeoutdf.loc[i, "RCH1_OUTFTEMP"] > sfeoutdf.loc[i, "RCH2_OUTFTEMP"], msg2
        assert sfeoutdf.loc[i, "RCH2_OUTFTEMP"] > sfeoutdf.loc[i, "RCH3_OUTFTEMP"], msg2

    # there should be no "external" energy outflow for the first or second
    # reaches
    msg3 = (
        "should not be external energy outflow from reach 1 in any of the "
        "stress periods"
    )
    msg4 = (
        "should not be external energy outflow from reach 2 in any of the "
        "stress periods"
    )
    assert sfeoutdf.loc[:, "RCH1_OUTFENER"].sum() == 0, msg3
    assert sfeoutdf.loc[:, "RCH2_OUTFENER"].sum() == 0, msg4

    # there should be non-zero energy exported from the final reach in every
    # stress period
    msg5 = "export of energy from last reach should be negative"
    for i in np.arange(sfeoutdf.shape[0]):
        assert sfeoutdf.loc[i, "RCH3_OUTFENER"] < 0, msg5


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
