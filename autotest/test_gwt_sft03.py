# Test for evaluating change in storage in transient GWT model.
# The check is for SFT stream concentrations with and without the use of
# the STORAGE keyword in the options block of SFR. When the STORAGE keyword
# is used, the SFT-calculated concentrations are correct in the first
# stress period. When omitted, SFT concentrations are wrong in the first
# stress period. Checks use both rectangular and trapazoidal channel
# geometry. First stress period is transient.

import os

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["sft_rect", "sft_rectfail", "sft_trap", "sft_trapfail"]


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
strt_gw_conc = 4.0

# NPF parameters
ss = 0.00001
sy = 0.20
hani = 1
laytyp = 1
k11 = 500.0
# SFR/SFT
rhk = 0.0
strm_conc = 18.0
rlen = delr
surf_Q_in = 8.64  # 86400 m^3/d = 1 m^3/s = 35.315 cfs


# Package boundary conditions
rwid_rect = [9.0, 9.0, 9.0]
rwid_trap = [9.0, 10.0, 20]
# Channel geometry: trapezoidal
x_sec_tab1 = get_xy_pts(
    [0.0, 2.0, 4.0, 5.0, 7.0, 9.0],
    [0.66666667, 0.33333333, 0.0, 0.0, 0.33333333, 0.66666667],
    rwid_trap[0],
)

x_sec_tab2 = get_xy_pts(
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    [0.5, 0.25, 0.0, 0.0, 0.25, 0.5],
    rwid_trap[1],
)

x_sec_tab3 = get_xy_pts(
    [0.0, 4.0, 8.0, 12.0, 16.0, 20.0],
    [0.33333333, 0.16666667, 0.0, 0.0, 0.16666667, 0.33333333],
    rwid_trap[2],
)
x_sec_tab = [x_sec_tab1, x_sec_tab2, x_sec_tab3]


# Transport related parameters
porosity = sy  # porosity (unitless)

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
    gwtname = "gwt-" + name

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
    gw_conc = strt_gw_conc
    if chd_on:
        chdlist1 = [
            [(0, 0, 0), chdelev1, gw_conc],
            [(0, nrow - 1, 0), chdelev1, gw_conc],
            [(0, 0, ncol - 1), chdelev2, gw_conc],
            [(0, nrow - 1, ncol - 1), chdelev2, gw_conc],
        ]
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist1,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="CHD",
            auxiliary="CONCENTRATION",
            filename=f"{gwfname}.chd",
        )

    # Instantiate streamflow routing package
    # Determine the middle row and store in rMid (account for 0-base)
    rMid = 1
    # sfr data
    nreaches = ncol
    roughness = 0.035
    rbth = 1.0
    strmbd_hk = rhk
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

    init_stgs = []
    if "rect" in name:
        rwid = rwid_rect
        init_stgs = [
            [0, 100.250908380705],
            [1, 100.000908380705],
            [2, 99.7509083807047],
        ]
    elif "trap" in name:
        rwid = rwid_trap
        init_stgs = [
            [0, 100.2533702860543],
            [1, 100.0022366903203],
            [2, 99.75148429736241],
        ]

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
    for t in np.arange(2):
        sfrbndx = []
        for i in np.arange(nreaches):
            if i == 0:
                sfrbndx.append([i, "INFLOW", surf_Q_in])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname}.sfr.obs.csv": [
            ("rch1_width", "wet-width", 1),
            ("rch2_width", "wet-width", 2),
            ("rch3_width", "wet-width", 3),
        ],
        "digits": 20,
        "print_input": True,
        "filename": name + ".sfr.obs",
    }

    budpth = f"{gwfname}.sfr.cbc"
    if "rect" in name:
        # omit 'crosssections=' argument since using rectangular geometry
        if "fail" not in name:
            flopy.mf6.ModflowGwfsfr(
                gwf,
                storage=True,
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
                initialstages=init_stgs,
                perioddata=sfr_perioddata,
                observations=sfr_obs,
                pname="SFR",
                filename=f"{gwfname}.sfr",
            )
        else:
            # by omitting storage, the initial sft storage calculation will be off
            flopy.mf6.ModflowGwfsfr(
                gwf,
                # storage=True,
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
                # initialstages=init_stgs,
                perioddata=sfr_perioddata,
                observations=sfr_obs,
                pname="SFR",
                filename=f"{gwfname}.sfr",
            )

    elif "trap" in name:
        if "fail" not in name:
            flopy.mf6.ModflowGwfsfr(
                gwf,
                storage=True,
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
                initialstages=init_stgs,
                perioddata=sfr_perioddata,
                observations=sfr_obs,
                pname="SFR",
                filename=f"{gwfname}.sfr",
            )
        else:
            flopy.mf6.ModflowGwfsfr(
                gwf,
                # storage=True,
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
                # initialstages=init_stgs,
                perioddata=sfr_perioddata,
                observations=sfr_obs,
                pname="SFR",
                filename=f"{gwfname}.sfr",
            )

    # ----------------------------------------------------
    # Setup the GWT model for simulating solute transport
    # ----------------------------------------------------
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    # Instantiating solver for GWT
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
    sim.register_ims_package(imsgwt, [gwtname])

    # Instantiating DIS for gwt
    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        pname="DIS",
        filename=f"{gwtname}.dis",
    )

    # Instantiate Mobile Storage and Transfer package
    flopy.mf6.ModflowGwtmst(
        gwt,
        save_flows=True,
        porosity=porosity,
        pname="MST",
        filename=f"{gwtname}.mst",
    )

    # Instantiate Mass Transport Initial Conditions package
    flopy.mf6.ModflowGwtic(gwt, strt=strt_gw_conc)

    # Instantiate Advection package
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # Instantiate Dispersion package (also handles conduction)
    flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        pname="DSP",
        filename=f"{gwtname}.dsp",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    # [b/c at least one boundary back is active (SFR), ssm must be on]
    sourcerecarray = [("CHD", "AUX", "CONCENTRATION")]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm")

    # Instantiate Streamflow Mass Transport package
    sftpackagedata = []
    for irno in range(ncol):
        t = (irno, strm_conc)
        sftpackagedata.append(t)

    sftperioddata = []
    for irno in range(ncol):
        if irno == 0:
            sftperioddata.append((irno, "INFLOW", strm_conc))

    # Instantiate SFT observation points
    sft_obs = {
        f"{gwtname}.sft.obs.csv": [
            ("rch1_outfconc", "concentration", 1),
            ("rch2_outfconc", "concentration", 2),
            ("rch3_outfconc", "concentration", 3),
            ("rch1_outfmass", "ext-outflow", 1),
            ("rch2_outfmass", "ext-outflow", 2),
            ("rch3_outfmass", "ext-outflow", 3),
        ],
        "digits": 8,
        "print_input": True,
        "filename": gwtname + ".sft.obs",
    }

    sft = flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=False,
        save_flows=True,
        print_input=False,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sftpackagedata,
        reachperioddata=sftperioddata,
        observations=sft_obs,
        flow_package_name="SFR",
        pname="SFT",
        filename=f"{gwtname}.sft",
    )

    # Instantiate Output Control package for transport
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("concentration", "ALL")],
        concentrationprintrecord=[("COLUMNS", 3, "WIDTH", 20, "DIGITS", 8, "GENERAL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwtname}.oc",
    )

    # Instantiate Gwf-gwt Exchange package
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{gwtname}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name

    # Retrieve simulated top width of each reach
    sft_pth0 = os.path.join(test.workspace, f"{gwtname}.sft.obs.csv")
    assert os.path.isfile(sft_pth0)
    sftoutdf = pd.read_csv(sft_pth0)

    # assert that each successive reach grows in surface area
    msg0 = "Concentrations in every reach for every SP should be " + str(strm_conc)
    msg1 = "Because STORAGE not specified in SFR OPTIONS block, \
            concentrations should not be as expected."
    if "fail" not in name:
        for i in np.arange(ncol):
            assert np.all(
                sftoutdf.loc[:, "RCH" + str(i + 1) + "_OUTFCONC"].to_numpy() == 18.0
            ), msg0
    else:
        for i in np.arange(ncol):
            assert np.all(
                sftoutdf.loc[:, "RCH" + str(i + 1) + "_OUTFCONC"].to_numpy() != 18.0
            ), msg1


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
