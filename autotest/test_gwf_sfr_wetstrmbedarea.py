# Test evap in SFR reaches (no interaction with gwf)

import math
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sfr-wetperim"]


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
Lx = 600.0
Ly = 300.0
nrow = 3
ncol = 6
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
        [101.50, 101.25, 101.00, 100.75, 100.50, 100.25],
        [101.25, 101.00, 100.75, 100.50, 100.25, 100.00],
        [101.50, 101.25, 101.00, 100.75, 100.50, 100.25],
    ]
)

botm = np.zeros(top.shape)
strthd = top - 1.0

# NPF parameters
k11 = 1
ss = 0.00001
sy = 0.20
hani = 1
laytyp = 1

# Package boundary conditions
surf_Q_in = [86400, 0]  # 86400 m^3/d = 1 m^3/s
sfr_evaprate = 0.1
streambed_K = 0.0
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
    if j < 2:
        rise = 1 / 3
        run = 2
        bot_wid = 1.0
    elif j < 4:
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


# time params
steady = {0: True, 1: False}
transient = {0: False, 1: True}
nstp = [1, 1]
tsmult = [1, 1]
perlen = [1, 1]

nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-3, 1e-4, 0.97


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    print(f"Building model...{name}")

    # generate names for each model
    gwfname_trapezoidal = "gwf-" + name

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
        modelname=gwfname_trapezoidal,
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
        filename=f"{gwfname_trapezoidal}.ims",
    )
    sim.register_ims_package(ims, [gwfname_trapezoidal])

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
        budget_filerecord=f"{gwfname_trapezoidal}.cbc",
        head_filerecord=f"{gwfname_trapezoidal}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate streamflow routing package
    # Determine the middle row and store in rMid (account for 0-base)
    rMid = 1
    # sfr data
    nreaches = ncol
    rlen = delr
    roughness = 0.035
    rbth = 1.0
    rhk = streambed_K
    strm_up = 100.25
    strm_dn = 99
    slope = (strm_up - strm_dn) / ((ncol - 1) * delr)
    ustrf = 1.0
    ndv = 0
    strm_incision = 1.0

    # use trapezoidal cross-section for channel geometry
    sfr_xsec_tab_nm1 = f"{gwfname_trapezoidal}.xsec.tab1"
    sfr_xsec_tab_nm2 = f"{gwfname_trapezoidal}.xsec.tab2"
    sfr_xsec_tab_nm3 = f"{gwfname_trapezoidal}.xsec.tab3"
    sfr_xsec_tab_nm = [sfr_xsec_tab_nm1, sfr_xsec_tab_nm2, sfr_xsec_tab_nm3]
    crosssections = []
    for n in range(nreaches):
        # 6 reaches, 3 cross section types
        crosssections.append([n, sfr_xsec_tab_nm[n // 2]])

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
            rwid[irch // 2],
            slope,
            top[rMid, irch] - strm_incision,
            rbth,
            rhk,
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
            sfrbndx.append([i, "EVAPORATION", sfr_evaprate])

        sfr_perioddata.update({t: sfrbndx})

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname_trapezoidal}.sfr.obs.csv": [
            ("rch1_depth", "depth", 1),
            ("rch2_depth", "depth", 2),
            ("rch3_depth", "depth", 3),
            ("rch4_depth", "depth", 4),
            ("rch5_depth", "depth", 5),
            ("rch6_depth", "depth", 6),
        ],
        "digits": 8,
        "print_input": True,
        "filename": name + ".sfr.obs",
    }

    budpth = f"{gwfname_trapezoidal}.sfr.cbc"
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
        filename=f"{gwfname_trapezoidal}.sfr",
    )

    return sim, None


def check_output(idx, test):
    # read flow results from model
    name = cases[idx]
    gwfname = "gwf-" + name

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
    for j, stg in enumerate(list(sfrstg[0])[1:]):
        wp = calc_wp(j, stg)
        wa = wp * delr
        msg = (
            "Wetted streambed area for reach "
            + str(j)
            + "in stress period 1 does not match explicitly-calculated answer"
        )
        assert np.isclose(wa, shared_area[0, j], atol=1e-4), msg

    msg = "Wetted streambed area of all reaches should be zero in stress period 2"
    for val in list(sfrstg[1])[1:]:
        assert val == 0.0, msg


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
