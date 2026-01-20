"""
Test evap in SFR reaches (no interaction with gwf)
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sfr-evap"]


# Model units
length_units = "m"
time_units = "days"

# model domain and grid definition
Lx = 500.0
Ly = 300.0
nrow = 3
ncol = 5
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
        [101.25, 101.00, 100.75, 100.50, 100.25],
        [101.00, 100.75, 100.50, 100.25, 100.00],
        [101.25, 101.00, 100.75, 100.50, 100.25],
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
surf_Q_in = 86400  # 1 m^3/s
sfr_evaprate = 0.1
streambed_K = 0.0

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
    gwfname_trapezoidal = "gwf-" + name + "-t"
    gwfname_rectangular = "gwf-" + name + "-r"

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

    gwft = flopy.mf6.ModflowGwf(
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
        gwft,
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
        gwft,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
    )

    # Instantiate storage package
    flopy.mf6.ModflowGwfsto(
        gwft,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state=steady,
        transient=transient,
    )

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwft, strt=strthd)

    # Instantiate output control package
    flopy.mf6.ModflowGwfoc(
        gwft,
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
    rwid = 9.0
    roughness = 0.035
    rbth = 1.0
    rhk = streambed_K
    strm_up = 100
    strm_dn = 99
    slope = (strm_up - strm_dn) / ((ncol - 1) * delr)
    ustrf = 1.0
    ndv = 0
    strm_incision = 1.0

    # use trapezoidal cross-section for channel geometry
    x_coord = [0.0, 2.0, 4.0, 5.0, 7.0, 9.0]
    x_xsec = [val / rwid for val in x_coord]
    y_xsec = [0.66666667, 0.33333333, 0.0, 0.0, 0.33333333, 0.66666667]
    x_sec_tab = [[x, h] for x, h in zip(x_xsec, y_xsec)]

    sfr_xsec_table_name = f"{gwfname_trapezoidal}.xsec.tab"
    crosssections = []
    for n in range(nreaches):
        crosssections.append([n, sfr_xsec_table_name])

    flopy.mf6.ModflowUtlsfrtab(
        gwft,
        nrow=len(x_xsec),
        ncol=2,
        table=x_sec_tab,
        filename=sfr_xsec_table_name,
        pname="sfrxsectable",
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
            rwid,
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

    sfrbndx = []
    for i in np.arange(nreaches):
        if i == 0:
            sfrbndx.append([i, "INFLOW", surf_Q_in])
        sfrbndx.append([i, "EVAPORATION", sfr_evaprate])

    sfr_perioddata = {0: sfrbndx}

    # Instantiate SFR observation points
    sfr_obs = {
        f"{gwfname_trapezoidal}.sfrobs": [
            ("rch1_in", "ext-inflow", 1),  # For now, these need to be 1-based
            ("rch1_rain", "rainfall", 1),
            ("rch1_evap", "evaporation", 1),
            ("rch1_runoff", "runoff", 1),
            ("rch1_flow", "sfr", 1),
            ("rch1_out", "downstream-flow", 1),
            ("rch1_xsec_area", "wet-area", 1),
            ("rch1_top_width", "wet-width", 1),
            ("rch1_depth", "depth", 1),
            ("rch1_stg", "stage", 1),
        ],
        "digits": 8,
        "print_input": True,
        "filename": name + ".sfr.obs",
    }

    budpth = f"{gwfname_trapezoidal}.sfr.cbc"
    flopy.mf6.ModflowGwfsfr(
        gwft,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400.0,
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

    # ---------------------------------------------------------------------------
    # Add a second GWF model that simulates a typical rectangular stream channel
    # ---------------------------------------------------------------------------

    gwfr = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname_rectangular,
        save_flows=True,
        newtonoptions="newton",
    )

    # Instantiating solver
    # (use same settings as above)
    ims2 = flopy.mf6.ModflowIms(
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
    sim.register_ims_package(ims2, [gwfname_rectangular])

    # Instantiate discretization package
    flopy.mf6.ModflowGwfdis(
        gwfr,
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
        gwfr,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=k11,
    )

    # Instantiate storage package
    flopy.mf6.ModflowGwfsto(
        gwfr,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state=steady,
        transient=transient,
    )

    # Instantiate initial conditions package
    flopy.mf6.ModflowGwfic(gwfr, strt=strthd)

    # Instantiate output control package
    flopy.mf6.ModflowGwfoc(
        gwfr,
        budget_filerecord=f"{gwfname_rectangular}.cbc",
        head_filerecord=f"{gwfname_rectangular}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate streamflow routing package
    # Uses same datasets established above except does not include x-sections
    budpth = f"{gwfname_rectangular}.sfr.cbc"
    flopy.mf6.ModflowGwfsfr(
        gwfr,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        budget_filerecord=budpth,
        mover=False,
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=sfr_perioddata,
        pname="SFR-2",
        filename=f"{gwfname_rectangular}.sfr",
    )

    return sim, None


def check_output(idx, test):
    # read flow results from model
    name = cases[idx]
    gwfname_t = "gwf-" + name + "-t"
    gwfname_r = "gwf-" + name + "-r"

    fname = gwfname_t + ".sfr.cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    sfrobj = flopy.utils.binaryfile.CellBudgetFile(fname, precision="double")
    sfrevap = sfrobj.get_data(text="EVAPORATION")

    # Extract evap
    sim_evap = []
    for i in range(ncol):
        sim_evap.append(sfrevap[-1][i][2])

    sim_evap = np.array(sim_evap)

    # Establish known answer:
    stored_strm_evap = np.array(
        [-62.17272623, -62.15731943, -62.14191043, -62.12649925, -62.11108587]
    )

    msg = "The SFR evaporation test with n-point x-section (trapezoid) is failing."
    assert np.allclose(stored_strm_evap, sim_evap, atol=1e-4), msg

    # Now check results from standard rectangular x-section setup
    # (not an n-point channel)
    fname2 = gwfname_r + ".sfr.cbc"
    fname2 = os.path.join(test.workspace, fname2)
    assert os.path.isfile(fname2)

    sfrobj = flopy.utils.binaryfile.CellBudgetFile(fname2, precision="double")
    sfrevap_r = sfrobj.get_data(text="EVAPORATION")

    # Extract evap
    sim_evap_r = []
    for i in range(ncol):
        sim_evap_r.append(sfrevap_r[-1][i][2])

    sim_evap_r = np.array(sim_evap_r)

    # Establish known answer:
    stored_strm_evap_r = np.array([-90.0] * 5)

    msg = "The SFR evaporation test with rectangular geometry is failing."
    assert np.allclose(stored_strm_evap_r, sim_evap_r, atol=1e-4), msg


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
