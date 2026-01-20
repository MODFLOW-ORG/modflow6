#
#  - Outer columns not active for unsaturated zone, but are present to host
#    constant head boundaries at the bottom of the model.
#
#     +-------+
#     |///////| = Inactive cell
#     +-------+
#
#     Model depiction:
#
#     +-------+-------+-------+
#     |///////|       |///////|  Layer 1
#     +-------+-------+-------+
#     |///////|       |///////|  Layer 2
#     +-------+-------+-------+
#     |///////|       |///////|  Layer 3
#     +-------+-------+-------+
#     |///////|       |///////|
#     + -- -- + -- -- + -- -- +
#     |///////|       |///////|  Layer x  (Middle portion of model not shown)
#     + -- -- + -- -- + -- -- +
#     |///////|       |///////|
#     +-------+-------+-------+
#     |       |       |       |  Layer 99
#     +-------+-------+-------+
#     |       |       |       |  Layer 100
#     +-------+-------+-------+

import math
import os

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pytest
from framework import TestFramework


# Analytical solution, from Barends (2010) Equation 5
def temp_analyt(t, z, t0, tinfil, v, d):
    if t == 0.0:
        temp = t0
    else:
        denom = 2.0 * math.sqrt(d * t)
        ztermm = (z - v * t) / denom
        ztermp = (z + v * t) / denom
        vterm = v * z / d
        if vterm < 100.0:
            # might need to adjust this limit
            temp = t0 + 0.5 * (tinfil - t0) * (
                math.erfc(ztermm) + math.exp(vterm) * math.erfc(ztermp)
            )
        else:
            zeta = 1.0 / (1.0 + 0.47047 * ztermp)
            polyterm = zeta * (0.3480242 + zeta * (-0.0958798 + zeta * 0.7478556))
            temp = t0 + 0.5 * (tinfil - t0) * (
                math.erfc(ztermm) + math.exp(vterm - ztermp**2) * polyterm
            )

    return temp


# Model units
length_units = "meters"
time_units = "days"

nlay, nrow, ncol = 101, 1, 3
nper = 2
perlen = [1.0e9, 100.0]
nstp = [1, 100]
tsmult = len(perlen) * [1.0]

delr = 1.0
delc = 1.0
delz = 0.1  # 10 cm
strt = 0.05
top = 10.0005
botm = [9.9995]  # Top layer is very thin for application of the boundary condition
for i in np.arange(1, nlay):
    bot = 10.0 - (i * delz)
    botm.append(round(bot, 1))

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97
steady = {0: False, 1: False}
transient = {0: True, 1: True}

idomain_u = [0, 1, 0]
idomain_l = [1, 1, 1]
idomain = []
for i in np.arange(nlay):
    if i < 99:
        idomain.append(idomain_u)
    else:
        idomain.append(idomain_l)

idomain = np.array(idomain)

strt_temp = 10.0
scheme = "UPSTREAM"
dispersivity = 0.0
prsity = 0.2
rhow = 1000.0
cpw = 4183.0
lhv = 2500.0
cps = 760.0
rhos = 1500.0

# transient uzf info
# iuzno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
uzf_pkdat = [[0, (0, 0, 1), 1, 1, 0.00001, 1, 0.0001, 0.20, 0.055, 4]]

# Continue building the UZF list of objects
for iuzno in np.arange(1, 101, 1):
    if iuzno < 99:
        ivertconn = iuzno + 1
    else:
        ivertconn = -1

    uzf_pkdat.append(
        [iuzno, (iuzno, 0, 1), 0, ivertconn, 0.01, 1, 0.0001, 0.20, 0.055, 4]
    )

iuz_cell_dict = {}
cell_iuz_dict = {}
for i, itm in enumerate(uzf_pkdat):
    iuz_cell_dict.update({itm[0]: (itm[1][0], itm[1][1], itm[1][2])})
    cell_iuz_dict.update({(itm[1][0], itm[1][1], itm[1][2]): itm[0]})

finf = 0.01
extdp = 0.0
pet = 0.0
extwc = 0.0
zero = 0.0
uzf_spd = {
    0: [[0, finf, pet, extdp, extwc, zero, zero, zero]],
    1: [[0, finf, pet, extdp, extwc, zero, zero, zero]],
}

cases = ["uze00"]


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(sim, time_units=time_units, nper=nper, perioddata=tdis_rc)

    gwfname = "gwf_" + name
    gwename = "gwe_" + name

    newtonoptions = ["NEWTON", "UNDER_RELAXATION"]
    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions=newtonoptions,
        save_flows=True,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=strt,
        filename=f"{gwfname}.ic",
    )

    # node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=100.0,
        k33=10,
        filename=f"{gwfname}.npf",
    )

    # aquifer storage
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-5,
        sy=prsity,
        steady_state=steady,
        transient=transient,
        filename=f"{gwfname}.sto",
    )

    # chd files
    chdval = 0.05
    chdspd = {0: [[(100, 0, 0), chdval, 10.0], [(100, 0, 2), chdval, 10.0]]}

    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary=["TEMPERATURE"],
        print_flows=True,
        stress_period_data=chdspd,
        pname="CHD-1",
        filename=f"{gwfname}.chd",
    )

    # Unsaturated-zone flow package
    flopy.mf6.ModflowGwfuzf(
        gwf,
        print_flows=True,
        save_flows=True,
        wc_filerecord=gwfname + ".uzfwc.bin",
        simulate_et=False,
        simulate_gwseep=False,
        linear_gwet=False,
        boundnames=False,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{gwfname}.uzf.bud",
        pname="UZF-1",
        filename=f"{gwfname}.uzf",
    )

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwfname}.oc",
    )

    # ----------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ----------------------------------
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    gwe.name_file.save_flows = True

    imsgwe = flopy.mf6.ModflowIms(
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
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
        pname="DIS",
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe,
        strt=strt_temp,
        pname="IC",
        filename=f"{gwename}.ic",
    )

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, pname="ADV", filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=False,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        kts=0.2700 * 86400,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 constant temperature boundary condition at
    ctpspd = {0: [[(0, 0, 1), 10.0]], 1: [[(0, 0, 1), 20.0]]}

    flopy.mf6.ModflowGwectp(
        gwe,
        save_flows=True,
        print_flows=True,
        stress_period_data=ctpspd,
        pname="CTP",
        filename=f"{gwename}.ctp",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    srctype = "AUX"
    auxname = "TEMPERATURE"
    pname = ["CHD-1"]
    # Inpput to SSM is: <pname> <srctype> <auxname>
    sources = [[itm, srctype, auxname] for itm in pname]

    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sources,
        pname="SSM",
        filename=f"{gwename}.ssm",
    )

    # Instantiating MODFLOW 6 energy transport source-sink mixing package
    uzepackagedata = [(iuz, 10.0) for iuz in range(nlay)]
    uzeperioddata = {
        0: [[0, "INFILTRATION", 10.0]],
        1: [[0, "INFILTRATION", 20.0]],
    }

    flopy.mf6.ModflowGweuze(
        gwe,
        flow_package_name="UZF-1",
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".uze.bin",
        budget_filerecord=gwename + ".uze.bud",
        packagedata=uzepackagedata,
        uzeperioddata=uzeperioddata,
        pname="UZE-1",
        filename=f"{gwename}.uze",
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        pname="OC",
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
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
    print("evaluating flow...")

    name = cases[idx]
    gwfname = "gwf_" + name
    gwename = "gwe_" + name
    ws = test.workspace

    # check some output...
    wc_fl = gwfname + ".uzfwc.bin"
    wcobj = flopy.utils.HeadFile(os.path.join(ws, wc_fl), text="water-content")
    wc = wcobj.get_alldata()

    fl2 = gwename + ".uze.bin"

    uzeobj = flopy.utils.HeadFile(os.path.join(ws, fl2), text="TEMPERATURE")
    temps = uzeobj.get_alldata()

    t = np.linspace(0.0, 100.0, 101)
    z = np.arange(0.05, 10.0, 0.1)
    z = np.insert(z, 0, 0.0)

    t0 = 10.0
    tinfil = 20.0

    q = finf  # infiltration rate
    rhos = 1500.0
    Cps = 760.0
    rhow = 1000.0
    Cpw = 4183.0
    rhowCpw = Cpw * rhow
    rhosCps = Cps * rhos

    Kts = 23328.0
    Ktw = 0.0

    steady_wc = wc[1, 0, 0, 1]
    Sw = steady_wc / prsity

    rhoCp_bulk = Sw * prsity * rhowCpw + (1 - prsity) * rhosCps
    Kt_bulk = Sw * prsity * Ktw + (1 - prsity) * Kts
    v = rhowCpw / rhoCp_bulk * q
    D = Kt_bulk / rhoCp_bulk

    # Put analytical solution in place
    analytical_sln = np.zeros((len(t), len(z)))
    for i, tm in enumerate(t):
        for j, depth in enumerate(z):
            temp = temp_analyt(tm, depth, t0, tinfil, v, D)
            analytical_sln[i, j] = temp

    # Run checks
    msg0 = (
        "Simulated solution no longer falling within"
        " default tolerance where it previously did"
    )
    # Compare day 1. For layer 20 and below, the defaults of allclose should work
    assert np.allclose(analytical_sln[1, 19:], temps[1, 0, 0, 19:]), msg0
    # Compare day 10. For layer 39 and below, the defaults of allclose should work
    assert np.allclose(analytical_sln[10, 38:], temps[10, 0, 0, 38:]), msg0
    # Compare day 50. For layer 84 and below, the defaults of allclose should work
    assert np.allclose(analytical_sln[50, 83:], temps[50, 0, 0, 83:]), msg0
    # Compare day 100, fits are generally good,
    # but do not pass allclose default settings

    # Ensure that the differences in the 1st day fall within established bounds
    msg1 = (
        "Simulated fits to analytical solution are "
        "falling outside established bounds on day 1"
    )
    assert np.max(analytical_sln[1, :18] - temps[1, 0, 0, :18]) <= 1.52921097880, msg1
    assert np.min(analytical_sln[1, :18] - temps[1, 0, 0, :18]) >= -0.32260871278, msg1

    # Ensure that the differences on day 10 fall within established bounds
    msg2 = (
        "Simulated fits to analytical solution are "
        "falling outside established bounds on day 10"
    )
    assert np.max(analytical_sln[10, :37] - temps[10, 0, 0, :37]) <= 0.15993441016, msg2
    assert np.min(analytical_sln[10, :37] - temps[10, 0, 0, :37]) >= -0.22298707253, (
        msg2
    )

    # Ensure that the differences on day 50 fall within established bounds
    msg3 = (
        "Simulated fits to analytical solution are "
        "falling outside established bounds on day 50"
    )
    assert np.max(analytical_sln[50, :82] - temps[50, 0, 0, :82]) <= 0.09327747258, msg3
    assert np.min(analytical_sln[50, :82] - temps[50, 0, 0, :82]) >= -0.21182907402, (
        msg3
    )

    # Ensure that the differences on day 50 fall within established bounds
    msg3 = (
        "Simulated fits to analytical solution are "
        "falling outside established bounds on day 50"
    )
    assert np.max(analytical_sln[50, :82] - temps[50, 0, 0, :82]) <= 0.09327747258, msg3
    assert np.min(analytical_sln[50, :82] - temps[50, 0, 0, :82]) >= -0.21182907402, (
        msg3
    )

    # Ensure that the differences on day 100 fall within established bounds
    msg4 = (
        "Simulated fits to analytical solution are "
        "falling outside established bounds on day 100"
    )
    assert np.max(analytical_sln[100] - temps[100]) <= 0.107, msg4
    assert np.min(analytical_sln[100] - temps[100]) >= -0.20763221276, msg4

    # If a plot is needed for visual inspection, change following if statement to "True"
    plot_results = False
    if plot_results:
        analytical_sln = np.zeros((len(t), len(z)))
        for i, tm in enumerate(t):
            for j, depth in enumerate(z):
                temp = temp_analyt(tm, depth, t0, tinfil, v, D)
                analytical_sln[i, j] = temp

        # first transient stress period
        line1 = plt.plot(analytical_sln[1], z, "-", color="red", label="Analytical")
        line2 = plt.plot(temps[1, 0, 0], z, "-.", color="blue", label="MODFLOW 6")
        # 10th transient stress period
        plt.plot(analytical_sln[10], z, "-", color="red")
        plt.plot(temps[10, 0, 0], z, "-.", color="blue")
        # 50th transient stress period
        plt.plot(analytical_sln[50], z, "-", color="red")
        plt.plot(temps[50, 0, 0], z, "-.", color="blue")
        # last stress period
        plt.plot(analytical_sln[100], z, "-", color="red")
        plt.plot(temps[100, 0, 0], z, "-.", color="blue")
        # add labels
        plt.text(11.0, 0.85, "1 day", fontsize=10)
        plt.text(12.0, 1.65, "10 days", fontsize=10)
        plt.text(14.0, 2.90, "50 days", fontsize=10)
        plt.text(16.0, 4.00, "100 days", fontsize=10)

        plt.gca().invert_yaxis()
        # For latex replace with: '$Temperature, ^{\circ}C$'
        plt.xlabel("$Temperature, C$")
        plt.ylabel("$Depth, m$")
        plt.minorticks_on()
        plt.axhline(y=0.0)
        plt.legend(loc="lower right", frameon=False)
        plt.savefig(os.path.join(ws, "fit_view.png"), format="png")


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
