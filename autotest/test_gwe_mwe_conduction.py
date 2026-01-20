# Test mwe package.  Looks at wetted area of well for calculating
# heat conduction exchange.  This test is related to test_gwf_maw06.py but
# with some further customization for testing GWE capabilities.
#  - Test has 0 flow conductance with gw cells in layers 1-3
#  - Test uses MAW to inject water into layer 4 (bottom layer)
#  - Test includes conductive exchng only between MWE feature and
#    layers 1 to 3
#  - Water extracted by normal WEL features in bottom corners of model
#  - Water table saturates only ~1/2 of top layer; therefore
#    conductive exchange

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework


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


def get_welbore_heat_flow(fname, srchStr):
    ener_Q = []
    with open(fname, "r") as f:
        for line in f:
            if srchStr in line:
                # Read an established format
                for i in np.arange(3):  # read & discard 3 lines
                    line = next(f)
                for i in np.arange(4):  # read & digest 4 lines of needed output
                    line = next(f)
                    m_arr = line.strip().split()
                    ener_Q.append([int(m_arr[0]), float(m_arr[2])])
                break

    return np.array(ener_Q)


def trenddetector(list_of_index, array_of_data, order=1):
    result = np.polyfit(list_of_index, list(array_of_data), order)
    slope = result[-2]
    return float(slope)


cases = ["mwe_01"]
mawstrt = 3.5

# Flow related parameters
lx = 70.0
ly = 70.0
nlay = 4
nrow = 7
ncol = 7
nper = 1
delc = ly / nrow
delr = lx / ncol
top = 4.0
botm = [3.0, 2.0, 1.0, 0.0]
strt = 3.5
transient = {0: True}

perlen = [10.0]
nstp = [10]
tsmult = [1.0]

Kh = [1.0, 1.0, 1e-6, 100]
Kv = [1.0, 1.0, 1e-6, 100]

Sy = 0.3
Ss = 0.0

# Transport related parameters
mixelm = 0  # Upstream vs TVD (Upstream selected)
strttemp = 1.0  # Initial temperature ($^{\circ}C$)
porosity = Sy  # porosity (unitless)
ktw = 0.5918  # Thermal Conductivity of Water ($W/m/^{\circ}C$)
kts = 0.2700  # Thermal Conductivity of Aquifer Solids ($W/m/^{\circ}C$)
rhow = 1000  # Density of water ($kg/m^3$)
rhos = 2650  # Density of the aquifer material ($kg/m^3$)
Cpw = 4180  # Heat Capacity of water ($J/kg/C$)
Cps = 880  # Heat capacity of the solids ($J/kg/C$)
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)
K_therm_maw = 1.5  # Thermal conductivity of the lakebed material ($W/m/C$)
wthkcnd = 0.01
mawradius = 0.1
mawbottom = 0.0

# Solver settings
nouter, ninner = 700, 10
hclose, rclose, relax = 1e-8, 1e-6, 0.97


def build_models(idx, test):
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # Instantiate MODFLOW 6 simulation
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # Instantiate Time Discretization package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # Instantiate Groundwater Flow model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions=newtonoptions)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="SIMPLE",
        under_relaxation_gamma=0.1,
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

    # Instantiate Discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        length_units="METERS",
        pname="DIS",
        filename=f"{gwfname}.dis",
    )

    # Instantiate Initial Conditions package
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=strt,
        pname="IC",
        filename=f"{gwfname}.ic",
    )

    # Instantiate Node Property Flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=Kh,
        k33=Kv,
        pname="NPF",
        filename=f"{gwfname}.npf",
    )

    # Instantiate Storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        sy=Sy,
        ss=Ss,
        iconvert=1,
        transient=transient,
        pname="STO",
        filename=f"{gwfname}.sto",
    )

    # Instantiate Well package (for extracting MAW-injected water)
    wellist = []
    for i in [0, nrow - 1]:
        for j in [0, ncol - 1]:
            wellist.append([(nlay - 1, i, j), -7.0, 0.01])

    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wellist,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL",
        auxiliary="TEMPERATURE",
        filename=f"{gwfname}.wel",
    )

    # Instantiate Multi-Aquifer Well package (injects water)
    mstrt = mawstrt
    mawcondeqn = "SPECIFIED"
    mawngwfnodes = nlay
    # <wellno> <radius> <bottom> <strt> <condeqn> <ngwfnodes>
    mawpackagedata = [[0, mawradius, mawbottom, mstrt, mawcondeqn, mawngwfnodes]]
    # <wellno> <icon> <cellid(ncelldim)> <scrn_top> <scrn_bot> <hk_skin> <radius_skin>
    conncond = [0.0, 0.0, 0.0, 1000.0]
    mawconnectiondata = [
        [0, icon, (icon, 3, 3), top, mawbottom, conncond[icon], -999.0]
        for icon in range(nlay)
    ]
    # <wellno> <mawsetting>
    mawperioddata = [[0, "STATUS", "ACTIVE"], [0, "RATE", 28.0]]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        head_filerecord=f"{gwfname}.maw.bin",
        budget_filerecord=f"{gwfname}.maw.cbc",
        packagedata=mawpackagedata,
        connectiondata=mawconnectiondata,
        perioddata=mawperioddata,
        pname="MAW-1",
    )
    opth = f"{gwfname}.maw.obs"
    obsdata = {
        f"{gwfname}.maw.obs.csv": [
            ("whead", "head", (0,)),
        ]
    }
    maw.obs.initialize(filename=opth, digits=20, print_input=True, continuous=obsdata)

    # Instantiate Output Control package
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    # Create GWE model
    # ----------------
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    gwe.name_file.save_flows = True

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

    # Instantiating MODFLOW 6 enregy transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
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

    # Instantiating MODFLOW 6 energy transport initial temperatures
    flopy.mf6.ModflowGweic(gwe, strt=strttemp, filename=f"{gwename}.ic")

    # Instantiate mobile storage and transfer package
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=porosity,
        heat_capacity_water=Cpw,
        density_water=rhow,
        latent_heat_vaporization=lhv,
        heat_capacity_solid=Cps,
        density_solid=rhos,
        pname="MST-1",
        filename=f"{gwename}.mst",
    )

    # Instantiating MODFLOW 6 energy transport advection package
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()

    # Instantiate advection package
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, pname="ADV", filename=f"{gwename}.adv")

    # Instantiate dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=True,
        ktw=0.5918,
        kts=0.2700,
        filename=f"{gwename}.dsp",
    )

    # Instantiate source/sink mixing package
    sourcerecarray = [
        ("WEL", "AUX", "TEMPERATURE"),
    ]
    flopy.mf6.ModflowGwessm(gwe, sources=sourcerecarray, filename=f"{gwename}.ssm")

    # Instantiating MODFLOW 6 transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 17, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwename}.oc",
    )

    # Instantiating MODFLOW 6 multi-well energy transport (mwe) package
    # <wellid>,<strt>, <filter pack cond>, <thickness filter pack>, <boundname>
    mwepackagedata = [(0, 1.0, K_therm_maw, wthkcnd, "well1")]

    mweperioddata = {0: [(0, "RATE", 40.0)]}

    # note: for specifying lake number, use fortran indexing!
    mwe_obs = {
        f"{gwename}.mweobs": [
            ("MweTemp", "temperature", 1),
        ]
    }

    flopy.mf6.ModflowGwemwe(
        gwe,
        flow_package_name="MAW-1",
        budget_filerecord=gwename + ".mwe.bud",
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        packagedata=mwepackagedata,
        mweperioddata=mweperioddata,
        observations=mwe_obs,
        pname="MWE-1",
        filename=f"{gwename}.mwe",
    )

    # Instantiate GWF-GWE exchange
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    top_local = [4.0, 3.0, 2.0, 1.0]
    botm_local = [3.0, 2.0, 1.0, 0.0]

    # calculate volume of water and make sure it is conserved
    name = cases[idx]
    gwfname = "gwf-" + name
    fname = gwfname + ".maw.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.HeadFile(fname, text="HEAD")
    stage = bobj.get_alldata().flatten()

    name = cases[idx]
    gwfname = "gwf-" + name
    fname = gwfname + ".maw.cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    bobj = flopy.utils.CellBudgetFile(fname, precision="double")
    gwfarea = bobj.get_data(text="GWF")

    # Retrieve simulated temperature for the multi-aquifer well
    gwename = "gwe-" + name
    fname = gwename + ".mweobs"
    mwtemp_file = os.path.join(test.workspace, fname)
    assert os.path.isfile(mwtemp_file)
    mwtemp = np.genfromtxt(mwtemp_file, names=True, delimiter=",")
    mwtemp = mwtemp["MWETEMP"].astype(float).reshape((mwtemp.size, 1))

    # Retrieve gw temperatures
    fname = gwename + ".ucn"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    gwtempobj = flopy.utils.HeadFile(fname, precision="double", text="TEMPERATURE")
    gwe_temps = gwtempobj.get_alldata()

    # Calculate conductive exchange external to MF6 and compare to MF6 values
    # Evaluates first time step only
    wellbore_cnd_time1 = []
    for i in np.arange(nlay):
        if stage[0] > top_local[i]:
            thk = top_local[i] - botm_local[i]
        elif stage[0] > botm_local[i]:
            thk = stage[0] - botm_local[i]
        else:
            thk = 0

        # Check that MF6 (GWF) wellbore wetted area matches explicitly calc

        wa = 2 * mawradius * np.pi * thk
        welborecnd = K_therm_maw * wa / wthkcnd
        gw_temp = gwe_temps[0, i, 3, 3]
        deltaT = mwtemp[0][0] - gw_temp

        wellbore_cnd_time1.append(welborecnd * deltaT)

    # Retrieve budget
    fname = os.path.join(test.workspace, gwename + ".lst")
    srchStr = (
        "MWE-1 BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   1"
    )
    T_in, T_out, in_bud_lst, out_bud_lst = get_bud(fname, srchStr)
    assert np.isclose(T_in, T_out, atol=0.1), (
        "There is a heat budget discrepancy where there shouldn't be"
    )

    msg1 = "Conductive heat exchanges calculated explicitly and by MF6 do not match"
    msg2 = (
        "Individually summing well bore 'heat flows' is not matching "
        "the global budget heat flow into the aquifer"
    )
    msg3 = "Groundwater should be warming, but isn't"

    # Get MF6 saved wellbore heat "flows"
    srchStr = "MWE PACKAGE (MWE-1) FLOW RATES   PERIOD      1   STEP        1"
    wbcnd_mf6 = get_welbore_heat_flow(fname, srchStr)

    # Check top 3 layers (4th layer handled different)
    for i in np.arange(nlay - 1):
        assert np.isclose(wbcnd_mf6[i, 1], round(wellbore_cnd_time1[i], 4)), msg1

    # Layer 4 "heat flow" includes convection and conduction, compare
    # "heat flow" from all layers to global budget line item 'IN: GWF'
    glob_bud_gw_in = out_bud_lst["GWF"] + out_bud_lst["WELLBORE-COND"]
    mwe_out = wbcnd_mf6.sum(axis=0)[1]
    assert np.isclose(mwe_out, glob_bud_gw_in, rtol=1e-9), msg2

    # Ensure that temperatures near the injection point are rising
    slp = trenddetector(np.arange(gwe_temps.shape[0]), gwe_temps[:, 3, 3, 3])
    assert slp > 0.0, msg3


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
