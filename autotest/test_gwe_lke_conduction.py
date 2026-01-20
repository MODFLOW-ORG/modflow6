# Simple single lake model.  Lake cut into top two layers of a 5 layer
# model.  Model is loosely based on the first example problem in
# Merritt and Konikow (2000) which also is one of the MT3D-USGS test
# problems.  This test developed to isolate lake-aquifer interaction;
# no SFR or other advanced packages.  Problem set up to have only
# conductive exchange with groundwater, then groundwater pass-through:
# that is, gw inflow on the left side, gw outflow on the
# right side of the lake.
#
# starting groundwater temperature: 4.0
# left chd boundary inflow temperature: 4.0
# starting lake temperature: 20.0
#

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


def trenddetector(list_of_index, array_of_data, order=1):
    result = np.polyfit(list_of_index, list(array_of_data), order)
    slope = result[-2]
    return float(slope.item())


cases = ["lke-conductn", "lke-conductm", "lke-conductr"]

#
# The last letter in the names above indicates the following
# n = "no lk/gw exchange"
# m = "mixed" (i.e., convection one direction, conductive gradient the other direction?)
# r = "reversed thermal gradients (warm gw, cold lake)

strt_gw_temp = [4.0, 4.0, 20.0]
strt_lk_temp = [20.0, 20.0, 4.0]
lak_leakance = [0.0, 1.0, 1.0]
strt_lk_stg = [33.75, 33.75, 33.75]
lkbdthkcnd = [0.0001, 0.0001, 0.0001]  # Thickness to consider for feature/gw conduction

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
nstp = [10]
tsmult = [1.0]
perlen = [5000]

# solver params
nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-3, 1e-4, 0.97

# Transport related parameters
al = 1  # longitudinal dispersivity ($m$)
ath1 = al  # horizontal transverse dispersivity
atv = al  # vertical transverse dispersivity
mixelm = 0  # Upstream vs TVD (Upstream selected)
porosity = 0.20  # porosity (unitless)
K_therm = 2.0  # Thermal conductivity  # ($W/m/C$)
rhow = 1000  # Density of water ($kg/m^3$)
rhos = 2650  # Density of the aquifer material ($kg/m^3$)
Cpw = 4180  # Heat Capacity of water ($J/kg/C$)
Cps = 880  # Heat capacity of the solids ($J/kg/C$)
lhv = 2454000.0  # Latent heat of vaporization ($J/kg$)
K_therm_lakebed = 1.5  # Thermal conductivity of the lakebed material ($W/m/C$)


#
# MODFLOW 6 flopy GWF & GWE simulation object (sim) is returned
#


def build_models(idx, test):
    global lak_lkup_dict

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
                chdlistl.append([(k, i, 0), chdl, strt_gw_temp[idx]])
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
                                ilak,  # <lakeno>
                                ilakconn,  # <iconn>
                                (k, i - 1, j),  # <cellid(ncelldim)>
                                "horizontal",  # <claktype>
                                lak_leakance[idx],  # <bedleak>
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
                                lak_leakance[idx],
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
                                lak_leakance[idx],
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
                                lak_leakance[idx],
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
                        lak_leakance[idx],
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                    ]
                    lakeconnectiondata.append(v)
                    lak_lkup_dict.update({ilakconn: (k, i, j)})

    strtStg = strt_lk_stg[idx]
    lakpackagedata = [[0, strtStg, len(lakeconnectiondata), strt_lk_temp[idx], "lake1"]]
    lak_pkdat_dict = {"filename": "lak_pakdata.in", "data": lakpackagedata}

    lakeperioddata = {0: []}

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

    # -----------------
    # Create GWE model
    # -----------------

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
        idomain=ibound,
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 energy transport initial temperatures
    strttemp = strt_gw_temp[idx]
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
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme, filename=f"{gwename}.adv")

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
        ("CHD-L", "AUX", "TEMPERATURE"),
        ("CHD-R", "AUX", "TEMPERATURE"),
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

    # Instantiating MODFLOW 6 lake energy transport (lke) package
    lkepackagedata = [(0, strt_lk_temp[idx], K_therm_lakebed, lkbdthkcnd[idx], "lake1")]

    # lkeperioddata = {0: [(0, "STATUS", "CONSTANT"), (0, "TEMPERATURE", 4.0)]}

    # note: for specifying lake number, use fortran indexing!
    lke_obs = {
        f"{gwename}.lakobs": [
            ("resTemp", "temperature", 1),
            ("resGwEnerExchng", "lke", "lake1"),
        ]
    }

    flopy.mf6.ModflowGwelke(
        gwe,  # Set time_conversion for use with Manning's eqn.
        flow_package_name="LAK-1",
        flow_package_auxiliary_name="TEMPERATURE",
        budget_filerecord=gwename + ".lke.bud",
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=False,
        print_temperature=True,
        packagedata=lkepackagedata,
        # lakeperioddata=lkeperioddata,
        observations=lke_obs,
        pname="LKE-1",
        filename=f"{gwename}.lke",
    )

    # GWF-GWE exchange
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

    # read flow results from model
    name = cases[idx]
    gwename = "gwe-" + name

    # Retrieve simulated temperature for the lake
    fname = gwename + ".lakobs"
    lktemp_file = os.path.join(test.workspace, fname)
    lktemp = np.genfromtxt(lktemp_file, names=True, delimiter=",")
    lktemp = lktemp["RESTEMP"].astype(float).reshape((lktemp.size, 1))

    # Retrieve groundwater temperatures
    fname = gwename + ".ucn"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    gwtempobj = flopy.utils.HeadFile(fname, precision="double", text="TEMPERATURE")
    gwe_temps = gwtempobj.get_alldata()

    # gw exchng (item 'GWF') should be zero in heat transport budget
    srchStr = (
        "LKE-1 BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   1"
    )
    fname = gwename + ".lst"
    fname = os.path.join(test.workspace, fname)

    # Retrieve budget
    T_in, T_out, in_bud_lst, out_bud_lst = get_bud(fname, srchStr)
    assert np.isclose(T_in, T_out, atol=0.1), (
        "There is a heat budget discrepancy where there shouldn't be"
    )

    msg1 = "Budget item 'GWF' should be 0.0 for this scenario"
    msg2 = (
        "Thermal conduction is occurring in the wrong direction based "
        "on the thermal gradient between the lake and groundwater system"
    )
    msg3 = (
        "There should be a cooling trend in the lake based on heat loss "
        "to the groundwater system"
    )
    msg4 = "There should be a warming trend in the groundwater adjacent to the lake"
    msg5 = (
        "Budget item 'GWF' should reflect heat entering the lake (via gw/sw exchange)"
    )
    msg6 = "Budget item 'GWF' should reflect heat exiting the lake (via gw/sw exchange)"

    if name[-1] == "n":
        assert in_bud_lst["GWF"] == 0.0, msg1
        assert out_bud_lst["GWF"] == 0.0, msg1

    if name[-1] != "n":
        assert in_bud_lst["GWF"] > 0.0, msg5
        assert out_bud_lst["GWF"] > 0.0, msg6

    # Determine gw/sfe temperature gradient direction
    if lktemp[0] > gwe_temps[0, 0, 0, 0]:
        # conduction will be from lake to gw cells
        assert in_bud_lst["LAKEBED-COND"] == 0.0, msg2
        assert out_bud_lst["LAKEBED-COND"] > 0.0, msg2

        slp = trenddetector(np.arange(len(lktemp)), lktemp)
        # Lake should be cooling through conductive exchange with cold gw
        assert slp < 0.0, msg3

        slp = trenddetector(np.arange(lktemp.shape[0]), gwe_temps[:, 1, 8, 11])
        # gw should be warming through conductive exchange with a warm lake
        assert slp > 0.0, msg4

    else:  # thermally reversed scenario (cold lake, warm gw)
        # conduction will be from gw cells to lake
        assert in_bud_lst["LAKEBED-COND"] > 0.0, msg2
        assert out_bud_lst["LAKEBED-COND"] == 0.0, msg2

        slp = trenddetector(np.arange(len(lktemp)), lktemp)
        # Lake should be warming through conductive exchange with warm gw
        assert slp > 0.0, msg3

        slp = trenddetector(np.arange(lktemp.shape[0]), gwe_temps[:, 1, 8, 11])
        # gw should be cooling through conductive exchange with a cold lake
        assert slp < 0.0, msg4


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
