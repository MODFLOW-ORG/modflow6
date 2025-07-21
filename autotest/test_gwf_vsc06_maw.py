"""
This test originally submitted as a bug report in issue #2323 on GitHub
(https://github.com/MODFLOW-ORG/modflow6/issues/2323).  Test sets up two
models, one with WEL and one with MAW with both using VSC to make sure
the MAW results roughly match the WEL results.  When heads from the two
models are roughly equal, MAW is interacting with VSC OK.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["vsc06"]

# time units
time_units = "DAYS"

# model data
h_above = 42  # head for layer 0
h_ufa = 47  # head for aquifer (layers 1-6 - zero based)
h_below = 48  # head for layer 7 - zero based
tds_above = 0.013  # TDS in layer 0
tds_ufa = 0.024  # TDS in layers 1-6
tds_below = 0.028  # TDS in layer 7
temp_above = 24.0  # Temperature in layer 0
temp_ufa = 25.0  # Temperature in layer 1-6
temp_below = 27.0  # Temperature in layer 7
Nlay = 8  # numjber of layers
N = 101  # number of rows/columns
L = 4000.0  # length of grid in x and y
k_above = 0.0007  # k for layer 0
k_fz = 300  # k for layers 1 (this is the flow zone)
k_ufa = 75  # k for layers 2 - 6 (rest of the aquifer)
k_below = 0.2  # k for layer 7
q_total = 668403  # flow for well in ft3/d
bot = np.array(
    [-560, -585, -625, -665, -705, -745, -785, -905], dtype=float
)  # layer bottoms
top = -425  # top of layer 0
ss = 0.00001
sy = 0.05
porosity = 0.30

# Calculated values to be used in the packages.
hvals = np.array(
    [h_above, h_ufa, h_ufa, h_ufa, h_ufa, h_ufa, h_ufa, h_below], dtype=float
)
start = hvals[:, np.newaxis, np.newaxis] * np.ones((Nlay, N, N))
cvals1 = np.array(
    [tds_above, tds_ufa, tds_ufa, tds_ufa, tds_ufa, tds_ufa, tds_ufa, tds_below],
    dtype=float,
)
start1 = cvals1[:, np.newaxis, np.newaxis] * np.ones((Nlay, N, N))
cvals2 = np.array(
    [
        temp_above,
        temp_ufa,
        temp_ufa,
        temp_ufa,
        temp_ufa,
        temp_ufa,
        temp_ufa,
        temp_below,
    ],
    dtype=float,
)
start2 = cvals2[:, np.newaxis, np.newaxis] * np.ones((Nlay, N, N))
k = np.array([k_above, k_fz, k_ufa, k_ufa, k_ufa, k_ufa, k_ufa, k_below], dtype=float)
k33 = np.array(
    [
        k_above / 2,
        k_fz / 10,
        k_ufa / 10,
        k_ufa / 10,
        k_ufa / 10,
        k_ufa / 10,
        k_ufa / 10,
        k_below / 2,
    ],
    dtype=float,
)
# TDS is normalized such that salt water (35000 mg/L) is set at 1;
# 0 mg/L is 0 in this model
Csalt = 1.0
Cfresh = 0.0
densesalt = 62.4  # using density units of lb/ft3
densefresh = 62.25
denseslp_tds = (densesalt - densefresh) / (Csalt - Cfresh)

# Temperature is in degrees C
Ccold = 10.0
Chot = 45.0
densecold = 62.4062
densehot = 61.8818
denseslp_temp = (densehot - densecold) / (Chot - Ccold)

# MAW well is placed in the same location across the same layers,
# with the same flow rates
maw_packagedata = [[0, 0.1, bot[6], h_ufa, "THIEM", 6, 0.0, 24.5]]
maw_connectiondata = [
    [0, 0, (1, int((N - 1) / 2), int((N - 1) / 2)), bot[0], bot[1], 1e10, 0.2],
    [0, 1, (2, int((N - 1) / 2), int((N - 1) / 2)), bot[1], bot[2], 1e10, 0.2],
    [0, 2, (3, int((N - 1) / 2), int((N - 1) / 2)), bot[2], bot[3], 1e10, 0.2],
    [0, 3, (4, int((N - 1) / 2), int((N - 1) / 2)), bot[3], bot[4], 1e10, 0.2],
    [0, 4, (5, int((N - 1) / 2), int((N - 1) / 2)), bot[4], bot[5], 1e10, 0.2],
    [0, 5, (6, int((N - 1) / 2), int((N - 1) / 2)), bot[5], bot[6], 1e10, 0.2],
]
maw_perioddata = {}
maw_perioddata[0] = [0, "Rate", 0.0]
maw_perioddata[1] = [0, "Rate", q_total]
maw_perioddata[2] = [0, "Rate", 0.0]
maw_perioddata[3] = [0, "Rate", -q_total]


# Setup constant head
chd_rec = []
for layer in range(1, Nlay - 1):
    h_layer = hvals[layer]
    c1_layer = cvals1[layer]
    c2_layer = cvals2[layer]
    for row_col in range(0, N):
        chd_rec.append(((layer, row_col, 0), h_layer, c1_layer, c2_layer))
        chd_rec.append(((layer, row_col, N - 1), h_layer, c1_layer, c2_layer))
        if row_col != 0 and row_col != N - 1:
            chd_rec.append(((layer, 0, row_col), h_layer, c1_layer, c2_layer))
            chd_rec.append(((layer, N - 1, row_col), h_layer, c1_layer, c2_layer))

for row in range(0, N):
    for col in range(0, N):
        chd_rec.append(((0, row, col), hvals[0], cvals1[0], cvals2[0]))
        chd_rec.append(
            ((Nlay - 1, row, col), hvals[Nlay - 1], cvals1[Nlay - 1], cvals2[Nlay - 1])
        )


# compute transmissivities and the percent of flow that should go through
# each layer from the well
thickness = (
    bot[:-1] - bot[1:]
)  # compute thickness of all layers except 0 (the well doesn't pull from layer 0)
k_vals = np.full_like(thickness, k_ufa)  # for each layer, assign the UFA k
k_vals[0] = k_fz  # layer 1, the first in this array is the flow zone, change the k
transmiss = k_vals * thickness  # compute transmissivities
transmiss = transmiss[
    :-1
]  # we don't need the last transmissivity because the well doesn't extend to layer 7
sum_transmiss = np.sum(transmiss)  # sum them up
flowpercent = (
    transmiss / sum_transmiss
)  # compute percentage of flow for each layer of the aquifer


# Add a trio of models - 1 GWF & 2 GWT models, the second of which uses the old
# 'parameter equivalents' approach to simulate temperature
def add_models(sim, gwfname, gwtname, gwtname2, ct=0, wel_on=False, maw_on=False):
    # Instantiate flow model
    model_nam_file = gwfname + ".nam"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=gwfname + ".nam")

    ims = flopy.mf6.modflow.mfims.ModflowIms(
        sim,
        pname="ims-" + str(ct),
        complexity="SIMPLE",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])

    # set up TDS model with solver
    gwt = flopy.mf6.ModflowGwt(
        sim,
        modelname=gwtname,
        model_nam_file=gwtname + ".nam",
    )
    ims2 = flopy.mf6.modflow.mfims.ModflowIms(
        sim,
        filename=gwtname + ".ims",
        linear_acceleration="BICGSTAB",
        pname="ims-" + str(ct + 1),
    )
    sim.register_ims_package(ims2, [gwtname])

    # Set up Temperature model with solver
    gwt2 = flopy.mf6.ModflowGwt(
        sim, modelname=gwtname2, model_nam_file=gwtname2 + ".nam"
    )
    ims3 = flopy.mf6.modflow.mfims.ModflowIms(
        sim,
        filename=gwtname2 + ".ims",
        linear_acceleration="BICGSTAB",
        pname="ims-" + str(ct + 2),
    )
    sim.register_ims_package(ims3, [gwtname2])

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        pname="sto",
        ss=ss,
        sy=sy,
        steady_state={0: True},
        transient={1: True},
    )

    # Set up dis packages
    delrow = delcol = L / (N - 1)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        nlay=Nlay,
        nrow=N,
        ncol=N,
        delr=delrow,
        delc=delcol,
        top=top,
        botm=bot,
        pname="dis-" + str(ct),
    )
    dis2 = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwt,
        nlay=Nlay,
        nrow=N,
        ncol=N,
        delr=delrow,
        delc=delcol,
        top=top,
        botm=bot,
        pname="dis-" + str(ct + 1),
    )
    dis3 = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwt2,
        nlay=Nlay,
        nrow=N,
        ncol=N,
        delr=delrow,
        delc=delcol,
        top=top,
        botm=bot,
        pname="dis-" + str(ct + 2),
    )

    # Set up starting heads
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, strt=start, pname="ic-" + str(ct))
    # Set up starting TDS
    ic2 = flopy.mf6.modflow.mfgwtic.ModflowGwtic(
        gwt,
        strt=start1,
        pname="ic-" + str(ct + 1),
    )
    # Set up starting Temperature (this uses the old "parameter equivalents" approach)
    ic3 = flopy.mf6.modflow.mfgwtic.ModflowGwtic(
        gwt2,
        pname="ic-" + str(ct + 2),
        strt=start2,
    )

    # Set up npf package
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=k,
        save_flows=True,
        save_specific_discharge=True,
        pname="npf-" + str(ct),
    )

    # Set up CHD package
    # CHD is applied to all of layer 0 and layer 7
    # CHD is applied to the edges of layers 1-6
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        auxiliary=["TDS", "TEMP"],
        maxbound=len(chd_rec),
        stress_period_data=chd_rec,
        save_flows=True,
        pname="chd-" + str(ct),
    )

    # Setup input for WEL package
    # A single well in the center of the model extending from layer 1 to 6
    # no pumping in first stress period
    # well injects during second stress period
    # no pumping in third stress period
    # well extracts during fourth stress period

    if wel_on:
        period_two = flopy.mf6.ModflowGwfwel.stress_period_data.empty(
            gwf, aux_vars=["TDS", "TEMP"], maxbound=6
        )
        # period_three = flopy.mf6.ModflowGwfwel.stress_period_data.empty(
        #     gwf, aux_vars=['TDS','TEMP'], maxbound=6
        # )
        # period_four = flopy.mf6.ModflowGwfwel.stress_period_data.empty(
        #     gwf, aux_vars=['TDS','TEMP'], maxbound=6
        # )
        for i in range(6):
            # flows are based on transmissivity ratios
            period_two[0][i] = (
                (i + 1, int((N - 1) / 2), int((N - 1) / 2)),
                q_total * flowpercent[i],
                0.0,
                24.5,
            )
            # period_three[0][i] = (
            #     (i+1, int((N-1)/2),
            #     int((N-1)/2)),
            #     0,
            #     0.0,
            #     24.5,
            # )
            # period_four[0][i] = (
            #     (i+1, int((N-1)/2), int((N-1)/2)),
            #     -q_total*flowpercent[i],
            #     0.0,
            #     24.5
            # )

        stress_period_data = {}
        stress_period_data[1] = period_two[0]
        # stress_period_data[2] = period_three[0]
        # stress_period_data[3] = period_four[0]

        wel = flopy.mf6.ModflowGwfwel(
            gwf,
            auxiliary=[("TDS", "TEMP")],
            maxbound=6,
            stress_period_data=stress_period_data,
            save_flows=True,
            pname="wel",
        )

    if maw_on:
        maw = flopy.mf6.ModflowGwfmaw(
            gwf,
            save_flows=True,
            print_flows=True,
            print_input=True,
            print_head=True,
            nmawwells=1,
            auxiliary=["TDS", "TEMP"],
            packagedata=maw_packagedata,
            connectiondata=maw_connectiondata,
            perioddata=maw_perioddata,
            pname="maw",
        )

    # Instantiate the BUY package
    buypd = [
        (0, denseslp_tds, 0.0, gwtname, "TDS"),
        (1, denseslp_temp, 25.0, gwtname2, "TEMP"),
    ]
    buy = flopy.mf6.ModflowGwfbuy(
        gwf,
        denseref=62.25,
        nrhospecies=2,
        packagedata=buypd,
    )

    # Add viscosity package.  Units are ft, days for flow and lbf/ft3 for
    # density, so viscosity units are lbf-d/ft2. An adjustment to the
    # TDS/viscosity slope is made because TDS is normalized to 0 for
    # freshwater and 1 for saltwater (35000 mg/L).
    viscpd = [
        (0, 1.69e-11, 0.0, gwtname, "TDS"),
        (1, 5.79e-12, 25.0, gwtname2, "TEMP"),
    ]
    vsc = flopy.mf6.ModflowGwfvsc(
        gwf,
        viscref=2.105e-10,
        temperature_species_name="TEMP",
        thermal_formulation="NONLINEAR",
        nviscspecies=2,
        packagedata=viscpd,
        viscosity_filerecord="visc_" + str(ct) + ".out",
        pname="vsc-" + str(ct),
    )

    # GWF output control
    headfile = gwfname + ".hds"
    head_filerecord = [headfile]
    budgetfile = gwfname + ".cbb"
    budget_filerecord = [budgetfile]
    saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
    printrecord = [("HEAD", "LAST")]
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        saverecord=saverecord,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        printrecord=printrecord,
        pname="oc-" + str(ct),
    )

    # -----------------------------
    # GWT - solute transport model
    # -----------------------------

    # Setup solute and temperature transport models
    adv1 = flopy.mf6.ModflowGwtadv(
        gwt,
        scheme="UPSTREAM",
        pname="adv-" + str(ct + 1),
    )
    mst1 = flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        pname="mst-" + str(ct + 1),
    )
    dsp1 = flopy.mf6.ModflowGwtdsp(
        gwt,
        alh=2.5,
        atv=2.5,
        ath1=2.5,
        diffc=0.0,
        pname="dsp-" + str(ct + 1),
    )
    # ssm package
    sourcerecarray_gwt = []
    sourcerecarray_gwte = []
    if wel_on:
        sourcerecarray_gwt = [
            ["chd-" + str(ct), "AUX", "TDS"],
            ["wel", "AUX", "TDS"],
        ]
        sourcerecarray_gwte = [
            ["chd-" + str(ct), "AUX", "TEMP"],
            ["wel", "AUX", "TEMP"],
        ]
    elif maw_on:
        sourcerecarray_gwt = [
            ["chd-" + str(ct), "AUX", "TDS"],
            ["maw", "AUX", "TDS"],
        ]
        sourcerecarray_gwte = [
            ["chd-" + str(ct), "AUX", "TEMP"],
            ["maw", "AUX", "TEMP"],
        ]

    ssm1 = flopy.mf6.ModflowGwtssm(
        gwt,
        sources=sourcerecarray_gwt,
        pname="ssm-" + str(ct + 1),
    )

    oc2 = flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=gwtname + ".ucn",
        saverecord=[("CONCENTRATION", "ALL")],
        pname="oc-" + str(ct + 1),
    )

    # ---------------------------------------------------------
    # GWT - energy transport model using parameter equivalents
    # ---------------------------------------------------------

    adv2 = flopy.mf6.ModflowGwtadv(
        gwt2,
        scheme="UPSTREAM",
        pname="adv-" + str(ct + 2),
    )
    mst2 = flopy.mf6.ModflowGwtmst(
        gwt2,
        porosity=porosity,
        pname="mst-" + str(ct + 2),
    )
    dsp2 = flopy.mf6.ModflowGwtdsp(
        gwt2,
        alh=2.5,
        atv=2.5,
        ath1=2.5,
        pname="dsp-" + str(ct + 2),
        diffc=0.0,
    )
    ssm2 = flopy.mf6.ModflowGwtssm(
        gwt2,
        sources=sourcerecarray_gwte,
        pname="ssm-" + str(ct + 2),
    )

    # Output control for the "temperature" that uses the old
    # parameter equivalents approach
    oc3 = flopy.mf6.ModflowGwtoc(
        gwt2,
        concentration_filerecord=gwtname2 + ".ucn",
        saverecord=[("CONCENTRATION", "ALL")],
        pname="oc-" + str(ct + 2),
    )

    return sim


def build_models(idx, test):  # test):
    # Base simulation and model name and workspace
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(
        sim_name="vsc_maw_test",
        exe_name="mf6",
        version="mf6",
        sim_ws=ws,
    )

    # Set up time discretization
    perioddata = [
        (1.0, 1, 1.0),
        (90.0, 1, 1.0),
    ]
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units=time_units,
        nper=len(perioddata),
        perioddata=perioddata,
    )

    gwfname1 = "gwf-wel-" + cases[idx]
    gwtname1 = "gwt1-" + cases[idx]
    gwtname2 = "gwt2-" + cases[idx]
    sim = add_models(sim, gwfname1, gwtname1, gwtname2, ct=1, wel_on=True)

    gwfname2 = "gwf-maw-" + cases[idx]
    gwtname3 = "gwt3-" + cases[idx]
    gwtname4 = "gwt4-" + cases[idx]
    sim = add_models(
        sim,
        gwfname2,
        gwtname3,
        gwtname4,
        ct=4,
        maw_on=True,
    )

    # Exchange files
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname1,
        exgmnameb=gwtname1,
        filename=f"{gwfname1}.gwfgwt",
        pname="gwfgwt1",
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname1,
        exgmnameb=gwtname2,
        filename=f"{gwfname1}.gwfgwte",
        pname="gwfgwt2",
    )

    # Exchange files for the second GWF model
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname2,
        exgmnameb=gwtname3,
        filename=f"{gwfname2}.gwfgwt",
        pname="gwfgwt1",
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname2,
        exgmnameb=gwtname4,
        filename=f"{gwfname2}.gwfgwte",
        pname="gwfgwt2",
    )
    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read flow results from model
    name = cases[idx]
    gwfname1 = "gwf-wel-" + cases[idx]
    gwtname1 = "gwt1-" + cases[idx]
    gwtname2 = "gwt2-" + cases[idx]

    gwfname2 = "gwf-maw-" + cases[idx]
    gwtname3 = "gwt3-" + cases[idx]
    gwtname4 = "gwt4-" + cases[idx]

    fpth1 = os.path.join(test.workspace, f"{gwfname1}.hds")
    fpth2 = os.path.join(test.workspace, f"{gwfname2}.hds")

    try:
        # load heads
        gwf1obj = flopy.utils.HeadFile(fpth1, precision="double")
        hds1 = gwf1obj.get_alldata()
    except:
        assert False, f'could not load head data from "{fpth1}"'

    try:
        # load heads
        gwf2obj = flopy.utils.HeadFile(fpth2, precision="double")
        hds2 = gwf2obj.get_alldata()
    except:
        assert False, f'could not load head data from "{fpth2}"'

    msg0 = (
        "Comparing heads in 2 GWF models using WEL and MAW, respectively, plus "
        "viscosity adjustments based on solute and temperature effects, should "
        "give the same result based on model setup, but doesn't. Layer "
    )

    lay = 2
    assert np.allclose(hds1[-1, lay - 1], hds2[-1, lay - 1], atol=0.03), msg0 + str(lay)
    lay = 3
    assert np.allclose(hds1[-1, lay - 1], hds2[-1, lay - 1], atol=0.012), msg0 + str(
        lay
    )
    lay = 7
    assert np.allclose(hds1[-1, lay - 1], hds2[-1, lay - 1], atol=7e-3), msg0 + str(lay)


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
