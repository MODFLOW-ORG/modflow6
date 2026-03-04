"""
Test problem for GWT with MVT active.  Another autotest will check this autotest
in parallel.

Test that the MVT successfully transfers solute from one reach to the next.  It
is expected that is will work fine in serial, the check is will it work in
parallel.

 Model configuration

                      Model 1                          Model 2
                     ---------                        ---------

               +---------+---------+            +---------+---------+
              /         /         /|           /         /         /|
SFR/SFT -> ==/=========/=========/=> MVR/MVT ======================>|
            /         /         /  +         /         /         /  +
           +---------+---------+  /         +---------+---------+  /
           |         |         | /          |         |         | /
           |         |         |/           |         |         |/
           +---------+---------+            +---------+---------+


"""

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

cases = ["mvt"]

# model units
length_units = "meters"
time_units = "days"

# static parameters
nrow = 1
ncol = 2
nlay = 1
top = 1.0
botm = 0.0
delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)
ss = 1e-6  # Specific storage
sy = 0.20  # Specific Yield
prsity = 0.20  # Porosity
nper = 1  # Number of periods
perlen = [1000.0]  # Simulation time ($days$)
nstp = [1]  # 10 day transient time steps
steady = {0: False}
transient = {0: True}
strthd = 0.8
strt_conc = 100.0
scheme = "UPSTREAM"

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible
icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
strt_temp1 = 4.0
strt_temp2 = 34.0
dispersivity = 0.0  # dispersion (remember, 1D model)

# SFR related parameters
surf_Q_in = 1.0
nreaches = ncol
rlen = delr
rwid = 9.0
roughness = 0.035
rbth = 0.1
rhk = 0.0
ustrf = 1.0
ndv = 0
strm_incision = 0.1

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-10, 1.0
ttsmult = 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen[i], nstp[i], ttsmult))


# Functions for creating MODFLOW 6 GWF/GWT models with an SFR/SFT &
# MVR/MVT connection


def build_gwf_model(sim, gwfname, mod_num):
    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    # Instantiating MODFLOW 6 discretization package
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
        idomain=idomain,
        pname="DIS-" + str(mod_num),
        filename=f"{gwfname}.dis",
    )

    # Instantiating MODFLOW 6 storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=ss,
        sy=sy,
        iconvert=iconvert,
        steady_state=steady,
        transient=transient,
        pname="STO-" + str(mod_num),
        filename=f"{gwfname}.sto",
    )

    # Instantiating node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        pname="NPF-" + str(mod_num),
        filename=f"{gwfname}.npf",
    )

    # Instantiating initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=strthd,
        pname="IC-" + str(mod_num),
        filename=f"{gwfname}.ic",
    )

    # Instantiate streamflow routing
    packagedata = []

    for i in np.arange(ncol):
        ncon = 2
        if i in [0, ncol - 1]:
            ncon = 1
        packagedata.append(
            [
                i,  # irch
                (0, 0, i),  # lay, row, col
                delr,  # reach length
                rwid,  # reach width
                0.001,  # slope
                top - (i * strm_incision),  # streambed elevation
                rbth,  # streambed thickness
                rhk,  # streambed K
                roughness,  # Manning's "n"
                ncon,  # number of connections
                ustrf,  # no idea
                ndv,  # number of diversions
            ]
        )

    con_data = []
    for irno in range(ncol):
        if irno == 0:
            t = (irno, -(irno + 1))
        elif irno == ncol - 1:
            t = (irno, irno - 1)
        else:
            t = (irno, irno - 1, -(irno + 1))
        con_data.append(t)

    sfrbndx = []
    for i in np.arange(len(packagedata)):
        if i == 0:
            sfrbndx.append([i, "INFLOW", surf_Q_in])

    sfr_perioddata = {0: sfrbndx}

    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        mover=True,
        nreaches=len(packagedata),
        packagedata=packagedata,
        connectiondata=con_data,
        perioddata=sfr_perioddata,
        pname="SFR-" + str(mod_num),
        filename=f"{gwfname}.sfr",
    )

    # Instantiating output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC-" + str(mod_num),
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, gwf


def build_gwt_model(sim, gwtname, mod_num):
    # Instantiating GWT model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, model_nam_file=f"{gwtname}.nam")

    gwt.name_file.save_flows = True

    # Instantiating discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
        pname="DIS-" + str(mod_num),
        filename=f"{gwtname}.dis",
    )

    # Instantiating initial concentrations
    flopy.mf6.ModflowGwtic(
        gwt,
        strt=strt_conc - ((mod_num - 1) * strt_conc),
        pname="IC-" + str(mod_num),
        filename=f"{gwtname}.ic",
    )

    # Instantiating advection package
    flopy.mf6.ModflowGwtadv(
        gwt, scheme=scheme, pname="ADV-" + str(mod_num), filename=f"{gwtname}.adv"
    )

    # Instantiating dispersion package
    flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        alh=dispersivity,
        ath1=dispersivity,
        pname="CND-" + str(mod_num),
        filename=f"{gwtname}.cnd",
    )

    # Instantiating mass storage package
    flopy.mf6.ModflowGwtmst(
        gwt,
        save_flows=True,
        porosity=prsity,
        pname="MST-" + str(mod_num),
        filename=f"{gwtname}.mst",
    )

    # Instantiate source-sink mixing package
    flopy.mf6.ModflowGwtssm(
        gwt, sources=[[]], pname="SSM-" + str(mod_num), filename=f"{gwtname}.ssm"
    )

    # Instantiate streamflow transport package
    sftpackagedata = []
    for irno in range(ncol):
        t = (irno, strt_conc - ((mod_num - 1) * strt_conc * 0.5), f"myreach{irno + 1}")
        sftpackagedata.append(t)

    sftperioddata = [
        (0, "STATUS", "ACTIVE"),
        (1, "STATUS", "ACTIVE"),
        (0, "INFLOW", strt_conc - ((mod_num - 1) * strt_conc)),
    ]

    sft_obs = {
        (gwtname + ".sft.obs.csv",): [
            (f"sft-{i + 1}-conc", "CONCENTRATION", i + 1) for i in range(ncol)
        ]
        + [(f"sft-{i + 1}-extinfq", "EXT-INFLOW", i + 1) for i in range(ncol)],
    }
    sft_obs["digits"] = 10
    sft_obs["print_input"] = True
    sft_obs["filename"] = gwtname + ".sft.obs"

    sft = flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        packagedata=sftpackagedata,
        reachperioddata=sftperioddata,
        observations=sft_obs,
        flow_package_name="SFR-" + str(mod_num),
        pname="SFT-" + str(mod_num),
        filename=f"{gwtname}.sft",
    )

    # Instantiate heat transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        pname="OC",
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, gwt


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating time discretization
    flopy.mf6.ModflowTdis(
        sim,
        nper=nper,
        perioddata=tdis_rc,
        time_units=time_units,
        filename=f"{name}.tdis",
    )

    gwf_models = []
    gwt_models = []

    for i in [1, 2]:
        # generate names for each model
        gwfname = "gwf" + str(i) + "-" + name
        gwtname = "gwt" + str(i) + "-" + name

        # gwf
        sim, gwf = build_gwf_model(sim, gwfname, i)
        # gwt
        sim, gwt = build_gwt_model(sim, gwtname, i)

        gwf_models.append(gwf.name)
        gwt_models.append(gwt.name)

        # Instantiating flow-transport exchange mechanism
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            pname="gwfgwt-" + str(i),
            filename=f"{gwtname}.gwfgwt" + str(i),
        )

    # Instantiating solver for both flow models
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf_models[0]])
    sim.register_ims_package(imsgwf, [gwf_models[1]])

    imsgwt = flopy.mf6.ModflowIms(
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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt_models[0]])
    sim.register_ims_package(imsgwt, [gwt_models[1]])

    # now the simulation-level GWF-GWF and GWT-GWT movers
    # GWF-GWF

    # add a gwf-gwf exchange
    gwfgwf_data = [
        ((0, 0, ncol - 1), (0, 0, 0), 1, delr / 2.0, delr / 2.0, delc, 0.0, delr)
    ]

    mvr_filerecord = f"{gwf_models[0]}-{gwf_models[1]}.exg.mvr"
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwf_models[0],
        exgmnameb=gwf_models[1],
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvr_filerecord=mvr_filerecord,
        filename=f"{gwf_models[0]}-{gwf_models[1]}.exg",
    )

    # simulation GWF-GWF Mover
    maxmvr, maxpackages = 1, 2
    mvrpack_sim = [[gwf_models[0], "SFR-1"], [gwf_models[1], "SFR-2"]]
    mvrspd = [
        [gwf_models[0], "SFR-1", ncol - 1, gwf_models[1], "SFR-2", 0, "FACTOR", 1.00]
    ]
    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackages,
        packages=mvrpack_sim,
        perioddata=mvrspd,
        filename=mvr_filerecord,
    )

    # GWT-GWT
    mvt_filerecord = f"{gwt_models[0]}-{gwt_models[1]}.exg.mvr"
    gwtgwt = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwt_models[0],
        exgmnameb=gwt_models[1],
        gwfmodelname1=gwf_models[0],
        gwfmodelname2=gwf_models[1],
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvt_filerecord=mvt_filerecord,
        filename=f"{gwt_models[0]}-{gwt_models[1]}.exg",
    )

    # simulation GWT-GWT Mover
    gwtgwt.mvt.initialize(filename=mvt_filerecord)

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwtname = "gwt" + str(2) + "-" + name

    ws = test.workspace
    fname = os.path.join(ws, gwtname + ".lst")
    assert os.path.isfile(fname), "File not found: " + fname

    srch_str1 = "SFT PACKAGE (SFT-2) CONCENTRATION FOR EACH CONTROL VOLUME"
    srch_str2 = " MYREACH1 "

    with open(fname, "r") as f:
        for line in f:
            if srch_str1 in line:
                while srch_str2 not in line:
                    line = next(f)
                break

    # process the last line read
    m_arr = line.strip().split()

    # if water is successfully moved from model 1 to model 2, the
    # stream concentration should be cut in half
    assert float(m_arr[-1]) == strt_conc / 2, (
        "SFT concentration for model2, reach 1 not equal to 50.0"
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
