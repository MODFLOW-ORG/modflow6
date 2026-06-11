"""
Test problem for GWT with MVT active.  Another autotest will check this autotest
in parallel.

Test that the MVT successfully transfers solute from two separate streams. First
check is that the solute from stream 1 goes from model 1 to model 2 and second,
check that the solute in stream 2 goes from model 2 to model 1.  The idea is to
test that when run in parallel mode MVT will keep the direction of the solute
move correct.


 Model configuration

                      Model 1                          Model 2
                     ---------                        ---------

               +---------+---------+            +---------+---------+
               |         |         |    --->    |         |         |
   SFR/SFT -> =|=========|=========|=> MVR/MVT =|=========|=========|==> Outlet 1
               |         |         |            |         |         |
               +---------+---------+            +---------+---------+
               |         |         |   <---     |         |         |
 Outlet 2 <- <=|=========|=========|= MVR/MVT <=|=========|=========|== SFR/SFT
               |         |         |            |         |         |
               +---------+---------+            +---------+---------+
               |         |         |    --->    |         |         |
   SFR/SFT -> =|=========|=========|=> MVR/MVT =|=========|=========|==> Outlet 3
               |         |         |            |         |         |
               +---------+---------+            +---------+---------+

  The test is: does MVT keep track of which is the receiver
                     and which is the provider when there is multi-directional
                     Movers and Transport Movers

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
nrow = 3
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
perlen = [10000.0]  # Simulation time ($days$)
nstp = [1]  # 10 day transient time steps
steady = {0: False}
transient = {0: True}
strthd = 0.8
strt_conc = [100.0, 200.0, 300.0]
scheme = "UPSTREAM"

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible
icelltype = 1  # Cell conversion type (>1: unconfined)

# Set some static transport related model parameter values
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


def add_strm_rch(num_segs, num_rchs, left_model):
    packagedata = []
    con_data = []
    sfrbndx = []
    ifno = 0

    # generate a list for number of rchs each "segment" will have
    rchs = []
    for dum in range(num_rchs):
        rchs.append(dum)

    # iterate over each individual stream
    for s in range(num_segs):  # number of distinct streams in the model
        # print("segment: " + str(s))
        srchs = rchs
        if s % 2 == 1:  # when evaluates to odd, stream flow is right to left
            srchs = rchs[::-1]

        # iterate over each reach within each stream
        for ct, i in enumerate(srchs):
            # set variables to go into packagedata
            cellid = (0, s, i)
            ncon = 2
            if i in [0, num_rchs - 1]:
                ncon = 1

            packagedata.append(
                [
                    ifno,  # irch
                    cellid,  # lay, row, col
                    delr,  # reach length
                    rwid,  # reach width
                    0.001,  # slope
                    top - strm_incision,  # streambed elevation
                    rbth,  # streambed thickness
                    rhk,  # streambed K
                    roughness,  # Manning's "n"
                    ncon,  # number of connections
                    ustrf,  # no idea
                    ndv,  # number of diversions
                ]
            )

            if ct == 0:
                t = (ifno, -(ifno + 1))
            elif ct == num_rchs - 1:
                t = (ifno, ifno - 1)
            else:
                t = (ifno, ifno - 1, -(ifno + 1))
            con_data.append(t)

            # set inflow > 0 if first reach of stream. For flow direction,
            # see header information
            print("which model: " + str(left_model))
            if left_model:
                if s % 2 == 0 and ct == 0:
                    sfrbndx.append([ifno, "inflow", surf_Q_in])
            else:
                if s % 2 == 1 and ct == 0:
                    sfrbndx.append([ifno, "inflow", surf_Q_in])

            # increment
            ifno += 1

    return packagedata, con_data, sfrbndx


def get_only_digit(s):
    # Filter for all digit characters
    digits_list = [char for char in s if char.isdigit()]  #

    # Check if exactly one digit was found
    if len(digits_list) == 1:
        # Return the digit as an integer (or string, if preferred)
        return int(digits_list[0])
    elif len(digits_list) == 0:
        return "No digit found"
    else:
        return "Multiple digits found"


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
    flow_dir = True
    if mod_num != 1:
        flow_dir = False

    packagedata, con_data, sfrbndx = add_strm_rch(nrow, ncol, flow_dir)

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
        strt=strt_conc[mod_num - 1],
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
    for irno in range(nrow * ncol):
        t = (irno, strt_conc[mod_num - 1], f"myreach{irno + 1}")
        sftpackagedata.append(t)

    sftperioddata = []

    for r in range(nrow):
        for c in range(ncol):
            if mod_num == 1:
                if r % 2 == 0 and c == 0:
                    sftperioddata.append(
                        (r * ncol, "INFLOW", strt_conc[r]),
                    )
            else:
                if r % 2 == 1 and c == 1:
                    sftperioddata.append(
                        (r * ncol, "INFLOW", strt_conc[r]),
                    )

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


def add_flow_model_solvers(sim, gwf_models, name):
    # Instantiate solver for the flow model
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
        filename=f"gwf-{name}.ims",
    )
    for idx in range(len(gwf_models)):
        sim.register_ims_package(imsgwf, [gwf_models[idx]])

    return sim


def add_transport_model_solvers(sim, gwt_models, name):
    # Instantiate solver for the transport models
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
        filename=f"gwt-{name}.ims",
    )
    for idx in range(len(gwt_models)):
        sim.register_ims_package(imsgwt, [gwt_models[idx]])

    return sim


def add_a_gwf_gwf_conn(sim, modnam_up, modnam_dn):
    gwfgwf_data = []
    for rw in range(nrow):
        gwfgwf_data.append(
            ((0, rw, ncol - 1), (0, rw, 0), 1, delr / 2.0, delr / 2.0, delc, 0.0, delr)
        )

    mvr_filerecord = f"{modnam_up}-{modnam_dn}.exg.mvr"
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=modnam_up,
        exgmnameb=modnam_dn,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvr_filerecord=mvr_filerecord,
        filename=f"{modnam_up}-{modnam_dn}.exg",
    )

    # simulation GWF-GWF Mover
    up_mod_id = get_only_digit(modnam_up)
    dn_mod_id = get_only_digit(modnam_dn)

    maxmvr, maxpackages = 3, 2
    mvrpack_sim = [
        [modnam_up, "SFR-" + str(up_mod_id)],
        [modnam_dn, "SFR-" + str(dn_mod_id)],
    ]

    mvrspd = []
    for mod in [1, 2]:
        for rw in range(nrow):
            if mod == 1 and rw in [0, 2]:
                mvrspd.append(
                    [
                        modnam_up,
                        "SFR-" + str(up_mod_id),
                        rw * ncol + (ncol - 1),
                        modnam_dn,
                        "SFR-" + str(dn_mod_id),
                        rw * ncol,
                        "FACTOR",
                        1.00,
                    ]
                )
            elif mod == 2 and rw == 1:
                mvrspd.append(
                    [
                        modnam_dn,
                        "SFR-" + str(dn_mod_id),
                        rw * ncol + (ncol - 1),
                        modnam_up,
                        "SFR-" + str(up_mod_id),
                        rw * ncol,
                        "FACTOR",
                        1.00,
                    ]
                )

    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackages,
        packages=mvrpack_sim,
        perioddata=mvrspd,
        filename=mvr_filerecord,
    )

    return sim, gwfgwf_data


def add_a_gwt_gwt_conn(sim, modnam_up, modnam_dn, gwfgwf_data):
    mvt_filerecord = f"{modnam_up}-{modnam_dn}.exg.mvr"
    gwf_modnam_up = modnam_up.replace("gwt", "gwf")
    gwf_modnam_dn = modnam_dn.replace("gwt", "gwf")
    gwtgwt = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        nexg=len(gwfgwf_data),
        exgmnamea=modnam_up,
        exgmnameb=modnam_dn,
        gwfmodelname1=gwf_modnam_up,
        gwfmodelname2=gwf_modnam_dn,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        mvt_filerecord=mvt_filerecord,
        filename=f"{modnam_up}-{modnam_dn}.exg",
    )

    # simulation GWT-GWT Mover
    gwtgwt.mvt.initialize(filename=mvt_filerecord)

    return sim


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

    # Instantiate solver for the flow models
    sim = add_flow_model_solvers(sim, gwf_models, name)

    # Instantiate solver for the transport models
    sim = add_transport_model_solvers(sim, gwt_models, name)

    # Next, add the simulation-level GWF-GWF and GWT-GWT movers:
    # GWF-GWF
    sim, gwfgwf_data = add_a_gwf_gwf_conn(sim, gwf_models[0], gwf_models[1])

    # GWT-GWT
    sim = add_a_gwt_gwt_conn(sim, gwt_models[0], gwt_models[1], gwfgwf_data)

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]

    # check both transport models
    for modnum in [1, 2]:
        gwtname = "gwt" + str(modnum) + "-" + name

        ws = test.workspace
        fname = os.path.join(ws, gwtname + ".lst")
        assert os.path.isfile(fname), "File not found: " + fname

        srch_str1 = (
            " SFT PACKAGE (SFT-"
            + str(modnum)
            + ") CONCENTRATION FOR EACH CONTROL VOLUME "
        )

        ifno = 1
        with open(fname, "r") as f:
            # advance to the line with the first srch string
            for line in f:
                if srch_str1 in line:
                    break

            # now cycle through the reach output and confirm
            for r in range(nrow):
                for c in range(ncol):
                    srch_str2 = " MYREACH" + str(ifno) + " "

                    while srch_str2 not in line:
                        line = next(f)

                    # process the last line read
                    m_arr = line.strip().split()

                    print()
                    # if water is successfully moved from models 1 & 2 to model 3,
                    #  the stream concentration should be preserved across models
                    idx = r * ncol + c
                    if idx in [0, 1]:
                        check_val = strt_conc[0]
                    elif idx in [2, 3]:
                        check_val = strt_conc[1]
                    elif idx in [4, 5]:
                        check_val = strt_conc[2]

                    assert float(m_arr[-1]) == check_val, (
                        "SFT concentration for model"
                        + str(modnum)
                        + ", reach "
                        + str(idx)
                        + ", not equal to "
                        + str(check_val)
                    )
                    ifno += 1


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
