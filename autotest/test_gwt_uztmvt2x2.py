"""
2x2 test problem for GWE

Test the rejected infiltration mover transport values

Model configuration:  2 SFR reaches exist in the first row.
                      2 UZF objects exist in the second row
                      2 MVR/MVT connections transfer rejected infiltration
                        from UZF to SFR

                   Col 1   Col 2
                 +-------+-------+
        Row 1  /       /       /|
  SFR         /       /       / |
Channel-> =========================
            /       /  ^    /   |
           +-------+---+---+    |
          /       /    |  /|    +
 Row 2   /       /  rej. / |   /
        /       /  inf. /  |  /
       /       /       /   | /
      +-------+-------+    |/
      |       |       |    +
      |       |       |   /
      |       |       |  /
      |       |       | /
      |       |       |/
      +-------+-------+

    Profile View:
    ------------
              +---- rej inf.
              v  +---------+
       +---------+         |
       |         |         |
       |         +---------+
       +---------+  Row 2
          Row 1
"""

# Imports

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

scheme = "UPSTREAM"
# scheme = "TVD"

cases = ["uztmvt"]  # 2-cell model, horizontally connected with staggered alignment

nrow = 2
ncol = 2
nlay = 1
top = np.array([[[1.0, 1.0], [1.5, 1.5]]], dtype=float)
bot = np.array([[[0.0, 0.0], [0.5, 0.5]]], dtype=float)
strthd = np.array([[[0.01, 0.01], [0.51, 0.51]]], dtype=float)

# Model units
length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

delr = 1.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
k11 = 1.0e3  # Horizontal hydraulic conductivity ($m/d$)
vks = 0.01  # vertical hydraulic conductivity of the unsaturated zone
ss = 1e-6  # Specific storage
sy = 0.30  # Specific Yield
prsity = 0.20  # Porosity
nper = 2  # Number of periods
perlen = [1, 1]  # Simulation time ($days$)
nstp = [1, 1]  # 10 day transient time steps
ttsmult = 1.0
steady = {0: False}
transient = {0: True}

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
idomain = 1  # All cells included in the simulation
iconvert = 1  # All cells are convertible

icelltype = 1  # Cell conversion type (>1: unconfined)

chdlist = []
chdlist.append([(0, 0, 0), 0.51])
chdlist.append([(0, 0, 1), 0.51])

# Set some static transport related model parameter values
strt_conc = 10.0
strt_uz_conc = 1.0
dispersivity = 0.0

# UZF related parameters
thtr = 0.05  # Residual water content
thts = sy  # Saturated water content
thti = 0.05  # Initial water content (unsaturated zone)
eps = 7.1  # Brooks-Corey epsilon

# GWT related parameters
al = 0.0  # Longitudinal dispersivity ($m$)
rhob = 1.5  # Bulk density of layer 1 ($g/cm^3$)
Kd = 0.176  # Distribution coefficient ($cm^3/g$)

# UZF/UZT related input
surfdep = 1.0e-5
finf1 = 1.01
finf2 = 2.01
pet = 0.0
extdp = 0.5

uzf_pkdat = [
    [0, (0, 1, 0), 1, -1, surfdep, vks, thtr, thts, thti, eps],
    [1, (0, 1, 1), 1, -1, surfdep, vks, thtr, thts, thti, eps],
]
uzf_spd = {
    0: [
        [0, finf1, pet, extdp, thtr, 0.0, 0.0, 0.0],
        [1, finf1, pet, extdp, thtr, 0.0, 0.0, 0.0],
    ],
    1: [
        [0, finf2, pet, extdp, thtr, 0.0, 0.0, 0.0],
        [1, finf2, pet, extdp, thtr, 0.0, 0.0, 0.0],
    ],
}
concCell = []
concCell.append(11.1)
concCell.append(22.2)
uzt_pkdat = [(0, strt_uz_conc), (1, strt_uz_conc)]  # ifno, strt_conc
uzt_perdat = [
    (0, "INFILTRATION", concCell[0]),
    (1, "INFILTRATION", concCell[1]),
]

# SFR/SFT related input
conn_dat = [[0, -1], [1, 0]]
sfr_pkdat = []
nreaches = 2
rlen = 1.0
rwid = 1.0
roughness = 0.03
rbth = 0.1
rhk = 0.0
slope = 0.001
ustrf = 1.0
ndv = 0
# reach 1
rp = [0, (0, 0, 0), rlen, rwid, slope, top[0, 0, 0], rbth, rhk, roughness, 1, ustrf, 0]
sfr_pkdat.append(rp)
# reach 2
rp = [1, (0, 0, 1), rlen, rwid, slope, top[0, 0, 1], rbth, rhk, roughness, 1, ustrf, 0]
sfr_pkdat.append(rp)

sfr_perdat = {0: [0, "INFLOW", 1.0]}
sft_pkdat = [[0, 0.0], [1, 0.0]]
sft_perdat = {
    0: [[0, "STATUS", "ACTIVE"], [0, "INFLOW", 0.0]],
    1: [[0, "STATUS", "ACTIVE"], [0, "INFLOW", 1.0]],
}

# MVR input
mvr_pkdat = []
mvr_pkdat.append(["UZF-1", 0, "SFR-1", 0, "FACTOR", 1.0])
mvr_pkdat.append(["UZF-1", 1, "SFR-1", 1, "FACTOR", 1.0])
mvr_packages = [
    ("UZF-1",),
    ("SFR-1",),
    # ("DRN-1",),
]

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-10, 1.0

# Set up temporal data used by TDIS file
tdis_rc = []
for i in np.arange(nper):
    tdis_rc.append((perlen[i], nstp[i], ttsmult))

# ### Create MODFLOW 6 GWT
#
# No GWF, only Heat conduction simulated


def build_models(idx, test):
    # Base MF6 GWF model type
    ws = test.workspace
    name = cases[idx]

    print(f"Building MF6 model...{name}")

    # generate names for each model
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    # Instantiating MODFLOW 6 solver for flow model
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
    sim.register_ims_package(imsgwf, [gwfname])

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
        botm=bot,
        idomain=1,
        pname="DIS-1",
        filename=f"{gwfname}.dis",
    )

    # Instantiating MODFLOW 6 storage package
    # flopy.mf6.ModflowGwfsto(
    #    gwf,
    #    ss=ss,
    #    sy=sy,
    #    iconvert=iconvert,
    #    steady_state=steady,
    #    transient=transient,
    #    pname="STO",
    #    filename="{}.sto".format(gwfname),
    # )

    # CHD
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist,
        pname="CHD",
        filename=f"{gwfname}.chd",
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        pname="NPF-1",
        filename=f"{gwfname}.npf",
    )

    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=ss,
        sy=sy,
        transient={0: True},
        pname="STO",
        filename=f"{gwfname}.sto",
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=strthd,
        filename=f"{gwfname}.ic",
    )

    # UZF
    flopy.mf6.ModflowGwfuzf(
        gwf,
        mover=True,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=False,
        simulate_et=False,
        unsat_etwc=False,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{gwfname}.uzf.bud",
        pname="UZF-1",
        filename=f"{gwfname}.uzf",
    )

    # SFR
    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=86400.0,
        budget_filerecord=f"{gwfname}.sfr.bud",
        mover=True,
        nreaches=nreaches,
        packagedata=sfr_pkdat,
        connectiondata=conn_dat,
        perioddata=sfr_perdat,
        pname="SFR-1",
        filename=f"{gwfname}.sfr",
    )

    # MVR
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=len(mvr_pkdat),
        budget_filerecord=f"{gwfname}.mvr.bud",
        maxpackages=len(mvr_packages),
        print_flows=True,
        packages=mvr_packages,
        perioddata=mvr_pkdat,
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        pname="OC-1",
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # ----------------------------------------------------
    # Instantiating MODFLOW 6 GWT model
    # ----------------------------------------------------
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, model_nam_file=f"{gwtname}.nam")
    gwt.name_file.save_flows = True
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
    sim.register_ims_package(imsgwt, [gwt.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=bot,
        idomain=1,
        pname="DIS-1",
        filename=f"{gwtname}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGwtic(gwt, strt=strt_conc, pname="IC-1", filename=f"{gwtname}.ic")

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGwtadv(
        gwt, scheme=scheme, pname="ADV-2", filename=f"{gwtname}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        alh=dispersivity,
        alv=dispersivity,
        ath1=dispersivity,
        atv=dispersivity,
        pname="CND-1",
        filename=f"{gwtname}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGwtmst(
        gwt,
        save_flows=True,
        porosity=prsity,
        sorption="linear",
        bulk_density=rhob,
        distcoef=Kd,
        pname="MST-1",
        filename=f"{gwtname}.mst",
    )

    flopy.mf6.modflow.ModflowGwtuzt(
        gwt,
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".uzt.bin",
        budget_filerecord=gwtname + ".uzt.bud",
        packagedata=uzt_pkdat,
        uztperioddata=uzt_perdat,
        flow_package_name="UZF-1",
        pname="UZT-1",
        filename=f"{gwtname}.uzt",
    )

    flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=False,
        flow_package_name="SFR-1",
        print_concentration=True,
        save_flows=True,
        print_flows=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sft_pkdat,
        reachperioddata=sft_perdat,
        pname="SFT-1",
        filename=f"{gwtname}.sft",
    )

    flopy.mf6.modflow.ModflowGwtmvt(
        gwt,
        save_flows=True,
        budget_filerecord=gwtname + ".mvt.bud",
        filename=f"{gwtname}.mvt",
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        pname="OC-2",
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    sourcerecarray = [[]]
    flopy.mf6.ModflowGwessm(gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm")

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        pname="GWFGWT",
        filename=f"{gwtname}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    print("evaluating results...")

    # read transport results from GWE model
    name = cases[idx]
    gwfname = "gwf-" + name
    gwtname = "gwt-" + name

    # overland flow
    fpth = os.path.join(test.workspace, f"{gwfname}.mvr.bud")

    try:
        # load mover flows
        fobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        txtnames = fobj.get_unique_record_names()
        times = fobj.get_times()
        mvrdat = {}
        for tm in times:
            for txtname in txtnames:
                # "MOVER-FLOW" should be the only text name
                flowX = fobj.get_data(totim=tm, text=txtname)
                mvrdat.update({tm: flowX})
    except:
        assert False, f'could not load flow data from "{fpth}"'

    fpth = os.path.join(test.workspace, f"{gwtname}.mvt.bud")
    try:
        # load transport mover results
        mobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        txtnames = mobj.get_unique_record_names()
        times = fobj.get_times()
        mvtdat = {}
        for tm in times:
            for txtname in txtnames:
                massX = mobj.get_data(totim=tm, text=txtname)
                mvtdat.update({tm: massX})
    except:
        assert False, f'could not load transport data from "{fpth}"'

    # Mover amounts are representative of rejected infiltration only
    # (there are no other MVR connections established). A known amount of
    # water will be transferred based on the infiltration and VKS: 1 unit
    # of water in stress period 1 and 2 units in the 2nd stress period
    msg0 = "Rejected infiltration transfer amount not as expected"
    for x in np.arange(len(mvrdat)):
        for y in np.arange(len(mvrdat[x + 1])):
            if len(mvrdat[x + 1][y]) == 0:
                continue
            else:
                for z in np.arange(len(mvrdat[x + 1][y])):
                    assert np.isclose(abs(mvrdat[x + 1][y][z][-1]), x + 1.0), msg0

    # Transport mover (MVT) amounts are known quantities
    msg1 = "Rejected infiltration transfer mass amount not as expected"
    for x in np.arange(len(mvtdat)):
        for y in np.arange(len(mvtdat[x + 1])):
            if len(mvtdat[x + 1][y]) == 0:
                continue
            else:
                for z in np.arange(len(mvtdat[x + 1][y])):
                    assert np.isclose(
                        mvtdat[x + 1][y][z][-1], (x + 1.0) * concCell[z]
                    ), msg1

    # Stream mass flows should be known
    fpth = os.path.join(test.workspace, f"{gwtname}.sft.bud")
    try:
        # load mover flows
        sfrobj = flopy.utils.CellBudgetFile(fpth, precision="double")
        txtnames = sfrobj.get_unique_record_names()
        txtnames = [nm.decode("utf-8").strip() for nm in txtnames]
        times = sfrobj.get_times()
        JAdat = {}
        fromMvrDat = {}
        for tm in times:
            for txtname in txtnames:
                if txtname == "FLOW-JA-FACE":
                    dat = sfrobj.get_data(totim=tm, text=txtname)
                    JAdat.update({tm: dat})
                if txtname == "FROM-MVR":
                    dat = sfrobj.get_data(totim=tm, text=txtname)
                    fromMvrDat.update({tm: dat})

    except:
        assert False, f'could not load flow data from "{fpth}"'

    msg2 = "Mass received by SFR ('FROM-MVR') not as expected"
    for x in np.arange(len(fromMvrDat)):
        for y in np.arange(len(fromMvrDat[x + 1][0])):
            if fromMvrDat[x + 1][0][y][-1] == concCell[y]:
                continue
            else:
                for z in np.arange(len(mvtdat[x + 1][y])):
                    assert np.isclose(
                        mvtdat[x + 1][y][z][-1], (x + 1.0) * concCell[z]
                    ), msg2


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
