"""
Uses a sloped surface to test MVR/MVE connections. Connection types include
 1) DRN to SFR
 2) UZF to SFR
 3) DRN to UZF
 4) UZF to UZF

Stress period 1 shouldn't have any DRN to MVR flows
Stress period 1 should include UZF to UZF and UZF to SFR flows
Stress period 2 shouldn't have any type of TO-MVR flows
Stress period 3 should include DRN to UZF and DRN to SFR flows
Stress period 4 should include all types (listed above) of TO-MVR flows

             (Not to scale: shows model grid design in profile view)

     +-----+
     |     +-----+                        precip
     |     |     +-----+        runoff      |             MVR connections: UZF -> SFR
     |     |     |     +-----+   ~~~>       v                              UZF -> UZF
     |     |     |     |     +-----+   ~~~>                                DRN -> UZF
     |     |     |     |     |     +-----+   ~~~>                          DRN -> SFR
     |     |     |     |     |     |     +-----+   ~~~>
     |     |     |     |     |     |     |     +-----+   ~~~>   SFR channel
     |     |     |     |     |     |     |     |     +-----+   /
     |     |     |     |     |     |     |     |     | gw  +-| |-+
     |     |     |     |     |     |     |     |     |dis Q| +-+ |
     +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
     |     |     |     |     |     |     |     |     |     |     |
     |     |     |     |     |     |     |     |     |     |     |
     +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
     |     |     |     |     |     |     |     |     |     |     |
     |     |     |     |     |     |     |     |     |     |     |
     +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+

 A note about "gw dis Q" in the 9th column: DRN connections using
 MVR/MVE have been set up for all active UZF cells; however, only
 the right-most active DRN connection actually transfers groundwater
 discharge to land surface to SFR.  That is, none of the DRN -> UZF
 connections accounting for cascading gw discharge to land surface
 should ever be anything but 0.0 (which is verified by the checks)

"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

include_NWT = False

cases = ["mve-01"]

iuz_cell_dict = {}
cell_iuz_dict = {}

nlay, nrow, ncol = 3, 3, 10
nper = 4
perlen = [1.0] * nper
nstp = [1] * nper
tsmult = [1.0] * nper

delr = 1.0
delc = 1.0
strt = 25

nouter, ninner = 300, 300
hclose, rclose, relax = 1e-6, 1e-3, 0.97

top = np.stack(
    [[30.9, 30.8, 30.6, 30.6, 30.5, 30.4, 30.3, 30.2, 30.1, 30.0] for _ in range(3)],
    axis=0,
)
botm = [25.0, 10, 0.0]

icelltype = 1
K = 1000.0
K33 = 100.0
ss = 1e-6
sy = 0.3
strt_temp = 1.0
scheme = "TVD"
dispersivity = 0.0
ktw = 0.5918
kts = 0.2700
rhow = 1000.0
cpw = 4183.0
lhv = 2500.0
cps = 760.0
rhos = 1500.0
# Thermal conductivity of the streambed material ($W/m/C$)
K_therm_strmbed = [1.5, 1.75, 2.0]
rbthcnd = 0.0001
prsity = sy
drn_depth = 2.0
ddrn = 1.0
drn_cond = 1e5

ghbelv1 = 28.0
ghbelv2 = 28.25
ghbcond = 50000.0
ghbtemp = 1.0

ghbspd = {
    0: [
        [0, 0, 0, ghbelv1, ghbcond, ghbtemp],
        [0, 1, 0, ghbelv1, ghbcond, ghbtemp],
        [0, 2, 0, ghbelv1, ghbcond, ghbtemp],
        [1, 0, 0, ghbelv1, ghbcond, ghbtemp],
        [1, 1, 0, ghbelv1, ghbcond, ghbtemp],
        [1, 2, 0, ghbelv1, ghbcond, ghbtemp],
        [2, 0, 0, ghbelv1, ghbcond, ghbtemp],
        [2, 1, 0, ghbelv1, ghbcond, ghbtemp],
        [2, 2, 0, ghbelv1, ghbcond, ghbtemp],
    ],
    2: [
        [0, 0, 0, ghbelv2, ghbcond, ghbtemp],
        [0, 1, 0, ghbelv2, ghbcond, ghbtemp],
        [0, 2, 0, ghbelv2, ghbcond, ghbtemp],
        [1, 0, 0, ghbelv2, ghbcond, ghbtemp],
        [1, 1, 0, ghbelv2, ghbcond, ghbtemp],
        [1, 2, 0, ghbelv2, ghbcond, ghbtemp],
        [2, 0, 0, ghbelv2, ghbcond, ghbtemp],
        [2, 1, 0, ghbelv2, ghbcond, ghbtemp],
        [2, 2, 0, ghbelv2, ghbcond, ghbtemp],
    ],
}

# ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
surfdep1 = 0.01
surfdep2 = 0.001
vks = 0.1
eps = 4.0
thti = 0.06
thtr = 0.05
thts = sy + thtr
iuzfbnd = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    ]
)

drn_pkdat = []
uzf_pkdat = []
uze_pkdat = []
uze_perdat = []
uzf_id_lkup = {}
ct = -1
for k in np.arange(nlay):
    landflg = 1
    if k > 0:
        landflg = 0

    for i in np.arange(nrow):
        for j in np.arange(ncol):
            if iuzfbnd[i, j] > 0:
                ct += 1
                if k < nlay - 1:
                    ivertconn = ct + iuzfbnd.sum()
                else:
                    ivertconn = -1

                # gwf
                uzf_pkdat.append(
                    [
                        ct,
                        (k, i, j),
                        landflg,
                        ivertconn,
                        surfdep1,
                        vks,
                        thtr,
                        thts,
                        thti,
                        eps,
                        "uzf" + str(ct + 1).zfill(2),
                    ]
                )
                # gwe
                uze_pkdat.append((ct, 0))
                uze_perdat.append([ct, "INFILTRATION", 1.0])
                # generate a lookup dictionary based on the top layer
                if k == 0:
                    drn_pkdat.append([(k, i, j), top[i, j] - drn_depth, drn_cond, ddrn])
                    uzf_id_lkup.update({(i, j): ct})


uze_perdat = {
    0: uze_perdat,
}

mvr_pkdat = []
for i in np.arange(nrow):
    for j in np.arange(ncol - 1):
        if iuzfbnd[i, j] > 0 and iuzfbnd[i, j + 1] > 0:
            mvr_pkdat.append(
                [
                    "UZF-1",
                    uzf_id_lkup[(i, j)],
                    "UZF-1",
                    uzf_id_lkup[(i, j + 1)],
                    "FACTOR",
                    1.0,
                ]
            )
            mvr_pkdat.append(
                [
                    "DRN-1",
                    uzf_id_lkup[(i, j)],
                    "UZF-1",
                    uzf_id_lkup[(i, j + 1)],
                    "FACTOR",
                    1.0,
                ]
            )
        elif iuzfbnd[i, j] > 0 and j + 1 == ncol - 1:
            mvr_pkdat.append(["UZF-1", uzf_id_lkup[(i, j)], "SFR-1", i, "FACTOR", 1.0])
            mvr_pkdat.append(["DRN-1", uzf_id_lkup[(i, j)], "SFR-1", i, "FACTOR", 1.0])

extdp = 3.0
extwc = 0.05
pet0 = 0.0
pet1 = 0.001
pet2 = 0.011
finf0 = 0.2
finf1 = 0.0
finf2 = 0.2
finf3 = 0.01
zero = 0.0
ntrail2 = 25
nsets2 = 80
uzf_spd = {
    0: [
        [0, finf0, pet0, extdp, extwc, zero, zero, zero],
        [1, finf0, pet0, extdp, extwc, zero, zero, zero],
        [2, finf0, pet0, extdp, extwc, zero, zero, zero],
        [3, finf0, pet0, extdp, extwc, zero, zero, zero],
        [4, finf0, pet0, extdp, extwc, zero, zero, zero],
        [5, finf0, pet0, extdp, extwc, zero, zero, zero],
        [6, finf0, pet0, extdp, extwc, zero, zero, zero],
        [7, finf0, pet0, extdp, extwc, zero, zero, zero],
        [8, finf1, pet0, extdp, extwc, zero, zero, zero],
        [9, finf1, pet0, extdp, extwc, zero, zero, zero],
        [10, finf1, pet0, extdp, extwc, zero, zero, zero],
        [11, finf1, pet0, extdp, extwc, zero, zero, zero],
        [12, finf1, pet0, extdp, extwc, zero, zero, zero],
        [13, finf1, pet0, extdp, extwc, zero, zero, zero],
        [14, finf1, pet0, extdp, extwc, zero, zero, zero],
        [15, finf1, pet0, extdp, extwc, zero, zero, zero],
    ],
    1: [],
    2: [],
    3: [
        [0, finf2, pet0, extdp, extwc, zero, zero, zero],
        [1, finf2, pet0, extdp, extwc, zero, zero, zero],
        [2, finf2, pet0, extdp, extwc, zero, zero, zero],
        [3, finf2, pet0, extdp, extwc, zero, zero, zero],
        [4, finf2, pet0, extdp, extwc, zero, zero, zero],
        [5, finf2, pet0, extdp, extwc, zero, zero, zero],
        [6, finf2, pet0, extdp, extwc, zero, zero, zero],
        [7, finf2, pet0, extdp, extwc, zero, zero, zero],
        [8, finf1, pet0, extdp, extwc, zero, zero, zero],
        [9, finf1, pet0, extdp, extwc, zero, zero, zero],
        [10, finf1, pet0, extdp, extwc, zero, zero, zero],
        [11, finf1, pet0, extdp, extwc, zero, zero, zero],
        [12, finf1, pet0, extdp, extwc, zero, zero, zero],
        [13, finf1, pet0, extdp, extwc, zero, zero, zero],
        [14, finf1, pet0, extdp, extwc, zero, zero, zero],
        [15, finf1, pet0, extdp, extwc, zero, zero, zero],
    ],
}

rlen = delr
rwid = 1.0
rgrd = 0.001
rtp = 26.0
rbth = 0.1
rhk = 0.0
rman = 0.02
ncon = 2
ustrf = 1.0
ndv = 0
pak_data = []
for irno in range(nrow):
    ncon = 2
    if irno in [0, nrow - 1]:
        ncon = 1
    cellid = (0, irno, ncol - 1)
    t = (irno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, rman, ncon, ustrf, ndv)
    pak_data.append(t)

con_data = []
for irno in range(nrow):
    if irno == 0:
        t = (irno, -(irno + 1))
    elif irno == nrow - 1:
        t = (irno, irno - 1)
    else:
        t = (irno, irno - 1, -(irno + 1))
    con_data.append(t)

p_data = [
    (0, "INFLOW", 0.0),
]


def build_mf6_model(idx, ws):
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]
    gwfname = "gwf-" + name

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions="NEWTON", save_flows=True
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
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=K, k33=K33)

    # aquifer storage
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=ss, sy=sy, transient=True)

    # ghb files
    flopy.mf6.ModflowGwfghb(
        gwf,
        print_flows=True,
        auxiliary="TEMPERATURE",
        stress_period_data=ghbspd,
        pname="GHB-1",
    )

    flopy.mf6.ModflowGwfdrn(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        auxiliary="DDRN",
        auxdepthname="DDRN",
        stress_period_data=drn_pkdat,
        mover=True,
        pname="DRN-1",
    )

    # transient uzf info
    uzf_obs = {
        f"{gwfname}.uzfobs": [
            ("uzf01_dpth=0.5", "water-content", "uzf01", 0.5),
            ("uzf01_dpth=1.5", "water-content", "uzf01", 1.5),  # Relies on boundnames
            ("uzf01_dpth=2.5", "water-content", "uzf01", 2.5),
            ("uzf01_dpth=3.5", "water-content", "uzf01", 3.5),
            ("uzf01_dpth=4.49", "water-content", "uzf01", 4.49),
        ]
    }
    flopy.mf6.ModflowGwfuzf(
        gwf,
        mover=True,
        print_flows=True,
        save_flows=True,
        wc_filerecord=gwfname + ".uzfwc.bin",
        simulate_et=True,
        simulate_gwseep=False,
        linear_gwet=True,
        boundnames=True,
        observations=uzf_obs,
        ntrailwaves=ntrail2,
        nwavesets=nsets2,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{gwfname}.uzf.bud",
        pname="UZF-1",
        filename=f"{gwfname}.uzf",
    )

    flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        mover=True,
        stage_filerecord=name + ".sfr.stg",
        budget_filerecord=name + ".sfr.bud",
        nreaches=nrow,
        packagedata=pak_data,
        pname="SFR-1",
        connectiondata=con_data,
        perioddata=p_data,
    )

    packages = [
        ("UZF-1",),
        ("SFR-1",),
        ("DRN-1",),
    ]
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=len(mvr_pkdat),
        budget_filerecord=f"{gwfname}.mvr.bud",
        maxpackages=len(packages),
        print_flows=True,
        packages=packages,
        perioddata=mvr_pkdat,
    )

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwfname}.oc",
    )

    # ----------
    # GWE model
    # ----------

    gwename = "gwe-" + name
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
        ktw=ktw * 86400,
        kts=kts * 86400,
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

    # Source-sink mixing package to support GHB input
    srctype = "AUX"
    auxname = "TEMPERATURE"
    pname = ["GHB-1"]
    # Inpput to SSM is: <pname> <srctype> <auxname>
    sources = [[itm, srctype, auxname] for itm in pname]

    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sources,
        pname="SSM",
        filename=f"{gwename}.ssm",
    )

    # Instantiate Streamflow Energy Transport package
    sfepackagedata = []
    for irno in range(nrow):
        t = (irno, strt_temp, K_therm_strmbed[irno], rbthcnd)
        sfepackagedata.append(t)

    sfeperioddata = {0: [(0, "INFLOW", strt_temp)]}

    flopy.mf6.modflow.ModflowGwesfe(
        gwe,
        boundnames=False,
        save_flows=True,
        print_input=False,
        print_flows=False,
        print_temperature=True,
        temperature_filerecord=gwename + ".sfe.bin",
        budget_filerecord=gwename + ".sfe.bud",
        packagedata=sfepackagedata,
        reachperioddata=sfeperioddata,
        flow_package_name="SFR-1",
        pname="SFE-1",
        filename=f"{gwename}.sfe",
    )

    # Instantiating MODFLOW 6 unsaturated zone energy transport
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
        packagedata=uze_pkdat,
        uzeperioddata=uze_perdat,
        pname="UZE-1",
        filename=f"{gwename}.uze",
    )

    # mover transport package
    fname = f"{gwename}.mve.bud"
    mve = flopy.mf6.modflow.ModflowGwemve(
        gwe, print_flows=True, budget_filerecord=fname
    )

    # Instantiate Output Control package for transport
    flopy.mf6.ModflowGweoc(
        gwe,
        temperature_filerecord=f"{gwename}.ucn",
        budget_filerecord=f"{gwename}.bud",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        temperatureprintrecord=[("COLUMNS", 3, "WIDTH", 20, "DIGITS", 8, "GENERAL")],
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

    return sim


def build_models(idx, test):
    # Start by building the MF6 model
    sim = build_mf6_model(idx, test.workspace)

    return sim, None


def check_output(idx, test):
    ws = test.workspace
    name = cases[idx]
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    # Get the model budget items
    fname = os.path.join(ws, gwfname + ".cbc")
    assert os.path.isfile(fname)
    modobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=True)

    # Get the MVR results from GWF
    fname = os.path.join(ws, gwfname + ".mvr.bud")
    assert os.path.isfile(fname)
    mvrobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=True)

    # Get the MVE results from GWE
    fname = os.path.join(ws, gwename + ".mve.bud")
    assert os.path.isfile(fname)
    mveobj = flopy.utils.CellBudgetFile(fname, precision="double", verbose=False)

    ckstpkper = mveobj.get_kstpkper()

    drnbud = modobj.get_data(text="DRN")
    drn2mvr = modobj.get_data(text="DRN-TO-MVR")
    # uzf2mvr = modobj.get_data(text="UZF-GWD TO-MVR")
    mvrdat = mvrobj.get_data(text="MOVER-FLOW")
    mvedat = mveobj.get_data(text="MVE-FLOW")

    msg0 = "Accumulated cascading runoff is not as expected"
    msg1 = "Rejected infiltration being passed to MVR where it  should not be happening"
    msg2 = (
        "The accumulated cascading runoff that is finally passed to SFR "
        "is not as expected"
    )
    msg3 = "There should be no UZF -> SFR MVR flows in rows 2 and 3"
    msg4 = (
        "There should be no cascading flow resulting from gw discharge "
        "to land surface during the first stress period, but at least "
        "one is not equal to 0.0"
    )
    msg5 = (
        "GW discharge simulated with DRN from the right-most active DRN "
        "cells should be 0 in the first stress period, but at least one "
        "is not equal to 0.0"
    )
    msg6 = (
        "GW discharge simulated by DRN for the right-most active DRN "
        "cells should not be 0.0 in stress period 3, but at least one "
        "is equal to 0.0"
    )
    msg7 = "Accumulated cascaded energy is not accumulating as expected"
    msg8 = (
        "Energy associated with rejected infiltration is being passed "
        "to MVE where it should not be happening"
    )
    msg9 = (
        "The accumulated cascading energy associated with runoff that "
        "is finally passed to SFR is not as expected"
    )
    msg10 = "There should be no UZE -> SFE MVE flows in rows 2 and 3"
    msg11 = (
        "There should be no cascading energy flow resulting from gw "
        "discharge to land surface, but at least one is not equal "
        "to 0.0"
    )
    msg12 = (
        "There should be no energy being passed by DRN from the "
        "right-most active DRN cells during the first stress period, "
        "but at least one is not equal to 0.0"
    )
    msg13 = (
        "Energy accompanying gw discharge simulated by DRN for the "
        "right-most active DRN cells should not be 0.0 in stress period "
        "3, but at least one is equal to 0.0"
    )

    accum_runoff = 0.0
    accum_energy = 0.0
    for i, current_kstpkper in enumerate(ckstpkper):
        mvr_rawdat = mvrobj.get_data(kstpkper=current_kstpkper)
        mve_rawdat = mveobj.get_data(kstpkper=current_kstpkper)
        if i == 0:
            for idx, (itm, itm_e) in enumerate(zip(mvr_rawdat, mve_rawdat)):
                # idx == 0 checks the UZF -> UZF MVR/MVE connections
                if idx == 0:
                    for ct, val in enumerate(itm):
                        if ct < 7:
                            assert np.isclose(
                                itm[ct][-1], -1 * (finf0 - vks) + accum_runoff
                            ), msg0
                            assert np.isclose(
                                itm_e[ct][-1], (finf0 - vks) * rhow * cpw + accum_energy
                            ), msg7
                            accum_runoff += -1 * (finf0 - vks)
                            accum_energy += (finf0 - vks) * rhow * cpw
                        elif ct == 8:
                            accum_runoff += -1 * (finf0 - vks)
                            accum_energy += (finf0 - vks) * rhow * cpw
                        else:
                            assert itm[ct][-1] == 0, msg1
                            assert itm_e[ct][-1] == 0, msg8

                # idx == 1 checks the UZF -> SFR MVR connections
                if idx == 1:
                    for ct, val in enumerate(itm):
                        if ct == 0:
                            assert np.isclose(itm[ct][-1], accum_runoff), msg2
                            assert np.isclose(itm_e[ct][-1], accum_energy), msg9
                        else:
                            assert itm[ct][-1] == 0, msg3
                            assert itm_e[ct][-1] == 0, msg10

                # idx == 6 checks the DRN -> UZF MVR connections
                # (groundwater discharge that cascades downhill)
                if idx == 6:
                    for ct, val in enumerate(itm):
                        assert itm[ct][-1] == 0, msg4
                        assert itm_e[ct][-1] == 0, msg11

                # idx == 7 check the DRN -> SFR MVR connections
                # (gw discharge from the lowest cell to SFR)
                if idx == 7:
                    for ct, val in enumerate(itm):
                        assert itm[ct][-1] == 0, msg5
                        assert itm_e[ct][-1] == 0, msg12

        # Expect DRN -> MVR flows in the third stress period
        if i == 2:
            accum_runoff = 0
            accum_energy = 0
            for idx, (itm, itm_e) in enumerate(zip(mvr_rawdat, mve_rawdat)):
                # idx == 0 checks the UZF -> UZF MVR connections
                if idx == 0:
                    for ct, val in enumerate(itm):
                        if ct < 7:
                            assert np.isclose(
                                itm[ct][-1], -1 * (finf0 - vks) + accum_runoff
                            ), msg0
                            assert np.isclose(
                                itm_e[ct][-1], (finf0 - vks) * rhow * cpw + accum_energy
                            ), msg7
                            accum_runoff += -1 * (finf0 - vks)
                            accum_energy += (finf0 - vks) * rhow * cpw
                        elif ct == 8:
                            accum_runoff += -1 * (finf0 - vks)
                            accum_energy += (finf0 - vks) * rhow * cpw
                        else:
                            assert itm[ct][-1] == 0, msg1
                            assert itm_e[ct][-1] == 0, msg8

                # idx == 1 checks the UZF -> SFR MVR connections
                if idx == 1:
                    for ct, val in enumerate(itm):
                        if ct == 0:
                            assert np.isclose(itm[ct][-1], accum_runoff), msg2
                            assert np.isclose(itm_e[ct][-1], accum_energy), msg9
                        else:
                            assert itm[ct][-1] == 0, msg3
                            assert itm_e[ct][-1] == 0, msg10

                # idx == 6 checks the DRN -> UZF MVR connections
                # (groundwater discharge that cascades downhill)
                if idx == 6:
                    for ct, val in enumerate(itm):
                        assert itm[ct][-1] == 0, msg4
                        assert itm_e[ct][-1] == 0, msg11

                # idx == 7 check the DRN -> SFR MVR connections
                # (gw discharge from the lowest cell to SFR)
                if idx == 7:
                    for ct, val in enumerate(itm):
                        assert itm[ct][-1] < 0, msg6
                        assert itm_e[ct][-1] > 0, msg13

    print("Finished running checks")


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
