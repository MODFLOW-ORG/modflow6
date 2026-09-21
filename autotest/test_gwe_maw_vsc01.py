"""
MODFLOW 6 test for MAW + MWE + VSC.

Purpose
-------
Stripped-down test to check whether temperature-dependent viscosity is
accounted for by MODFLOW 6 in the GWF solution and in the MAW well connection.
Within a single simulation, there are 3 GWF/GWE models with 1, 2, and 5 layers,
respectively.  This test checks to make sure that regardless of the layering
that is used, the calculated head inside the respective MAW wells (there is
1 injection well and 1 extraction well) are always the same for a uniform
temperature profile.  If the test ensuring that the MAW well heads are the
same regardless of layering passes, the next test is that for the extraction
well the heads in the MAW well is different between runs with and without
VSC active.  Heads in the injection MAW well are similarly checked.  There are
two tests, one for a confined setting and another for an unconfined setting.

MAW uses CONDEQN = THIEM
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["maw_ucf", "maw_cfn"]  # keep short because of MF6 file names

exe_name = "mf6.exe"
length_units = "meters"
time_units = "days"

Lx = 2500.0
Ly = 2500.0
dx = 100.0
dy = 100.0

ncol = int(Lx / dx)
nrow = int(Ly / dy)
delr = np.full(ncol, dx, dtype=float)
delc = np.full(nrow, dy, dtype=float)

nlayers = [1, 2, 5, 1, 2, 5]
ztop = 0.0
head_ini = [ztop - 50.0, ztop + 30.0]
reservoir_thickness = 100.0

ss = 1.0e-5
sy = 0.15

# Temperature / heat / VSC
T_ini_C = 60.0
T_inj_C = 20.0
T_ref_C = 20.0
mu_ref_Pa_s = 0.001  # 1 cP

# C. Voss nonlinear thermal viscosity parameters, matching many MF6 VSC examples
thermal_a2 = 10.0
thermal_a3 = 248.37
thermal_a4 = 133.15

porosity = 0.25
Cw_water = 4000.0
rho_water = 1000.0
heat_capacity_solid = 1000.0
density_solid = 2500.0
ktw = 0.6
kts = 2.5

# Time: long enough to see the cool injection temperature in the injector cell.
nper = 1
perlen = [365.0]
nstp = [1]
tsmult = [1.0]
tdis_data = list(zip(perlen, nstp, tsmult))

# Wells. Python/FloPy indices are 0-based.
wells = [
    {
        "name": "PROD",
        "role": "production",
        "row": nrow // 2,
        "col": ncol // 2 - int(ncol * 0.2),
        "q_m3_d": -2400.0,
        "T_rate_C": None,  # production does not get an imposed MWE RATE temperature
        "T_aux_C": T_ini_C,  # local reservoir temperature as AUX
    },
    {
        "name": "INJ",
        "role": "injection",
        "row": nrow // 2,
        "col": ncol // 2 + int(ncol * 0.2),
        "q_m3_d": 2400.0,
        "T_rate_C": T_inj_C,  # injection: 20 degrees through MWE RATE
        "T_aux_C": T_inj_C,  # also MAW AUX TEMP = 20, as in the Laurant logic
    },
]

# Solver values
outer_dvclose = 1.0e-6
outer_maximum = 100
under_relaxation = "NONE"
inner_maximum = 300
inner_dvclose = 1.0e-6
rcloserecord = 1.0e-5


def add_gwf_model(sim, idx, gwfname, gwename, nlay):
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    layer_thickness = reservoir_thickness / nlay
    botm_1d = np.linspace(ztop - layer_thickness, ztop - reservoir_thickness, nlay)
    botm = np.zeros((nlay, nrow, ncol), dtype=float)
    for k in range(nlay):
        botm[k, :, :] = botm_1d[k]

    idomain = np.ones((nlay, nrow, ncol), dtype=int)

    # Hydraulic parameters at reference viscosity.
    k_ref_m_d = 10.0
    k_layers = np.full((nlay, nrow, ncol), k_ref_m_d, dtype=float)
    kv_layers = np.full((nlay, nrow, ncol), k_ref_m_d, dtype=float)

    # DIS packages
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=ztop,
        botm=botm,
        idomain=idomain,
        xorigin=0.0,
        yorigin=0.0,
        pname="dis",
        filename=gwfname + ".dis",
    )

    flopy.mf6.ModflowGwfic(gwf, strt=head_ini[idx])

    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=0,
        k=k_layers,
        k33=kv_layers,
        save_specific_discharge=False,
        pname="npf",
        filename=gwfname + ".npf",
    )

    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=0,
        ss=ss,
        sy=sy,
        steady_state={0: True},
        transient={0: False},
        pname="sto",
        filename=gwfname + ".sto",
    )

    # CHD only on the outer boundary. Laurant used top CHD; in a 1-layer model
    # that would fix all cells and therefore cause dH_cell = 0.
    chd_spd = []
    for k in range(nlay):
        for r in range(nrow):
            for c in range(ncol):
                if r in [0, nrow - 1] or c in [0, ncol - 1]:
                    chd_spd.append([(k, r, c), head_ini[idx], T_ini_C])

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: chd_spd},
        auxiliary=["TEMP"],
        save_flows=True,
        pname="CHD",
        filename=gwfname + ".chd",
    )

    # MAW/MWE only in the production run
    maw_pdata = []
    maw_conn = []

    maw_spd = {0: []}

    for i, w in enumerate(wells):
        row = int(w["row"])
        col = int(w["col"])
        q = float(w["q_m3_d"])
        T_aux = float(w["T_aux_C"])

        # MAW packagedata as in Laurant: ifno, radius, bottom, strt, condeqn,
        # ngwfnodes, auxiliary TEMP, boundname.
        maw_pdata.append(
            [
                i,
                0.10,
                botm[-1, 0, 0],
                head_ini[idx],
                "THIEM",
                nlay,
                T_aux,
                w["name"],
            ]
        )

        # Connectiondata with THIEM. No custom conductance.
        # For THIEM, hk_skin and radius_skin are dummy values. The official
        # MF6 examples use -999.0/-999.0; this script follows that convention.

        for k in range(nlay):
            top_k = ztop if k == 0 else botm[k - 1, 0, 0]
            bot_k = botm[k, 0, 0]

            maw_conn.append(
                [
                    i,
                    k,
                    (k, row, col),
                    top_k,
                    bot_k,
                    -999.0,
                    -999.0,
                ]
            )

        maw_spd[0].append([i, "RATE", q])
        maw_spd[0].append([i, "AUXILIARY", "TEMP", T_aux])

    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        no_well_storage=True,
        nmawwells=len(wells),
        packagedata=maw_pdata,
        connectiondata=maw_conn,
        perioddata=maw_spd,
        auxiliary=["TEMP"],
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_head=True,
        print_flows=True,
        head_filerecord=gwfname + ".maw.hds",
        budget_filerecord=gwfname + ".maw.bud",
        pname="maw-" + str(nlay),
        filename=gwfname + ".maw",
    )

    # MAW observations following Laurant: head, MAW exchange, conductance.
    maw_obs_list = []
    for i, w in enumerate(wells):
        ifno = i + 1
        maw_obs_list.append((f"H_w{i}", "head", ifno))
        maw_obs_list.append((f"Q_w{i}_c0", "maw", ifno, 1))
        maw_obs_list.append((f"G_w{i}_c0", "conductance", ifno, 1))

    maw.obs.initialize(
        filename=gwfname + ".maw.obs",
        print_input=True,
        continuous={gwfname + ".maw.obs.csv": maw_obs_list},
    )

    if "nv" not in gwfname:
        # VSC, linked to GWE TEMP.
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=mu_ref_Pa_s,
            temperature_species_name="TEMP",
            thermal_formulation="NONLINEAR",
            thermal_a2=thermal_a2,
            thermal_a3=thermal_a3,
            thermal_a4=thermal_a4,
            viscosity_filerecord=gwfname + ".vis",
            nviscspecies=1,
            packagedata=[
                (0, 0.0, T_ref_C, gwename, "TEMP"),
            ],
            pname="VSC",
            filename=gwfname + ".vsc",
        )

    # Output control
    suffix = gwfname.strip().split("_")[-1]
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="flow_" + suffix + ".cbc",
        head_filerecord="flow_" + suffix + ".hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        filename=gwfname + ".oc",
    )

    return sim, gwf


def add_gwe_model(sim, gwename, nlay):
    gwe = flopy.mf6.ModflowGwe(
        sim,
        modelname=gwename,
        save_flows=True,
        model_nam_file=f"{gwename}.nam",
    )

    layer_thickness = reservoir_thickness / nlay
    botm_1d = np.linspace(ztop - layer_thickness, ztop - reservoir_thickness, nlay)
    botm = np.zeros((nlay, nrow, ncol), dtype=float)
    for k in range(nlay):
        botm[k, :, :] = botm_1d[k]

    flopy.mf6.ModflowGwedis(
        gwe,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=ztop,
        botm=botm,
        xorigin=0.0,
        yorigin=0.0,
        filename=gwename + ".dis",
    )

    flopy.mf6.ModflowGweic(
        gwe,
        strt=T_ini_C,
        filename=gwename + ".ic",
    )

    flopy.mf6.ModflowGweest(
        gwe,
        porosity=porosity,
        heat_capacity_water=Cw_water,
        density_water=rho_water,
        heat_capacity_solid=heat_capacity_solid,
        density_solid=density_solid,
        pname="EST",
        filename=gwename + ".est",
    )

    flopy.mf6.ModflowGwecnd(
        gwe,
        alh=10.0,
        ath1=1.0,
        ath2=1.0,
        ktw=ktw,
        kts=kts,
        pname="CND",
        filename=gwename + ".cnd",
    )

    flopy.mf6.ModflowGweadv(
        gwe, scheme="upstream", pname="ADV", filename=gwename + ".adv"
    )

    # Only CHD inflow gets TEMP through SSM. MAW temperature is handled through MWE.
    flopy.mf6.ModflowGwessm(
        gwe,
        sources=[("CHD", "AUX", "TEMP")],
        pname="SSM",
        filename=gwename + ".ssm",
    )

    # MWE. Only injection gets an imposed RATE temperature.
    mwe_pdata = []
    mwe_spd = {0: []}
    for i, w in enumerate(wells):
        mwe_pdata.append([i, T_ini_C, 0.0, 1.0])
        if float(w["q_m3_d"]) > 0.0:
            mwe_spd[0].append([i, "RATE", float(w["T_rate_C"])])

    mwe = flopy.mf6.ModflowGwemwe(
        gwe,
        flow_package_name="maw-" + str(nlay),
        flow_package_auxiliary_name="TEMP",
        packagedata=mwe_pdata,
        mweperioddata=mwe_spd,
        save_flows=True,
        print_input=True,
        print_temperature=True,
        print_flows=True,
        temperature_filerecord=gwename + ".mwe.bin",
        budget_filerecord=gwename + ".mwe.bud",
        budgetcsv_filerecord=gwename + ".mwe.bud.csv",
        pname="mwe-" + str(nlay),
        filename=gwename + ".mwe",
    )

    mwe_obs_list = []
    for i, w in enumerate(wells):
        ifno = i + 1
        mwe_obs_list.append((f"T_w{i}", "temperature", ifno))
        mwe_obs_list.append((f"E_w{i}_c0", "mwe", ifno, 1))

    mwe.obs.initialize(
        filename=gwename + ".mwe.obs",
        print_input=True,
        continuous={gwename + ".mwe.obs.csv": mwe_obs_list},
    )

    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord="energy_" + str(nlay) + "l.cbc",
        temperature_filerecord="energy_" + str(nlay) + "l.ucn",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
        filename=gwename + ".oc",
    )

    return sim, gwe


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(
        sim_name=cases[idx],
        sim_ws=ws,
        exe_name="mf6",
        version="mf6",
    )

    flopy.mf6.ModflowTdis(
        sim,
        time_units=time_units,
        nper=nper,
        perioddata=tdis_data,
    )

    # "nv" in name indicates no viscosity
    suffix = [
        "_1l",
        "_2l",
        "_5l",
        "_1lnv",
        "_2lnv",
        "_5lnv",
    ]
    for i in np.arange(len(suffix)):
        gwfname = "gwf-" + name + suffix[i]
        gwename = "gwe-" + name + suffix[i]

        # For THIEM, no skin is used.
        # c_spec = c_ref * mu_ref_pa_s / mu[idx]
        # hk_skin = c_spec if condeqn[i].upper() == "SPECIFIED" else 0.0

        nlay = nlayers[i]
        sim, gwf = add_gwf_model(sim, idx, gwfname, gwename, nlay)

        ims_gwf = flopy.mf6.ModflowIms(
            sim,
            print_option="SUMMARY",
            # complexity="MODERATE",
            outer_dvclose=outer_dvclose,
            inner_dvclose=inner_dvclose,
            rcloserecord=rcloserecord,
            outer_maximum=outer_maximum,
            inner_maximum=inner_maximum,
            linear_acceleration="BICGSTAB",
            filename=gwfname + ".ims",
        )
        sim.register_ims_package(ims_gwf, [gwf.name])

        if "nv" not in gwename:
            sim, gwe = add_gwe_model(sim, gwename, nlay)

            ims_gwe = flopy.mf6.ModflowIms(
                sim,
                print_option="SUMMARY",
                # complexity="MODERATE",
                outer_dvclose=outer_dvclose,
                inner_dvclose=inner_dvclose,
                rcloserecord=f"{rcloserecord} strict",
                outer_maximum=outer_maximum,
                inner_maximum=inner_maximum,
                linear_acceleration="BICGSTAB",
                filename=gwename + ".ims",
            )
            sim.register_ims_package(ims_gwe, [gwe.name])

            # Setup the GWF-GWE exchange
            flopy.mf6.ModflowGwfgwe(
                sim,
                exgtype="GWF6-GWE6",
                exgmnamea=gwfname,
                exgmnameb=gwename,
                filename="gwfgwe" + str(i) + ".exg",
            )

    return sim, None


def check_output(idx, test):
    ws = test.workspace
    # get all the gwf .lst files in the directory
    srchstr = ") HEADS FOR EACH CONTROL VOLUME "

    matching_files = [str(f) for f in ws.glob("*.lst") if "gwf" in f.name]

    prod_results_vsc = []
    prod_results_novsc = []
    inj_results_vsc = []
    inj_results_novsc = []

    for file in matching_files:
        with open(file) as f:
            for line in f:
                if srchstr in line:
                    fname = Path(file).name
                    digit = next(char for char in fname if char.isdigit())

                    for i in range(4):
                        line = next(f)

                    m_arr = line.strip().split()
                    if "PROD" in m_arr[0]:
                        if "nv" in fname:
                            prod_results_novsc.append(float(m_arr[-1]))
                        else:
                            prod_results_vsc.append(float(m_arr[-1]))

                    elif "INJ" in m_arr[0]:
                        if "nv" in fname:
                            inj_results_novsc.append(float(m_arr[-1]))
                        else:
                            inj_results_vsc.append(float(m_arr[-1]))

                    # read one more line
                    line = next(f)
                    m_arr = line.strip().split()
                    if "PROD" in m_arr[0]:
                        if "nv" in fname:
                            prod_results_novsc.append(float(m_arr[-1]))
                        else:
                            prod_results_vsc.append(float(m_arr[-1]))

                    elif "INJ" in m_arr[0]:
                        if "nv" in fname:
                            inj_results_novsc.append(float(m_arr[-1]))
                        else:
                            inj_results_vsc.append(float(m_arr[-1]))

    msg0 = "layering should not cause difference in results"
    assert min(prod_results_vsc) == max(prod_results_vsc), msg0
    assert min(prod_results_novsc) == max(prod_results_novsc), msg0
    assert min(inj_results_vsc) == max(inj_results_vsc), msg0
    assert min(inj_results_novsc) == max(inj_results_novsc), msg0

    msg1 = "the with and without viscosity runs should give different results"
    assert all(np.not_equal(prod_results_vsc, prod_results_novsc)), msg1
    assert all(np.not_equal(inj_results_vsc, inj_results_novsc)), msg1


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
