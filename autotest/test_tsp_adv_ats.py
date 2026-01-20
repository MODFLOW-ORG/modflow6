"""
Test the ATS_PERCEL option in the ADV Package for controlling the time step.
This problem is based on test_gwt_adv01, which is just 1D transport using
a well injecting with a concentration of 1 and a constant-head boundary on the
right side allowing water to leave.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwtadvats01", "gweadvats01"]

nlay, nrow, ncol = 1, 1, 100
nper = 1
perlen = [5.0]
nstp = [200]
tsmult = [1.0]
steady = [True]
delr = 1.0
delc = 1.0
top = 1.0
botm = [0.0]
strt = 1.0
hk = 1.0
laytyp = 0
rhos = 2700.0
rhow = 1000.0
Cps = 703.7
Cpw = 4183.0

c = {0: [[(0, 0, 99), 0.0000000]]}
w = {0: [[(0, 0, 0), 1.0, 1.0]]}

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def add_solute_tsp(sim, name, gwtname):
    # create gwt model
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        ats_outer_maximum_fraction=0.0,
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

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(
        gwt, scheme="upstream", ats_percel=0.5, filename=f"{gwtname}.adv"
    )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        "conc_obs.csv": [
            ("(1-1-10)", "CONCENTRATION", (0, 0, 9)),
            ("(1-1-50)", "CONCENTRATION", (0, 0, 49)),
        ],
        "flow_obs.csv": [
            ("c10-c11", "FLOW-JA-FACE", (0, 0, 9), (0, 0, 10)),
            ("c50-c51", "FLOW-JA-FACE", (0, 0, 49), (0, 0, 50)),
            ("c99-c100", "FLOW-JA-FACE", (0, 0, 98), (0, 0, 99)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname="conc_obs",
        filename=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    return sim


def add_energy_tsp(sim, name, gwename):
    # create gwt model
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    gwe.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
    imsgwe = flopy.mf6.ModflowIms(
        sim,
        ats_outer_maximum_fraction=0.0,
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

    dis = flopy.mf6.ModflowGwtdis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwename}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGweic(gwe, strt=0.0, filename=f"{gwename}.ic")

    # advection
    adv = flopy.mf6.ModflowGweadv(
        gwe, scheme="upstream", ats_percel=0.5, filename=f"{gwename}.adv"
    )

    # mass storage and transfer
    est = flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.1,
        heat_capacity_water=Cpw,
        density_water=rhow,
        heat_capacity_solid=Cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwessm(
        gwe, sources=sourcerecarray, filename=f"{gwename}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("TEMPERATURE", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        "temp_obs.csv": [
            ("(1-1-10)", "TEMPERATURE", (0, 0, 9)),
            ("(1-1-50)", "TEMPERATURE", (0, 0, 49)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwe,
        pname="temp_obs",
        filename=f"{gwename}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    return sim


def build_models(idx, test):
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # set dt0, dtmin, dtmax, dtadj, dtfailadj
    dt0 = 0.01
    dtmin = 1.0e-5
    dtmax = perlen
    dtadj = 2.0
    dtfailadj = 5.0
    ats_filerecord = name + ".ats"
    atsperiod = [(0, dt0, dtmin, dtmax[i], dtadj, dtfailadj) for i in range(nper)]
    tdis.ats.initialize(
        maxats=len(atsperiod),
        perioddata=atsperiod,
        filename=ats_filerecord,
    )

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        ats_outer_maximum_fraction=0.0,
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
    sim.register_ims_package(imsgwf, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=laytyp,
        k=hk,
        k33=hk,
        save_specific_discharge=True,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        pname="CHD-1",
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(w),
        stress_period_data=w,
        save_flows=False,
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    if "gwt" in name:
        gwtname = "gwt_" + name
        sim = add_solute_tsp(sim, name, gwtname)

        # GWF GWT exchange
        gwfgwt = flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            filename=f"{name}.gwfgwt",
        )

    elif "gwe" in name:
        gwename = "gwe_" + name
        sim = add_energy_tsp(sim, name, gwename)

        # GWF GWE exchange
        gwfgwe = flopy.mf6.ModflowGwfgwe(
            sim,
            exgtype="GWF6-GWE6",
            exgmnamea=gwfname,
            exgmnameb=gwename,
            filename=f"{name}.gwfgwe",
        )

    return sim, None


def check_gwt_output(gwtname, test):
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    print("made it here?")
    print(fpth)
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
        times = cobj.times
    except:
        assert False, f'could not load data from "{fpth}"'

    # Verify number of time steps; there should be 101 of them.
    assert len(times) == 101, f"Found {len(times)} time steps.  Expecting 101."

    # This is the answer to this problem.  These concentrations are for
    # time step 200.
    cres1 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999997e-01,
                9.99999991e-01,
                9.99999971e-01,
                9.99999914e-01,
                9.99999761e-01,
                9.99999372e-01,
                9.99998435e-01,
                9.99996286e-01,
                9.99991577e-01,
                9.99981712e-01,
                9.99961893e-01,
                9.99923632e-01,
                9.99852532e-01,
                9.99725120e-01,
                9.99504599e-01,
                9.99135431e-01,
                9.98536850e-01,
                9.97595635e-01,
                9.96158712e-01,
                9.94026505e-01,
                9.90948130e-01,
                9.86619748e-01,
                9.80687319e-01,
                9.72754814e-01,
                9.62398489e-01,
                9.49187176e-01,
                9.32707801e-01,
                9.12594513e-01,
                8.88559134e-01,
                8.60420154e-01,
                8.28127324e-01,
                7.91779115e-01,
                7.51630867e-01,
                7.08092322e-01,
                6.61714306e-01,
                6.13165405e-01,
                5.63200494e-01,
                5.12623768e-01,
                4.62249349e-01,
                4.12862664e-01,
                3.65185517e-01,
                3.19847250e-01,
                2.77363614e-01,
                2.38124183e-01,
                2.02388273e-01,
                1.70288648e-01,
                1.41841739e-01,
                1.16962748e-01,
                9.54838854e-02,
                7.71740354e-02,
                6.17583229e-02,
                4.89363652e-02,
                3.83983188e-02,
                2.98381826e-02,
                2.29641338e-02,
                1.75059339e-02,
                1.32196416e-02,
                9.89000005e-03,
                7.33093269e-03,
                5.38459977e-03,
                3.91944360e-03,
                2.82760119e-03,
                2.02199855e-03,
                1.43337156e-03,
                1.00739149e-03,
                7.02013580e-04,
                4.85116958e-04,
                3.32465664e-04,
                2.25991387e-04,
                1.52379541e-04,
                1.01928496e-04,
                6.76460984e-05,
                4.45462926e-05,
                2.91101871e-05,
                1.88792800e-05,
                1.21527525e-05,
                7.76522212e-06,
                4.92565188e-06,
                3.10201677e-06,
                1.93969988e-06,
                1.20440812e-06,
                7.42676511e-07,
                4.54831064e-07,
                2.76669882e-07,
                1.67174989e-07,
                1.00349240e-07,
                5.98446532e-08,
                3.54600737e-08,
            ]
        ]
    ]
    cres1 = np.array(cres1)
    diff = cres1 - conc
    print(f"{diff.min()=} {diff.max()=}")
    assert np.allclose(cres1, conc, atol=0.03), (
        "simulated concentrations do not match with known solution: "
        f"{diff.min()=} {diff.max()=}"
    )


def check_gwe_output(gwename, test):
    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    try:
        tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
        temp = tobj.get_data()
        times = tobj.times
    except:
        assert False, f'could not load data from "{fpth}"'

    # Verify number of time steps; there should be 101 of them.
    assert len(times) == 101, f"Found {len(times)} time steps.  Expecting 101."

    # This is the answer to this problem.  These temperatures are for
    # time step 200.
    tres1 = [
        [
            [
                9.999151944e-01,
                9.991561416e-01,
                9.957253385e-01,
                9.852855303e-01,
                9.612267633e-01,
                9.164423891e-01,
                8.463063960e-01,
                7.512651490e-01,
                6.375137933e-01,
                5.153686634e-01,
                3.962365131e-01,
                2.896398445e-01,
                2.014154568e-01,
                1.334081536e-01,
                8.429617551e-02,
                5.090190174e-02,
                2.942801589e-02,
                1.631898616e-02,
                8.696001156e-03,
                4.460691183e-03,
                2.206323719e-03,
                1.053931687e-03,
                4.869516105e-04,
                2.179244364e-04,
                9.459177996e-05,
                3.987257926e-05,
                1.634113891e-05,
                6.518691628e-06,
                2.533746783e-06,
                9.605469452e-07,
                3.554919904e-07,
                1.285514974e-07,
                4.545895740e-08,
                1.573243379e-08,
                5.332454436e-09,
                1.771400518e-09,
                5.771017791e-10,
                1.845044867e-10,
                5.792144119e-11,
                1.786470543e-11,
                5.416398110e-12,
                1.615126084e-12,
                4.739099550e-13,
                1.368930367e-13,
                3.894541782e-14,
                1.091705682e-14,
                3.016515444e-15,
                8.219118876e-16,
                2.209156368e-16,
                5.859548158e-17,
                1.534217316e-17,
                3.966772505e-18,
                1.013099873e-18,
                2.556602974e-19,
                6.376711964e-20,
                1.572441116e-20,
                3.834538402e-21,
                9.249656326e-22,
                2.207600521e-22,
                5.214377130e-23,
                1.219192845e-23,
                2.822460104e-24,
                6.470870604e-25,
                1.469494511e-25,
                3.306217086e-26,
                7.371193120e-27,
                1.628805914e-27,
                3.567840659e-28,
                7.748596208e-29,
                1.668770248e-29,
                3.564496198e-30,
                7.552622494e-31,
                1.587680860e-31,
                3.311770283e-32,
                6.855711936e-33,
                1.408647074e-33,
                2.873222323e-34,
                5.818513572e-35,
                1.170004097e-35,
                2.336416973e-36,
                4.633984263e-37,
                9.129607438e-38,
                1.786873029e-38,
                3.474790692e-39,
                6.714365487e-40,
                1.289345817e-40,
                2.460749895e-41,
                4.668135881e-42,
                8.803201667e-43,
                1.650444997e-43,
                3.076562467e-44,
                5.702623597e-45,
                1.051155543e-45,
                1.926995214e-46,
                3.513603911e-47,
                6.372645060e-48,
                1.149785438e-48,
                2.063852098e-49,
                3.685859442e-50,
                6.549834922e-51,
            ]
        ]
    ]
    tres1 = np.array(tres1)
    diff = tres1 - temp
    print(f"{diff.min()=} {diff.max()=}")
    assert np.allclose(tres1, temp, atol=0.03), (
        "simulated temperatures do not match with known solution: "
        f"{diff.min()=} {diff.max()=}"
    )


def check_output(idx, test):
    name = cases[idx]
    if "gwt" in name:
        gwtname = "gwt_" + name
        check_gwt_output(gwtname, test)
    elif "gwe" in name:
        gwename = "gwe_" + name
        check_gwe_output(gwename, test)


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
