"""
Test the advection schemes in the gwt advection package for a one-dimensional
model grid of square cells.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["nodvscale", "dvscale"]
options = [False, True]


def build_models(idx, test):
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

    c = {0: [[(0, 0, 99), 0.0000000, 100.0, 1.0]]}
    w = {0: [[(0, 0, 0), 1.0, 1e6]]}

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

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

    # ghb files
    chd = flopy.mf6.ModflowGwfghb(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        auxiliary="CONCENTRATION",
        pname="GHB-1",
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

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(
        sim,
        dependent_variable_scaling=options[idx],
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
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
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwtname}.adv")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    sourcerecarray = [
        ("WEL-1", "AUX", "CONCENTRATION"),
        ("GHB-1", "AUX", "CONCENTRATION"),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_results(test):
    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwt_names = sorted([name for name in sim.model_names if "gwt" in name])
    gwt_models = [sim.get_model(name) for name in gwt_names]

    conc = []
    for gwt in gwt_models:
        conc += gwt.output.concentration().get_data().squeeze().tolist()
    conc = np.array(conc)

    # fmt: off
    answer = np.array(
        (
            999999.9999999998, 999999.9999999993, 999999.999999999, 
            999999.9999999987, 999999.9999999942, 999999.999999957, 
            999999.9999997029, 999999.9999982064, 999999.9999904616, 
            999999.9999546632, 999999.9998050263, 999999.9992336852, 
            999999.9972244684, 999999.9906713305, 999999.9707310678, 
            999999.9138348516, 999999.7609262704, 999999.3723585819, 
            999998.43547871, 999996.2855859515, 999991.5773208102, 
            999981.7123843236, 999961.8928301098, 999923.6324732795, 
            999852.5319768365, 999725.1198872108, 999504.5989628591, 
            999135.430600611, 998536.8504703942, 997595.6348173638, 
            996158.7122537374, 994026.5045786785, 990948.1297478123, 
            986619.7481674431, 980687.3192955254, 972754.8144039324,             
            962398.4885732416, 949187.1756216577, 932707.8010452086, 
            912594.5131006193, 888559.1340068352, 860420.1536043559, 
            828127.3237138918, 791779.1151860203, 751630.8666756897, 
            708092.3216244867, 661714.3062438574, 613165.4050369008, 
            563200.4942114081, 512623.76815131743, 462249.34899546724, 
            412862.6635485552, 365185.5172132672, 319847.2497548047, 
            277363.6139511342, 238124.1830633804, 202388.27279060465, 
            170288.64812453242, 141841.7393687374, 116962.74798231329,
            95483.88541870046, 77174.03536447317, 61758.32289946244, 
            48936.36523015194, 38398.31877068738, 29838.182569768476, 
            22964.13380236391, 17505.93388555909, 13219.641597950602, 
            9890.000052793866, 7330.932693801976, 5384.599772878566, 
            3919.4436018501124, 2827.601194946717, 2021.9985541774554, 
            1433.371557988715, 1007.39149495739, 702.0135796414273, 
            485.11695773752564, 332.4656643469317, 225.9913872069924, 
            152.37954128308388, 101.92849566206363, 67.64609839669565, 
            44.546292620269156, 29.11018711319829, 18.87927997479085, 
            12.152752522872406, 7.7652221167346935, 4.925651876357927, 
            3.1020167664270706, 1.939699883174438, 1.20440811572549, 
            0.7426765112199145, 0.4548310641558006, 0.276669882183486, 
            0.1671749890963343, 0.10034923991324796, 0.05984465316349971, 
            0.03546007366567149,
            )
    )

    assert np.allclose(conc, answer), (
        "Results for the transport model are not close to the defined answer."
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_results(t),
    )
    test.run()
