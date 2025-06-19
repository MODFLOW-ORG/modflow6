import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "moc3d01a",
    "moc3d01b",
    "moc3d01c",
    "moc3d01d",
    "moc3d01e",
    "moc3d01f",
    "moc3d01g",
    "moc3d01h",
]
diffc = [0.0, 0.01, 0.0, 0.1, 0.0, 0.0, 0.0, 0]
alphal = [0.1, 0.0, 1.0, 0.0, 0.1, 0.1, 0.1, 0.1]
retardation = [None, None, None, None, 40.0, 4.0, 2.0, None]
perlens = 4 * [120.0] + 3 * [240.0] + [120.0]
decay = 7 * [None] + [0.01]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 122, 1
    nper = 1
    perlen = perlens[idx]  # [120.]
    perlen = [perlen]
    nstp = [240]
    tsmult = [1.0]
    steady = [True]
    delr = 0.1
    delc = 0.1
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 0.01
    laytyp = 0
    # ss = 0.
    # sy = 0.1

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-8, 1e-6, 1.0

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
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
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
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
    )
    # storage
    # sto = flopy.mf6.ModflowGwfsto(gwf, save_flows=False,
    #                              iconvert=laytyp[idx],
    #                              ss=ss[idx], sy=sy[idx],
    #                              steady_state={0: True, 2: True},
    #                              transient={1: True})

    # chd files
    c = {0: [[(0, 121, 0), 0.0000000]]}
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=c, save_flows=False, pname="CHD-1"
    )

    # wel files
    w = {0: [[(0, 0, 0), 0.001, 1.0]]}
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
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
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

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
    strt = np.zeros((nlay, nrow, ncol))
    strt[0, 0, 0] = 0.0
    ic = flopy.mf6.ModflowGwtic(gwt, strt=strt, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="tvd", filename=f"{gwtname}.adv")

    # dispersion
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        diffc=diffc[idx],
        alh=alphal[idx],
        alv=alphal[idx],
        ath1=0.0,
        atv=0.0,
        filename=f"{gwtname}.dsp",
    )

    # constant concentration
    # cncs = {0: [[(0, 0, 0), 1.0]]}
    # cnc = flopy.mf6.ModflowGwtcnc(gwt, maxbound=len(cncs),
    #                              stress_period_data=cncs,
    #                              save_flows=False,
    #                              pname='CNC-1')

    # storage
    porosity = 0.1

    rtd = retardation[idx]
    sorption = None
    kd = None
    if rtd is not None:
        rhob = 1.0
        kd = (rtd - 1.0) * porosity / rhob

    decay_rate = decay[idx]
    first_order_decay = False
    if decay_rate is not None:
        first_order_decay = True

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        first_order_decay=first_order_decay,
        decay=decay_rate,
        decay_sorbed=decay_rate,
        sorption=sorption,
        distcoef=kd,
    )

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
        saverecord=[("CONCENTRATION", "ALL")],
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


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        station = [(0, 0, 0), (0, 40, 0), (0, 110, 0)]
        tssim = cobj.get_ts(station)[::10]
    except:
        assert False, f'could not load data from "{fpth}"'

    tsresab = [
        [5.00000000e-01, 3.03340361e-01, 7.99861483e-18, 6.87424430e-43],
        [5.50000000e00, 9.40912181e-01, 2.47486965e-10, -2.35610314e-34],
        [1.05000000e01, 9.89534465e-01, 1.20547721e-06, 5.05194948e-26],
        [1.55000000e01, 9.97735140e-01, 1.83398462e-04, 8.49328574e-22],
        [2.05000000e01, 9.99460043e-01, 4.34064097e-03, 1.32141059e-18],
        [2.55000000e01, 9.99863529e-01, 3.33283202e-02, 4.85619862e-16],
        [3.05000000e01, 9.99964115e-01, 1.24494336e-01, 6.07985010e-14],
        [3.55000000e01, 9.99990289e-01, 2.89407918e-01, 1.05411307e-11],
        [4.05000000e01, 9.99997314e-01, 4.90869811e-01, 7.82976879e-10],
        [4.55000000e01, 9.99999244e-01, 6.76144792e-01, 2.90810581e-08],
        [5.05000000e01, 9.99999783e-01, 8.14370893e-01, 6.05340027e-07],
        [5.55000000e01, 9.99999937e-01, 9.02516345e-01, 7.82329321e-06],
        [6.05000000e01, 9.99999982e-01, 9.52377775e-01, 6.74435581e-05],
        [6.55000000e01, 9.99999995e-01, 9.78079028e-01, 4.12193060e-04],
        [7.05000000e01, 9.99999998e-01, 9.90391770e-01, 1.87927427e-03],
        [7.55000000e01, 1.00000000e00, 9.95955737e-01, 6.66985043e-03],
        [8.05000000e01, 1.00000000e00, 9.98354094e-01, 1.91012548e-02],
        [8.55000000e01, 1.00000000e00, 9.99348795e-01, 4.55057277e-02],
        [9.05000000e01, 1.00000000e00, 9.99748429e-01, 9.25739475e-02],
        [9.55000000e01, 1.00000000e00, 9.99904760e-01, 1.64507712e-01],
        [1.00500000e02, 1.00000000e00, 9.99964562e-01, 2.60492805e-01],
        [1.05500000e02, 1.00000000e00, 9.99987008e-01, 3.74082729e-01],
        [1.10500000e02, 1.00000000e00, 9.99995298e-01, 4.94921031e-01],
        [1.15500000e02, 1.00000000e00, 9.99998317e-01, 6.11836684e-01],
    ]

    tsrescd = [
        [5.00000000e-01, 1.63903900e-01, 1.98576084e-08, 2.23155914e-20],
        [5.50000000e00, 5.71069200e-01, 7.70702960e-04, 1.32230703e-12],
        [1.05000000e01, 7.10756331e-01, 1.50144518e-02, 4.41433414e-09],
        [1.55000000e01, 7.90734103e-01, 6.02106964e-02, 5.74646888e-07],
        [2.05000000e01, 8.42900447e-01, 1.32630784e-01, 1.43891695e-05],
        [2.55000000e01, 8.79242901e-01, 2.19274259e-01, 1.38109607e-04],
        [3.05000000e01, 9.05612271e-01, 3.09294834e-01, 7.25311048e-04],
        [3.55000000e01, 9.25287101e-01, 3.96016251e-01, 2.55099736e-03],
        [4.05000000e01, 9.40271299e-01, 4.75983392e-01, 6.79158982e-03],
        [4.55000000e01, 9.51863681e-01, 5.47772190e-01, 1.47977972e-02],
        [5.05000000e01, 9.60945331e-01, 6.11125470e-01, 2.77897366e-02],
        [5.55000000e01, 9.68133171e-01, 6.66408631e-01, 4.66123366e-02],
        [6.05000000e01, 9.73870806e-01, 7.14286611e-01, 7.16152557e-02],
        [6.55000000e01, 9.78484019e-01, 7.55539173e-01, 1.02652526e-01],
        [7.05000000e01, 9.82216253e-01, 7.90959134e-01, 1.39163570e-01],
        [7.55000000e01, 9.85252079e-01, 8.21299240e-01, 1.80292746e-01],
        [8.05000000e01, 9.87733167e-01, 8.47247169e-01, 2.25014235e-01],
        [8.55000000e01, 9.89769412e-01, 8.69416424e-01, 2.72242050e-01],
        [9.05000000e01, 9.91446840e-01, 8.88345988e-01, 3.20915892e-01],
        [9.55000000e01, 9.92833338e-01, 9.04504498e-01, 3.70061119e-01],
        [1.00500000e02, 9.93982855e-01, 9.18296556e-01, 4.18825304e-01],
        [1.05500000e02, 9.94938532e-01, 9.30069757e-01, 4.66495807e-01],
        [1.10500000e02, 9.95735063e-01, 9.40121702e-01, 5.12503175e-01],
        [1.15500000e02, 9.96400487e-01, 9.48706581e-01, 5.56414806e-01],
    ]
    tsresab = np.array(tsresab)
    tsrescd = np.array(tsrescd)

    tsreslist = [tsresab, tsresab, tsrescd, tsrescd, None, None, None, None]
    tsres = tsreslist[idx]
    if tsres is not None:
        assert np.allclose(tsres, tssim), (
            "simulated concentrations do not match with known solution."
        )


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
