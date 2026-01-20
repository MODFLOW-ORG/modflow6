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
        [5.00000000e-01, 2.83603277e-01, 1.98913375e-16, 4.55149741e-41],
        [5.50000000e00, 9.17951825e-01, 2.69455937e-10, -2.07960785e-31],
        [1.05000000e01, 9.84228070e-01, 1.28040427e-06, 2.48438252e-26],
        [1.55000000e01, 9.96457822e-01, 1.90868536e-04, 3.41305818e-22],
        [2.05000000e01, 9.99136905e-01, 4.44854016e-03, 4.50363732e-19],
        [2.55000000e01, 9.99778638e-01, 3.37898628e-02, 1.68992365e-16],
        [3.05000000e01, 9.99941176e-01, 1.25336227e-01, 6.43299718e-14],
        [3.55000000e01, 9.99983954e-01, 2.90158229e-01, 1.10351780e-11],
        [4.05000000e01, 9.99995534e-01, 4.91078611e-01, 8.14269315e-10],
        [4.55000000e01, 9.99998737e-01, 6.75816022e-01, 3.00545939e-08],
        [5.05000000e01, 9.99999636e-01, 8.13798122e-01, 6.22095988e-07],
        [5.55000000e01, 9.99999894e-01, 9.01968239e-01, 7.99960358e-06],
        [6.05000000e01, 9.99999969e-01, 9.51974734e-01, 6.86618749e-05],
        [6.55000000e01, 9.99999991e-01, 9.77827407e-01, 4.18059967e-04],
        [7.05000000e01, 9.99999997e-01, 9.90251754e-01, 1.89993945e-03],
        [7.55000000e01, 9.99999999e-01, 9.95884247e-01, 6.72520577e-03],
        [8.05000000e01, 1.00000000e00, 9.98319951e-01, 1.92175257e-02],
        [8.55000000e01, 1.00000000e00, 9.99333379e-01, 4.57014056e-02],
        [9.05000000e01, 1.00000000e00, 9.99741747e-01, 9.28423497e-02],
        [9.55000000e01, 1.00000000e00, 9.99901965e-01, 1.64807939e-01],
        [1.00500000e02, 1.00000000e00, 9.99963427e-01, 2.60760778e-01],
        [1.05500000e02, 1.00000000e00, 9.99986559e-01, 3.74256690e-01],
        [1.10500000e02, 1.00000000e00, 9.99995124e-01, 4.94966053e-01],
        [1.15500000e02, 1.00000000e00, 9.99998251e-01, 6.11754824e-01],
    ]

    tsrescd = [
        [5.00000000e-01, 1.62245199e-01, 1.99581015e-08, 1.96866573e-20],
        [5.50000000e00, 5.69190185e-01, 7.72268573e-04, 1.32726540e-12],
        [1.05000000e01, 7.09416269e-01, 1.50275044e-02, 4.42507587e-09],
        [1.55000000e01, 7.89741082e-01, 6.02353184e-02, 5.75599816e-07],
        [2.05000000e01, 8.42144507e-01, 1.32658034e-01, 1.44051847e-05],
        [2.55000000e01, 8.78656387e-01, 2.19297521e-01, 1.38208708e-04],
        [3.05000000e01, 9.05151032e-01, 3.09311460e-01, 7.25613436e-04],
        [3.55000000e01, 9.24920151e-01, 3.96026115e-01, 2.55150051e-03],
        [4.05000000e01, 9.39976345e-01, 4.75987512e-01, 6.79179030e-03],
        [4.55000000e01, 9.51625144e-01, 5.47771861e-01, 1.47963272e-02],
        [5.05000000e01, 9.60751243e-01, 6.11121931e-01, 2.77842830e-02],
        [5.55000000e01, 9.67974420e-01, 6.66402923e-01, 4.65998752e-02],
        [6.05000000e01, 9.73740367e-01, 7.14279555e-01, 7.15925045e-02],
        [6.55000000e01, 9.78376414e-01, 7.55531394e-01, 1.02616451e-01],
        [7.05000000e01, 9.82127171e-01, 7.90951085e-01, 1.39111803e-01],
        [7.55000000e01, 9.85178098e-01, 8.21291242e-01, 1.80223847e-01],
        [8.05000000e01, 9.87671554e-01, 8.47239438e-01, 2.24927782e-01],
        [8.55000000e01, 9.89717967e-01, 8.69409099e-01, 2.72138588e-01],
        [9.05000000e01, 9.91403786e-01, 8.88339150e-01, 3.20796786e-01],
        [9.55000000e01, 9.92797229e-01, 9.04498187e-01, 3.69928359e-01],
        [1.00500000e02, 9.93952511e-01, 9.18290783e-01, 4.18681292e-01],
        [1.05500000e02, 9.94912987e-01, 9.30064513e-01, 4.66343157e-01],
        [1.10500000e02, 9.95713522e-01, 9.40116964e-01, 5.12344545e-01],
        [1.15500000e02, 9.96382294e-01, 9.48702318e-01, 5.56252766e-01],
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
