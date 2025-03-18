"""
Test isotherms.
"""

import os

import flopy
import numpy as np
import pytest
from binary_util import write_budget, write_head
from flopy.utils.gridutil import uniform_flow_field
from framework import TestFramework

cases = ["mst05a", "mst05b"]
isotherm = ["freundlich", "langmuir"]
distcoef = [0.3, 100.0]
sp2 = [0.7, 0.003]
xmax_plot = [1500, 500]
ymax_plot = [0.5, 1.0]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 101
    perlen = [160.0, 1340.0]
    nper = len(perlen)
    tslength = 10.0
    nstp = [p / tslength for p in perlen]
    tsmult = nper * [1.0]
    delr = 0.16
    delc = 0.16
    top = 1.0
    botm = 0.0
    velocity = 0.1
    porosity = 0.37
    bulk_density = 1.587
    dispersivity = 1.0
    source_concentration = 0.05
    specific_discharge = velocity * porosity
    inflow_rate = specific_discharge * delc * (top - botm)

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
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="SECONDS", nper=nper, perioddata=tdis_rc
    )

    # create gwt model
    gwtname = "gwt_" + name
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
    )
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt, xt3d_off=True, alh=dispersivity, ath1=dispersivity
    )
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        sorption=isotherm[idx],
        porosity=porosity,
        bulk_density=bulk_density,
        distcoef=distcoef[idx],
        sp2=sp2[idx],
    )

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # create a heads file with head equal top
    fname = os.path.join(ws, "myhead.hds")
    with open(fname, "wb") as fbin:
        kstp = 0
        totim = 0
        for kper in range(nper):
            totim += perlen[kper]
            write_head(
                fbin,
                top * np.ones((nrow, ncol)),
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
            )

    # create a budget file
    qx = specific_discharge
    qy = 0.0
    qz = 0.0
    shape = (nlay, nrow, ncol)
    spdis, flowja = uniform_flow_field(qx, qy, qz, shape, delr=delr, delc=delc)
    dt = np.dtype(
        [
            ("ID1", np.int32),
            ("ID2", np.int32),
            ("FLOW", np.float64),
            ("CONCENTRATION", np.float64),
        ]
    )
    wel = [
        np.array([(0 + 1, 0 + 1, inflow_rate, source_concentration)], dtype=dt),
        np.array([(0 + 1, 0 + 1, inflow_rate, 0.0)], dtype=dt),
    ]
    chd = np.array([(ncol - 1 + 1, ncol - 1 + 1, -inflow_rate, 0.0)], dtype=dt)

    dt = np.dtype(
        [
            ("ID1", np.int32),
            ("ID2", np.int32),
            ("FLOW", np.float64),
            ("SATURATION", np.float64),
        ]
    )
    sat = np.array([(i, i, 0.0, 1.0) for i in range(nlay * nrow * ncol)], dtype=dt)

    fname = os.path.join(ws, "mybudget.bud")
    with open(fname, "wb") as fbin:
        kstp = 0
        totim = 0
        for kper in range(nper):
            totim += perlen[kper]
            delt = perlen[kper] / nstp[kper]
            write_budget(
                fbin,
                flowja,
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
                delt=delt,
            )
            write_budget(
                fbin,
                spdis,
                text="      DATA-SPDIS",
                imeth=6,
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
                delt=delt,
            )
            write_budget(
                fbin,
                sat,
                text="        DATA-SAT",
                imeth=6,
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
                delt=delt,
            )
            write_budget(
                fbin,
                wel[kper],
                text="             WEL",
                imeth=6,
                text2id2="           WEL-1",
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
                delt=delt,
            )
            write_budget(
                fbin,
                chd,
                text="             CHD",
                imeth=6,
                text2id2="           CHD-1",
                kstp=kstp + 1,
                kper=kper + 1,
                pertim=perlen[kper],
                totim=totim,
                delt=delt,
            )
    fbin.close()

    # flow model interface
    packagedata = [
        ("GWFBUDGET", "mybudget.bud", None),
        ("GWFHEAD", "myhead.hds", None),
    ]
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=packagedata)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        "conc_obs.csv": [
            ("X008", "CONCENTRATION", (0, 0, 50)),
        ]
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname="conc_obs",
        filename=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    fpth = os.path.join(test.workspace, "conc_obs.csv")
    try:
        obs = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    cnorm = obs["X008"] / 0.05
    cnorm_max = [0.32842034, 0.875391418]
    msg = f"{cnorm_max[idx]} /= {cnorm.max()}"
    assert np.allclose(cnorm_max[idx], cnorm.max(), atol=0.001), msg

    savefig = False
    if savefig:
        import matplotlib.pyplot as plt

        fig = plt.figure()
        plt.plot(obs["time"], obs["X008"] / 0.05, "bo-")
        plt.xlim(0, xmax_plot[idx])
        plt.ylim(0, ymax_plot[idx])
        plt.xlabel("Time, in seconds")
        plt.ylabel("Normalized Concentration")
        plt.title(isotherm[idx])
        fname = os.path.join(test.workspace, "results.png")
        plt.savefig(fname)


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
