"""This is the reinjection problem described in the MT3D supplementary information."""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mwt_02"]


def build_models(idx, test):
    nlay = 1
    nrow = 31
    ncol = 46
    nper = 1
    delc = 10.0
    delr = 10.0
    delz = 10.0
    top = 10.0
    botm = 0.0

    perlen = [365.0]
    nstp = [100]
    tsmult = [1.0]

    Kh = 10.0
    Kv = 10.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    single_matrix = False
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

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

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    strt = 0.0
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=Kh,
        k33=Kv,
    )

    # chd files

    cinflow = 0.0
    chdlist1 = [[(0, i, 0), 1000.0 + 7.29e-4, cinflow] for i in range(nrow)]
    chdlist1 += [[(0, i, ncol - 1), 1000.0 + -4.5, 0.0] for i in range(nrow)]
    chd1 = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary=[
            ("CONCENTRATION"),
        ],
    )

    # MAW
    opth = f"{name}.maw.obs"
    # <ifno> <radius> <bottom> <strt> <condeqn> <ngwfnodes> [<aux(naux)>] [<boundname>]
    wellbottom = 0.0
    wellradius = 0.1
    ngwfnodes = 1
    concwell = 1000.0
    wellrecarray = [
        [iwell, wellradius, wellbottom, strt, "THIEM", ngwfnodes, concwell]
        for iwell in range(4)
    ]
    # <ifno> <icon> <cellid(ncelldim)> <scrn_top> <scrn_bot> <hk_skin> <radius_skin>
    wellconnectionsrecarray = [
        [0, 0, (0, 15, 15), 10.0, 0.0, 10.0, 0.1],
        [1, 0, (0, 15, 20), 10.0, 0.0, 10.0, 0.1],
        [2, 0, (0, 4, 15), 10.0, 0.0, 10.0, 0.1],
        [3, 0, (0, 26, 15), 10.0, 0.0, 10.0, 0.1],
    ]
    wellperiodrecarray = [[0, "rate", 1.0], [1, "rate", -1.0]]
    # wellperiodrecarray = [[0, 'rate', 0.], [1, 'rate', 0.]]

    mawon = True
    if mawon:
        maw = flopy.mf6.ModflowGwfmaw(
            gwf,
            mover=True,
            no_well_storage=True,
            print_input=True,
            print_head=True,
            print_flows=True,
            save_flows=True,
            packagedata=wellrecarray,
            connectiondata=wellconnectionsrecarray,
            perioddata=wellperiodrecarray,
            pname="MAW-1",
            auxiliary=["CONCENTRATION"],
        )

        packages = [
            ("maw-1",),
        ]
        perioddata = [
            ("MAW-1", 1, "MAW-1", 2, "factor", 0.5),
            ("MAW-1", 1, "MAW-1", 3, "factor", 0.5),
        ]
        mvr = flopy.mf6.ModflowGwfmvr(
            gwf,
            maxmvr=len(perioddata),
            budget_filerecord=f"{name}.mvr.bud",
            maxpackages=len(packages),
            print_flows=True,
            packages=packages,
            perioddata=perioddata,
        )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("HEAD", "LAST"),
            ("BUDGET", "LAST"),
        ],
    )

    # create gwt model
    transport = True
    if transport:
        gwtname = "gwt_" + name
        gwt = flopy.mf6.MFModel(
            sim,
            model_type="gwt6",
            modelname=gwtname,
            model_nam_file=f"{gwtname}.nam",
        )

        if not single_matrix:
            # ninner is set to 2 because the residual for the first time step
            # and 3rd inner iteration is zero, and causes a divide by zero.
            # This is a hack to get the problem to run to completion.
            nouter, ninner = 700, 2
            hclose, rclose, relax = 1e-8, 1e-6, 0.97
            imsgwt = flopy.mf6.ModflowIms(
                sim,
                print_option="ALL",
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

        # initial conditions
        ic = flopy.mf6.ModflowGwtic(gwt, strt=0.000, filename=f"{gwtname}.ic")

        # advection
        adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

        # dispersion
        diffc = 0.0
        dsp = flopy.mf6.ModflowGwtdsp(
            gwt,
            diffc=diffc,
            alh=10.0,
            ath1=3.0,
            atv=0.30,
            filename=f"{gwtname}.dsp",
        )

        # storage
        porosity = 0.30
        sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")

        mwt_obs = {
            (gwtname + ".mwt.obs.csv",): [
                ("mwt1mwt", "MWT", (0,), (0,)),
                ("mwt2mwt", "MWT", (1,), (0,)),
                ("mwt3mwt", "MWT", (2,), (0,)),
                ("mwt4mwt", "MWT", (3,), (0,)),
                ("mwt1conc", "CONCENTRATION", (0,)),
                ("mwt2conc", "CONCENTRATION", (1,)),
                ("mwt3conc", "CONCENTRATION", (2,)),
                ("mwt4conc", "CONCENTRATION", (3,)),
                ("mwt1stor", "STORAGE", (0,)),
                ("mwt2stor", "STORAGE", (1,)),
                ("mwt3stor", "STORAGE", (2,)),
                ("mwt4stor", "STORAGE", (3,)),
                ("mwt1cnst", "CONSTANT", (0,)),
                ("mwt2cnst", "CONSTANT", (1,)),
                ("mwt3cnst", "CONSTANT", (2,)),
                ("mwt4cnst", "CONSTANT", (3,)),
                ("mwt1fmvr", "FROM-MVR", (0,)),
                ("mwt2fmvr", "FROM-MVR", (1,)),
                ("mwt3fmvr", "FROM-MVR", (2,)),
                ("mwt4fmvr", "FROM-MVR", (3,)),
                ("mwt1rate", "RATE", (0,)),
                ("mwt2rate", "RATE", (1,)),
                ("mwt3rate", "RATE", (2,)),
                ("mwt4rate", "RATE", (3,)),
                ("mwt1rtmv", "RATE-TO-MVR", (0,)),
                ("mwt2rtmv", "RATE-TO-MVR", (1,)),
                ("mwt3rtmv", "RATE-TO-MVR", (2,)),
                ("mwt4rtmv", "RATE-TO-MVR", (3,)),
                ("mwt1fwrt", "FW-RATE", (0,)),
                ("mwt2fwrt", "FW-RATE", (1,)),
                ("mwt3fwrt", "FW-RATE", (2,)),
                ("mwt4fwrt", "FW-RATE", (3,)),
                ("mwt1frtm", "FW-RATE-TO-MVR", (0,)),
                ("mwt2frtm", "FW-RATE-TO-MVR", (1,)),
                ("mwt3frtm", "FW-RATE-TO-MVR", (2,)),
                ("mwt4frtm", "FW-RATE-TO-MVR", (3,)),
            ],
        }
        # append additional obs attributes to obs dictionary
        mwt_obs["digits"] = 12
        mwt_obs["print_input"] = True
        mwt_obs["filename"] = gwtname + ".mwt.obs"

        mwtpackagedata = [
            (0, 0.0, 99.0, 999.0, "inject"),
            (1, 0.0, 99.0, 999.0, "extract"),
            (2, 0.0, 99.0, 999.0, "reinject1"),
            (3, 0.0, 99.0, 999.0, "reinject2"),
        ]
        mwtperioddata = [
            (0, "RATE", 1000.0),
            (0, "CONCENTRATION", 1000.0),
            (0, "STATUS", "ACTIVE"),
            (1, "STATUS", "ACTIVE"),
            (2, "STATUS", "ACTIVE"),
            (3, "STATUS", "ACTIVE"),
        ]

        mwton = True
        if mwton:
            mwt = flopy.mf6.modflow.ModflowGwtmwt(
                gwt,
                boundnames=True,
                save_flows=True,
                print_input=True,
                print_flows=True,
                print_concentration=True,
                concentration_filerecord=gwtname + ".mwt.bin",
                budget_filerecord=gwtname + ".mwt.bud",
                packagedata=mwtpackagedata,
                mwtperioddata=mwtperioddata,
                observations=mwt_obs,
                pname="MAW-1",
                auxiliary=["aux1", "aux2"],
            )

            mvt = flopy.mf6.modflow.ModflowGwtmvt(gwt, print_flows=True)

        # sources
        sourcerecarray = [
            ("CHD-1", "AUX", "CONCENTRATION"),
        ]
        if not mwton:
            sourcerecarray.append(("MAW-1", "AUX", "CONCENTRATION"))
        ssm = flopy.mf6.ModflowGwtssm(
            gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
        )

        # output control
        oc = flopy.mf6.ModflowGwtoc(
            gwt,
            budget_filerecord=f"{gwtname}.cbc",
            concentration_filerecord=f"{gwtname}.ucn",
            concentrationprintrecord=[
                ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord=[
                ("CONCENTRATION", "ALL"),
            ],
            printrecord=[
                ("CONCENTRATION", "ALL"),
                ("BUDGET", "ALL"),
            ],
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


def plot_output(idx, test):
    name = test.name
    ws = test.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    fname = gwtname + ".ucn"
    fname = os.path.join(ws, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")  # , precision='double')
    conc = cobj.get_data()

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pmv.plot_grid()
    pmv.plot_bc(ftype="MAW")
    pmv.plot_bc(ftype="CHD")
    a = np.ma.masked_less(conc, 0.01)
    pa = pmv.plot_array(a, cmap="jet", alpha=0.25)
    cs = pmv.contour_array(conc, levels=[0.01, 0.1, 1, 10, 100], colors="y")
    plt.colorbar(pa, shrink=0.5)
    plt.draw()
    fname = os.path.join(ws, gwtname + ".png")
    plt.savefig(fname)

    return


def check_output(idx, test):
    # ensure concentrations were saved
    name = cases[idx]
    gwtname = "gwt_" + name
    fname = gwtname + ".mwt.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    # load and check the well concentrations
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    cmwt = cobj.get_data().flatten()
    answer = np.array([999.98345654, 18.67908708, 15.9497297, 15.94973001])
    assert np.allclose(cmwt, answer), f"{cmwt} {answer}"

    # make sure concentrations can be loaded
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    cobj = flopy.utils.HeadFile(fname, text="CONCENTRATION")
    caq = cobj.get_alldata()

    # make sure observations can be loaded
    fpth = os.path.join(test.workspace, gwtname + ".mwt.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'
    res = [
        tc["MWT1CONC"][-1],
        tc["MWT2CONC"][-1],
        tc["MWT3CONC"][-1],
        tc["MWT4CONC"][-1],
    ]
    res = np.array(res)
    answer = np.array([999.98345654, 18.67908708, 15.9497297, 15.94973001])
    assert np.allclose(res, answer), f"{res} {answer}"

    res = tc["MWT1RATE"]
    answer = np.ones(res.shape) * 1000.0
    assert np.allclose(res, answer), f"{res} {answer}"


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
    )
    test.run()
