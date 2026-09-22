"""
Test for the Richards based unsaturated zone package UZR.
"""

import os

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pytest
from framework import TestFramework
from gwf_test_utils import PLOT_UZR_TESTS, get_uzr_soil_data
from modflow_devtools.misc import is_in_ci

cases = ["uzr-trs-pic", "uzr-trs-nwt", "uzr-std-pic", "uzr-std-nwt"]
steadystate = [False, False, True, True]
newton = [False, True, False, True]


def build_models(idx, test):
    column_height = 42.0
    nlay, nrow, ncol = 42, 1, 1
    nper = 10
    perlen = nper * [360.0]  # s
    tsmult = nper * [1.0]
    nstp = nper * [1]
    if steadystate[idx]:
        nper = 1
        perlen = [sum(perlen)]
        tsmult = [1.0]
        nstp = [1]

    delr = 10.0  # cm
    delc = 10.0
    delz = column_height / nlay
    top = delz
    laytyp = 0
    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
    hk = 0.00944  # cm/s

    # saturated lowest cell:
    hp_upper = -20.7
    hp_lower = -61.5
    h_upper = hp_upper + botm[0] + 0.5 * delz
    h_lower = hp_lower + botm[-1] + 0.5 * delz
    strt = np.zeros((nlay, nrow, ncol))
    hstart = [botm[i] + 0.5 * delz + hp_lower for i in range(nlay)]
    strt[:, 0, 0] = hstart[:]

    nouter, ninner = 200, 300
    hclose, rclose, relax = 1e-6, 1e-6, 0.97

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
    nwtopts = ""
    if newton[idx]:
        nwtopts = "newton"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
        save_flows=True,
        newtonoptions=nwtopts,
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="None",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    if not newton[idx]:
        imsgwf.under_relaxation = ("DBD",)
        imsgwf.under_relaxation_theta = (0.9,)
        imsgwf.under_relaxation_kappa = (0.0001,)
        imsgwf.under_relaxation_gamma = (0.0,)

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
        gwf, save_specific_discharge=True, icelltype=laytyp, k=hk, k33=hk
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        ss=0.0001,
        sy=0.3,
        iconvert=0,
        transient=not steadystate[idx],
        steady_state=steadystate[idx],
    )

    # unsaturated zone Richards flow
    soil_data = get_uzr_soil_data("Celia1990-eq10-Haverkamp")
    uzr = flopy.mf6.ModflowGwfuzr(
        gwf,
        iunsat=1,
        phead_filerecord=f"{gwfname}.phd",
        storage_scheme="chord-slope",
        kr_averaging="harmonic",
        soil_model="Haverkamp",
        porosity=soil_data["porosity"],
        satres=soil_data["satres"],
        alphahvk=soil_data["alpha"],
        nhvk=soil_data["n"],
        betahvk=soil_data["beta"],
        khvk=soil_data["k"],
    )

    # constant head
    c = {0: [[(0, 0, 0), h_upper], [(nlay - 1, 0, 0), h_lower]]}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        print_flows=True,
        pname="CHD-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


saved_profiles = {}


def check_output(idx, test):
    model_name = "gwf_" + test.name
    fpth = os.path.join(test.workspace, f"{model_name}.phd")
    pfile = flopy.utils.HeadFile(fpth)
    times = pfile.get_times()

    saved_profiles[idx] = pfile.get_data(totim=times[-1])[:, 0, 0]

    fpth = os.path.join(test.workspace, f"{model_name}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    mg = grb.modelgrid
    nlay = mg.nlay
    dz = mg.delz.flatten()
    botm = mg.botm.flatten()
    depth = [-botm[ilay] - 0.5 * dz[ilay] for ilay in range(nlay)]

    figpth = os.path.join(test.workspace, f"pressure_head-{cases[idx]}.png")
    if PLOT_UZR_TESTS and not is_in_ci():
        plt.figure(figsize=(12, 6))

        for i, t in enumerate(times):
            pheads = pfile.get_data(totim=t)[:, 0, 0]
            label = f"{cases[idx]}: @t={times[-1]}"
            plt.plot(depth, pheads, label=label)

        # if last, plot all final profiles of previous in same graph
        if idx == len(cases) - 1:
            mrk = ["o", "x", "v", ">", "<", "^"]
            for iprev in range(len(cases) - 1):
                plt.plot(
                    depth,
                    saved_profiles[iprev],
                    marker=mrk[iprev],
                    linestyle="None",
                    label=f"{cases[iprev]}: @t={times[-1]}",
                )

        plt.xlim(0.0, 40.0)
        plt.ylim(-70.0, 0.0)
        plt.legend()
        plt.savefig(figpth)


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
