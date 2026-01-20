"""
Test the seepage face boundary for the Richards based unsaturated
zone package UZR. A transient 2D (xz) model with unequal constant heads
on either side.
"""

import os

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pytest
from framework import TestFramework
from gwf_test_utils import PLOT_UZR_TESTS, get_uzr_soil_data
from modflow_devtools.misc import is_in_ci

cases = ["npf", "npf-drn", "uzr", "uzr-drn", "uzr-spf"]
use_uzr = [False, False, True, True, True]
use_seepage = [False, False, False, False, True]
use_drain = [False, True, False, True, False]

refinement = 1
height = 200.0
width = 800.0
nlay, nrow, ncol = refinement * 40, 1, refinement * 40

delr = width / ncol  # cm
delc = delr
delz = height / nlay

water_tables = {}


def build_models(idx, test):
    nper = 1
    nstp = 48
    dtmin = 60
    dt = 3600
    perlen = nstp * dt  # s
    tsmult = 1.0
    tdis_rc = [(perlen, nstp, tsmult)]

    top = height
    icelltype = 1
    iconvert = 1
    newtonopts = "newton"

    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
    hk = 0.01  # cm/s

    # saturated lowest cell:
    h_left = 180.0  # cm
    h_right = 20.0  # cm

    strt = np.zeros((nlay, nrow, ncol)) + h_right

    nouter, ninner = 300, 300
    hclose, rclose, relax = 1e-4, 1e-4, 1.0

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    ats_filerecord = name + ".ats"
    atsperiod = [(0, dtmin, dtmin, dt, 2.0, 2.0)]
    tdis.ats.initialize(
        maxats=len(atsperiod),
        perioddata=atsperiod,
        filename=ats_filerecord,
    )

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        save_flows=True,
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
        newtonoptions=newtonopts,
    )

    # create iterative model solution and register the gwf model with it
    # imsgwf = flopy.mf6.ModflowIms(
    #     sim,
    #     print_option="SUMMARY",
    #     outer_dvclose=hclose,
    #     outer_maximum=nouter,
    #     inner_maximum=ninner,
    #     inner_dvclose=hclose,
    #     rcloserecord=rclose,
    #     linear_acceleration="BICGSTAB",
    #     scaling_method="NONE",
    #     reordering_method="NONE",
    #     relaxation_factor=relax,
    #     filename=f"{gwfname}.ims",
    # )
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        inner_dvclose=hclose,
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
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
        save_specific_discharge=True,
        icelltype=icelltype,
        k=hk,
        k33=hk,
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf, ss=0.0001, sy=0.3, iconvert=iconvert, transient={0: True}
    )

    # unsaturated zone Richards flow
    if use_uzr[idx]:
        soil_data = get_uzr_soil_data("Celia1990-eq10-Haverkamp")
        uzr = flopy.mf6.ModflowGwfuzr(
            gwf,
            iunsat=1,
            storage_scheme="chord-slope",
            kr_averaging="arithmetic",
            soil_model="Haverkamp",
            porosity=soil_data["porosity"],
            satres=soil_data["satres"],
            alphahvk=soil_data["alpha"],
            nhvk=soil_data["n"],
            betahvk=soil_data["beta"],
            khvk=soil_data["k"],
            phead_filerecord=f"{gwfname}.phd",
        )

    # boundary data
    chd_data = []
    spf_data = []
    drn_data = []

    # left boundary
    for ilay in range(nlay):
        z = height - 0.5 * delz - ilay * delz
        if z < h_left:
            chd_data.append([(ilay, 0, 0), h_left])
    # right boundary
    for ilay in range(nlay):
        z = height - 0.5 * delz - ilay * delz
        if z < h_right:
            chd_data.append([(ilay, 0, ncol - 1), h_right])
        else:
            spf_data.append([(ilay, 0, ncol - 1), 0.5 * delr, delc * delz])

            dcond = hk * delc * delz / (0.5 * delr)
            drn_data.append([(ilay, 0, ncol - 1), z, dcond])

    chd_spd = {0: chd_data}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chd_spd),
        stress_period_data=chd_spd,
        save_flows=True,
        print_flows=True,
        pname="CHD-1",
    )

    if use_seepage[idx]:
        spf_spd = {0: spf_data}
        spf = flopy.mf6.ModflowGwfspf(
            gwf,
            maxbound=len(spf_spd),
            stress_period_data=spf_spd,
            save_flows=True,
            print_flows=True,
            pname="SPF-1",
        )

    if use_drain[idx]:
        drn_spd = {0: drn_data}
        drn = flopy.mf6.ModflowGwfdrn(
            gwf,
            stress_period_data=drn_spd,
            print_input=True,
        )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    model_name = "gwf_" + test.name

    # plot heads and wt
    fpth = os.path.join(test.workspace, f"{model_name}.hds")
    hfile = flopy.utils.HeadFile(fpth)
    heads = hfile.get_data(idx=-1)
    wt = flopy.utils.postprocessing.get_water_table(heads)
    water_tables[cases[idx]] = wt

    heads[heads == -1e30] = np.nan
    x = [icol * width / ncol + delr for icol in range(ncol)]

    m = [
        "o",
        "x",
        "v",
        ">",
        "<",
        "^",
    ]
    if idx == len(cases) - 1:
        plt.figure()
        for idx, case in enumerate(cases):
            wtable = water_tables[case]
            plt.plot(x, wtable, marker=m[idx], label=case)
        figpth = os.path.join(test.workspace, "wt-all.png")
        plt.xlim(0.0, width)
        plt.ylim(0.0, height)
        plt.legend()
        plt.savefig(figpth)

    if PLOT_UZR_TESTS and not is_in_ci():
        levels = np.arange(0.0, 200.0, 10.0)
        plt.figure()
        ax = plt.gca()
        ax.imshow(heads[:, 0, :], extent=[0.0, width, height, 0])
        figpth = os.path.join(test.workspace, f"head-{cases[idx]}.png")
        plt.savefig(figpth)

        plt.figure()
        plt.plot(x, wt)
        figpth = os.path.join(test.workspace, f"wt-{cases[idx]}.png")
        plt.ylim(0.0, 200.0)
        plt.savefig(figpth)

        if use_uzr[idx]:
            # plot pressure heads and wt
            fpth = os.path.join(test.workspace, f"{model_name}.phd")
            pfile = flopy.utils.HeadFile(fpth)
            pheads = pfile.get_data(idx=-1)

            levels = np.arange(-200.0, 200.0, 10.0)
            plt.figure()
            ax = plt.gca()
            cs = ax.contour(pheads[:, 0, :], levels, extent=[0.0, width, height, 0])
            _ = ax.clabel(cs)
            figpth = os.path.join(test.workspace, f"uzr-phead-{cases[idx]}.png")
            plt.savefig(figpth)

            level_idx = np.where(abs(cs.get_array()) < 1e-10)[0][0]
            x = cs.get_paths()[level_idx].vertices[:, 0]
            p = cs.get_paths()[level_idx].vertices[:, 1]
            plt.figure()
            plt.plot(x, p)
            plt.ylim(0.0, height)
            figpth = os.path.join(test.workspace, f"uzr-wt-{cases[idx]}.png")
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
