"""
Test for the Richards based unsaturated zone package UZR.
"""

import os

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pytest
from framework import TestFramework
from gwf_test_utils import PLOT_UZR_TESTS, get_balance_error, get_uzr_soil_data
from modflow_devtools.misc import is_in_ci

cases = ["celia2-dense", "celia2-dt144", "celia2-dt720", "celia2-dt3600"]
dt = [10.0, 144.0, 720.0, 3600.0]

gwfname = "celia"

def build_models(idx, test):
    column_height = 105.0
    nlay, nrow, ncol = 42, 1, 1
    nper = 1
    perlen = [36 * 3600]  # s
    nstp = [int(perlen[0] / dt[idx])]
    tsmult = [1.0]
    delr = 10.0  # cm
    delc = 10.0
    delz = column_height / nlay
    top = delz
    laytyp = 0
    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
    hk = 0.00922  # cm/s

    # saturated lowest cell:
    hp_upper = -75.0
    hp_lower = -1000.0
    h_upper = hp_upper + botm[0] + 0.5 * delz
    h_lower = hp_lower + botm[-1] + 0.5 * delz
    strt = np.zeros((nlay, nrow, ncol))
    hstart = [botm[i] + 0.5 * delz + hp_lower for i in range(nlay)]
    strt[:, 0, 0] = hstart[:]

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
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    gwf.name_file.save_flows = True

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
        gwf, save_specific_discharge=True, icelltype=laytyp, k=hk, k33=hk
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf, ss=0.0001, sy=0.3, iconvert=0, transient={0: True}
    )

    # unsaturated zone Richards flow
    soil_data = get_uzr_soil_data("Celia1990-eq13-VanGenuchten")
    uzr = flopy.mf6.ModflowGwfuzr(
        gwf,
        iunsat=1,
        storage_scheme="modified-picard",
        kr_averaging="arithmetic",
        soil_model="VanGenuchten",
        porosity=soil_data["porosity"],
        satres=soil_data["satres"],
        alphavgn=soil_data["alpha"],
        nvgn=soil_data["n"],
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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, f"{gwfname}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    mg = grb.modelgrid
    nlay = mg.nlay
    dz = mg.delz.flatten()
    botm = mg.botm.flatten()

    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data(idx=-1).flatten()

    pheads = [heads[ilay] - botm[ilay] - 0.5 * dz[ilay] for ilay in range(nlay)]
    depth = [-botm[ilay] - 0.5 * dz[ilay] for ilay in range(nlay)]

    if PLOT_UZR_TESTS and not is_in_ci():
        plt.figure()
        plt.plot(depth, pheads)
        plt.savefig(f"pressure_head-{cases[idx]}.png")

    list_pth = os.path.join(test.workspace, f"{gwfname}.lst")
    error = get_balance_error(list_pth)
    assert abs(error) < 0.00001, (
        f"Cumulative balance error = {error} for {gwfname} too large"
    )


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
