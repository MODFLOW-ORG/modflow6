"""
Test the storage formulations for the Richards based unsaturated 
zone package UZR. Test flow in a single homogeneous column with
upper and lower Dirichlet boundary.
"""

import os
import matplotlib.pyplot as plt
import flopy
import numpy as np
import pytest
from modflow_devtools.misc import is_in_ci
from gwf_test_utils import get_uzr_soil_data, get_balance_error, PLOT_UZR_TESTS
from framework import TestFramework


cases = ["chord-slope", "mod-picard"]
scheme = ["chord-slope", "modified-picard"]


def build_models(idx, test):
    column_height = 42.0
    nlay, nrow, ncol = 42, 1, 1
    nper = 1
    perlen = [360.0]  # s
    nstp = [3]
    tsmult = [1.0]
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
    soil_data = get_uzr_soil_data("Celia1990-eq10-Haverkamp")
    uzr = flopy.mf6.ModflowGwfuzr(
        gwf,
        iunsat=1,
        storage_scheme=scheme[idx],
        kr_averaging="geometric",
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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):

    model_name = "gwf_" + test.name

    fpth = os.path.join(test.workspace, f"{model_name}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    mg = grb.modelgrid
    nlay = mg.nlay
    dz = mg.delz.flatten()
    botm = mg.botm.flatten()

    fpth = os.path.join(test.workspace, f"{model_name}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data(idx=-1).flatten()

    pheads = [heads[ilay] - botm[ilay] - 0.5 * dz[ilay] for ilay in range(nlay)]
    depth = [-botm[ilay] - 0.5 * dz[ilay] for ilay in range(nlay)]

    if PLOT_UZR_TESTS and not is_in_ci():
        plt.plot(depth, pheads)
        plt.xlim(0.0, 40.0)
        plt.ylim(-70.0, -10.0)
        plt.savefig(f"pressure_head-{cases[idx]}.png")

    list_pth = os.path.join(test.workspace, f"{model_name}.lst")
    error = get_balance_error(list_pth)
    assert (
        abs(error) < 0.00001
    ), f"Cumulative balance error = {error} for {model_name} too large"


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
