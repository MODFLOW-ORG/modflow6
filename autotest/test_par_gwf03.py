"""
Scaling parallel MODFLOW running a simple
(multi-)model setup on different partitionings
with constant head set at the lower-left corner.

a: 1 cpus, 1 model
b: 1 cpus, 4 models
c: 2 cpus, 4 models
d: 4 cpus, 4 models

The test is that for all configurations, the head
converges globally to the specified boundary value.
In general, the test can be used to compare parallel
vs. serial behavior on an identical problem.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["par_gwf03-a", "par_gwf03-b", "par_gwf03-c", "par_gwf03-d"]
ncpus = [1, 1, 2, 4]
domain_grid = [(1, 1), (2, 2), (2, 2), (2, 2)]
dis_shape = [(2, 100, 100), (2, 50, 50), (2, 50, 50), (2, 50, 50)]
delr = 100.0
delc = 100.0
head_initial = -1.0
cst_head_south_west = 435.0
hclose = 1.0e-8


def get_model_name(ix, iy):
    return f"model-{ix}-{iy}"


def get_simulation(idx, ws):
    name = cases[idx]
    nr_models_x = domain_grid[idx][0]
    nr_models_y = domain_grid[idx][1]

    nlay = dis_shape[idx][0]
    nrow = dis_shape[idx][1]
    ncol = dis_shape[idx][2]

    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for _ in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    rclose, relax = 1e-3, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        relaxation_factor=0.0,  # turn this off for comparison
    )

    # create models (and exchanges)
    for ix in range(nr_models_x):
        for iy in range(nr_models_y):
            add_model(sim, ix, iy, nr_models_x, nr_models_y, nlay, nrow, ncol)

    # add exchanges from west to east
    for iy in range(nr_models_y):
        for ix in range(nr_models_x - 1):
            name_west = get_model_name(ix, iy)
            name_east = get_model_name(ix + 1, iy)
            add_exchange_west_east(sim, name_west, name_east, nlay, nrow, ncol)

    # add exchange from south to north
    for ix in range(nr_models_x):
        for iy in range(nr_models_y - 1):
            name_south = get_model_name(ix, iy)
            name_north = get_model_name(ix, iy + 1)
            add_exchange_south_north(sim, name_south, name_north, nlay, nrow, ncol)

    return sim


def add_model(sim, ix, iy, nr_models_x, nr_models_y, nlay, nrow, ncol):
    # model spatial discretization
    shift_x = ix * ncol * delr
    shift_y = iy * nrow * delc
    model_name = get_model_name(ix, iy)

    # top/bot of the aquifer
    tops = [-100.0 * i for i in range(nlay + 1)]

    # hydraulic conductivity
    k11 = 10.0

    # initial head
    h_start = head_initial

    gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1 : nlay + 1],
        xorigin=shift_x,
        yorigin=shift_y,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{model_name}.hds",
        budget_filerecord=f"{model_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    if ix == 0 and iy == 0:
        # add SW corner BC
        sw_chd = [[(0, nrow - 1, 0), cst_head_south_west]]
        chd_spd_sw = {0: sw_chd}
        chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_sw)


def add_exchange_west_east(sim, name_west, name_east, nlay, nrow, ncol):
    exg_filename = f"we_{name_west}_{name_east}.gwfgwf"
    # exchangedata
    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (ilay, irow, ncol - 1),
            (ilay, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=name_west,
        exgmnameb=name_east,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename=exg_filename,
    )


def add_exchange_south_north(sim, name_south, name_north, nlay, nrow, ncol):
    exg_filename = f"sn_{name_south}_{name_north}.gwfgwf"

    # exchangedata
    angldegx = 90.0
    cdist = delc
    gwfgwf_data = [
        [
            (ilay, 0, icol),
            (ilay, nrow - 1, icol),
            1,
            delc / 2.0,
            delc / 2.0,
            delr,
            angldegx,
            cdist,
        ]
        for icol in range(ncol)
        for ilay in range(nlay)
    ]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=name_south,
        exgmnameb=name_north,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename=exg_filename,
    )


def build_models(idx, test):
    sim = get_simulation(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    mf6_sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    for mname in mf6_sim.model_names:
        m = mf6_sim.get_model(mname)
        hds = m.output.head().get_data().flatten()
        hds_compare = cst_head_south_west * np.ones_like(hds)
        assert np.allclose(hds, hds_compare, rtol=1.0e-6, atol=0.0001)


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=ncpus[idx],
    )
    test.run()
