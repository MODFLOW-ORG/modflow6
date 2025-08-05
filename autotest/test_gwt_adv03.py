"""
Test the advection schemes in the gwt advection package for a three-dimensional
model grid of triangular cells.  The cells are created by starting with a
regular grid of squares and then cutting every cell into a triangle, except
the first and last.
"""

import os

import flopy
import flopy.utils.cvfdutil
import numpy as np
import pytest
from framework import TestFramework

cases = ["adv03a", "adv03b", "adv03c", "adv03d"]
scheme = ["upstream", "central", "tvd", "utvd"]


def grid_triangulator(itri, delr, delc):
    nrow, ncol = itri.shape
    if np.isscalar(delr):
        delr = delr * np.ones(ncol)
    if np.isscalar(delc):
        delc = delc * np.ones(nrow)
    regular_grid = flopy.discretization.StructuredGrid(delc, delr)
    vertdict = {}
    icell = 0
    for i in range(nrow):
        for j in range(ncol):
            vs = regular_grid.get_cell_vertices(i, j)
            if itri[i, j] == 0:
                vertdict[icell] = [vs[0], vs[1], vs[2], vs[3], vs[0]]
                icell += 1
            elif itri[i, j] == 1:
                vertdict[icell] = [vs[0], vs[1], vs[3], vs[0]]
                icell += 1
                vertdict[icell] = [vs[3], vs[1], vs[2], vs[3]]
                icell += 1
            elif itri[i, j] == 2:
                vertdict[icell] = [vs[0], vs[2], vs[3], vs[0]]
                icell += 1
                vertdict[icell] = [vs[0], vs[1], vs[2], vs[0]]
                icell += 1
            else:
                raise Exception(f"Unknown itri value: {itri[i, j]}")
    verts, iverts = flopy.utils.cvfdutil.to_cvfd(vertdict)
    return verts, iverts


def cvfd_to_cell2d(verts, iverts):
    vertices = []
    for i in range(verts.shape[0]):
        x = verts[i, 0]
        y = verts[i, 1]
        vertices.append([i, x, y])
    cell2d = []
    for icell2d, vs in enumerate(iverts):
        points = [tuple(verts[ip]) for ip in vs]
        xc, yc = flopy.utils.cvfdutil.centroid_of_polygon(points)
        cell2d.append([icell2d, xc, yc, len(vs), *vs])
    return vertices, cell2d


def build_models(idx, test):
    nlay, nrow, ncol = 5, 10, 20
    nper = 1
    delr = 1.0
    delc = 1.0
    delz = 1.0
    top = 1.0
    botm = np.linspace(top - delz, top - nlay * delz, nlay)
    strt = 1.0
    hk = 1.0
    laytyp = 0
    porosity = 0.1
    qwell = 1.0
    specific_discharge = qwell / delr / delz
    timetoend = float(ncol) * delc * porosity / specific_discharge

    perlen = [timetoend]
    nstp = [50]
    tsmult = [1.0]
    steady = [True]

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
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    itri = np.zeros((nrow, ncol), dtype=int)
    itri[int(nrow / 2) :, 1 : ncol - 1] = 1
    itri[: int(nrow / 2), 1 : ncol - 1] = 2
    verts, iverts = grid_triangulator(itri, delr, delc)
    vertices, cell2d = cvfd_to_cell2d(verts, iverts)
    ncpl = len(cell2d)
    nvert = len(verts)

    # A grid array that has the cellnumber of the first triangular cell in
    # the original grid
    itricellnum = np.empty((nrow, ncol), dtype=int)
    icell = 0
    for i in range(nrow):
        for j in range(ncol):
            itricellnum[i, j] = icell
            if itri[i, j] != 0:
                icell += 2
            else:
                icell += 1

    chdlist = []
    wellist = []
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                if j == ncol - 1:
                    icellnum = itricellnum[i, j]
                    chdlist.append([(k, icellnum), 0.0])
                if j == 0:
                    icellnum = itricellnum[i, j]
                    wellist.append([(k, icellnum), qwell, 1.0])

    c = {0: chdlist}
    w = {0: wellist}

    disv = flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=nlay,
        ncpl=ncpl,
        nvert=nvert,
        top=top,
        botm=botm,
        vertices=vertices,
        cell2d=cell2d,
        filename=f"{gwfname}.disv",
    )

    # dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol,
    #                              delr=delr, delc=delc,
    #                              top=top, botm=botm,
    #                              idomain=np.ones((nlay, nrow, ncol), dtype=int),
    #                              filename='{}.dis'.format(gwfname))

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=laytyp,
        xt3doptions=[()],
        k=hk,
        k33=hk,
        save_specific_discharge=True,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        pname="CHD-1",
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

    disv = flopy.mf6.ModflowGwtdisv(
        gwt,
        nlay=nlay,
        ncpl=ncpl,
        nvert=nvert,
        top=top,
        botm=botm,
        vertices=vertices,
        cell2d=cell2d,
        filename=f"{gwtname}.disv",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme[idx], filename=f"{gwtname}.adv")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

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
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        "conc_obs.csv": [
            ("(1-2)", "CONCENTRATION", (0, 1)),
            ("(1-50)", "CONCENTRATION", (0, 49)),
        ],
        "flow_obs.csv": [
            ("c10-c11", "FLOW-JA-FACE", (0, 9), (0, 10)),
            ("c50-c51", "FLOW-JA-FACE", (0, 49), (0, 50)),
            ("c99-c100", "FLOW-JA-FACE", (0, 98), (0, 99)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname="conc_obs",
        filename=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
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
        times = cobj.get_times()
        tdistplot = times[int(len(times) / 5)]
        conc = cobj.get_data(totim=tdistplot)
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # the time equal to 1/5 of perlen.
    cres1 = [
        9.75305991e-01,
        9.52167956e-01,
        9.13498872e-01,
        8.56445999e-01,
        7.81408232e-01,
        6.92105105e-01,
        5.94528731e-01,
        4.95411393e-01,
        4.00833986e-01,
        3.15336441e-01,
        2.41600024e-01,
        1.80579344e-01,
        1.31891325e-01,
        9.42846759e-02,
        6.60695310e-02,
        4.54473575e-02,
        3.07275409e-02,
        2.04445131e-02,
        1.34007670e-02,
        8.66203225e-03,
        5.52640971e-03,
        3.48306266e-03,
        2.17022758e-03,
        1.33775921e-03,
        8.16313473e-04,
        4.93397804e-04,
        2.95551417e-04,
        1.75541689e-04,
        1.03427957e-04,
        6.04769621e-05,
        3.51080337e-05,
        2.02416209e-05,
        1.15945184e-05,
        6.60031880e-06,
        3.73515492e-06,
        2.10185144e-06,
        1.17640119e-06,
        4.52726328e-07,
    ]
    cres1 = np.array(cres1)

    cres2 = [
        9.91666434e-01,
        9.86953589e-01,
        9.80778200e-01,
        9.56402197e-01,
        8.95573236e-01,
        7.94392952e-01,
        6.64374858e-01,
        5.24436564e-01,
        3.92222699e-01,
        2.79262303e-01,
        1.90214211e-01,
        1.24507826e-01,
        7.86404693e-02,
        4.81012219e-02,
        2.85822809e-02,
        1.65449352e-02,
        9.35203645e-03,
        5.17291859e-03,
        2.80515816e-03,
        1.49375662e-03,
        7.82221788e-04,
        4.03333373e-04,
        2.05011749e-04,
        1.02829466e-04,
        5.09423767e-05,
        2.49472500e-05,
        1.20857205e-05,
        5.79592865e-06,
        2.75322526e-06,
        1.29620095e-06,
        6.05117190e-07,
        2.80245975e-07,
        1.28825482e-07,
        5.87741819e-08,
        2.66893805e-08,
        1.19138575e-08,
        5.57865554e-09,
        1.01007992e-09,
    ]
    cres2 = np.array(cres2)

    cres3 = [
        9.75305991e-01,
        9.61567354e-01,
        9.44318192e-01,
        8.99780324e-01,
        8.53444404e-01,
        7.57522910e-01,
        6.71042660e-01,
        5.31510001e-01,
        4.26626889e-01,
        2.94209725e-01,
        2.11149438e-01,
        1.26488690e-01,
        8.13031158e-02,
        4.27834908e-02,
        2.48703736e-02,
        1.16850355e-02,
        6.22057950e-03,
        2.65348164e-03,
        1.30994628e-03,
        5.15075141e-04,
        2.38494633e-04,
        8.75806378e-05,
        3.84075388e-05,
        1.33165430e-05,
        5.57286464e-06,
        1.84101848e-06,
        7.40051254e-07,
        2.34149161e-07,
        9.11456220e-08,
        2.77537863e-08,
        1.06492160e-08,
        3.18846190e-09,
        1.21166402e-09,
        3.76715755e-10,
        1.31716572e-10,
        3.23063284e-11,
        9.77588885e-12,
        -1.56168535e-12,
    ]
    cres3 = np.array(cres3)

    cres4 = [
        9.75305992e-01,
        9.60923712e-01,
        9.40345041e-01,
        8.98518356e-01,
        8.43452227e-01,
        7.50966209e-01,
        6.55085404e-01,
        5.23794573e-01,
        4.15360236e-01,
        2.94697679e-01,
        2.10720136e-01,
        1.35962806e-01,
        8.93878174e-02,
        5.38698838e-02,
        3.32165717e-02,
        1.88821558e-02,
        1.10910530e-02,
        5.99079841e-03,
        3.39811390e-03,
        1.75392691e-03,
        9.71354219e-04,
        4.81791201e-04,
        2.62200432e-04,
        1.25921717e-04,
        6.72680080e-05,
        3.15465073e-05,
        1.64347475e-05,
        7.56625500e-06,
        3.83544217e-06,
        1.73771238e-06,
        8.59572442e-07,
        3.83960616e-07,
        1.85884337e-07,
        8.19811804e-08,
        3.89673178e-08,
        1.69141522e-08,
        8.17467605e-09,
        1.47129715e-09,
    ]
    cres4 = np.array(cres4)

    # Compare the first row in the layer with the answer and compare the
    # last row in the bottom layer with the answer.  This will verify that
    # the results are one-dimensional even though the model is three
    # dimensional
    creslist = [cres1, cres2, cres3, cres4]
    ncellsperrow = cres1.shape[0]
    assert np.allclose(creslist[idx], conc[0, 0, 0:ncellsperrow]), (
        "simulated concentrations do not match with known solution.",
        creslist[idx],
        conc[0, 0, -ncellsperrow:],
    )
    assert np.allclose(creslist[idx], conc[0, 0, -ncellsperrow:]), (
        "simulated concentrations do not match with known solution.",
        creslist[idx],
        conc[0, 0, -ncellsperrow:],
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
