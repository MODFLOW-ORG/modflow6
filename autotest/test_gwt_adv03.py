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

cases = ["adv03a", "adv03b", "adv03c"]
scheme = ["upstream", "central", "tvd"]


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
    width = 5.0
    length = 20.0
    depth = 5.0
    nper = 1
    delr = length / ncol
    delc = width / nrow
    delz = depth / nlay
    top = 1.0
    botm = np.linspace(top - delz, top - nlay * delz, nlay)
    strt = 1.0
    hk = 1.0
    laytyp = 0
    porosity = 0.1
    velocity_x = 0.5
    qwell = velocity_x * porosity * delc * delz
    specific_discharge = 1.0 # concentration
    timetoend = length / velocity_x

    perlen = [timetoend]
    nstp = [100]
    tsmult = [1.0]

    nouter, ninner = 100, 300
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
        outer_maximum=nouter,
        inner_maximum=ninner,
        outer_dvclose=hclose,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        relaxation_factor=relax,
        linear_acceleration="BICGSTAB",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    itri = np.zeros((nrow, ncol), dtype=int)
    itri[:int(nrow/2), 1 : ncol - 1] = 1
    itri[int(nrow/2):, 1 : ncol - 1] = 2
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
                    wellist.append([(k, icellnum), qwell, specific_discharge])

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
        xt3doptions=True,
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
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 16, "GENERAL")],
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
        outer_maximum=nouter,
        inner_maximum=ninner,
        outer_dvclose=hclose,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
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
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

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
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 16, "GENERAL")],
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
        9.78263289e-01,
        9.64602117e-01,
        9.52100342e-01,
        9.05738549e-01,
        8.71038362e-01,
        7.60312327e-01,
        6.83635709e-01,
        4.97262462e-01,
        3.88629846e-01,
        2.10421365e-01,
        1.31466051e-01,
        4.84712709e-02,
        2.25083551e-02,
        5.29149297e-03,
        1.70775757e-03,
        2.35723950e-04,
        4.81264694e-05,
        3.32524217e-06,
        2.78104706e-07,
        -1.46890248e-08,
        -1.06009800e-08,
        -7.29893033e-09,
        -3.19570174e-09,
        -1.55626640e-09,
        -5.21166446e-10,
        -1.60129996e-10,
        -5.26729496e-11,
        -2.26445075e-11,
        -1.02289090e-11,
        -3.57047421e-12,
        -1.41048657e-12,
        -2.36086048e-13,
        4.40920030e-15,
        6.78407639e-14,
        3.77485549e-14,
        2.05570945e-14,
        1.09960343e-14,
        4.97448017e-15,
    ]
    cres3 = np.array(cres3)

    # Compare the first row in the layer with the answer and compare the
    # last row in the bottom layer with the answer.  This will verify that
    # the results are one-dimensional even though the model is three
    # dimensional
    creslist = [cres1, cres2, cres3]
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
