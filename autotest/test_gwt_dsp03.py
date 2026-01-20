"""
Test the dispersion schemes in the gwt dispersion package for a three-dimensional
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

cases = ["dsp03a", "dsp03b"]
xt3d = [False, True]


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

    ihalfrow = int(nrow / 2)
    itri = np.zeros((nrow, ncol), dtype=int)
    itri[:ihalfrow, 1 : ncol - 1] = 1
    itri[ihalfrow:, 1 : ncol - 1] = 2
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

    cnclist = []
    chdlist = []
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                if j == ncol - 1:
                    # right
                    icellnum = itricellnum[i, j]
                    chdlist.append([(k, icellnum), 0.0])
                if j == 0:
                    # left
                    icellnum = itricellnum[i, j]
                    chdlist.append([(k, icellnum), 0.0])
                    cnclist.append([(k, icellnum), 1.0])

    cncs = {0: cnclist}
    chds = {0: chdlist}

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
        maxbound=len(chds),
        stress_period_data=chds,
        save_flows=False,
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
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwtname}.adv")

    # dispersion
    xt3d_off = not xt3d[idx]
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=xt3d_off,
        diffc=100.0,
        alh=0.0,
        alv=0.0,
        ath1=0.0,
        atv=0.0,
        filename=f"{gwtname}.dsp",
    )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # constant concentration
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt, stress_period_data=cncs, save_flows=False, pname="CNC-1"
    )

    # sources
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[[]], filename=f"{gwtname}.ssm")

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
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
        times = cobj.get_times()
        tdistplot = times[int(len(times) / 5)]
        conc = cobj.get_data(totim=tdistplot)
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # the time equal to 1/5 of perlen.
    cres1 = [
        1.0,
        0.94926607,
        0.9290872,
        0.88120143,
        0.85752405,
        0.80606293,
        0.78075725,
        0.72828023,
        0.70262037,
        0.65094507,
        0.62582375,
        0.5762093,
        0.5522329,
        0.50552886,
        0.48309497,
        0.43985539,
        0.41921253,
        0.37977384,
        0.36106262,
        0.32559608,
        0.3088772,
        0.27742691,
        0.26270138,
        0.23521484,
        0.22244005,
        0.19879526,
        0.1878988,
        0.16792926,
        0.1588215,
        0.14234021,
        0.13492598,
        0.12174938,
        0.11593943,
        0.10591361,
        0.10163723,
        0.09467278,
        0.09189463,
        0.08706135,
    ]
    cres1 = np.array(cres1)

    cres2 = [
        1.0,
        0.92688534,
        0.89781328,
        0.84000798,
        0.81152174,
        0.75519378,
        0.72764958,
        0.67348252,
        0.64719735,
        0.59578687,
        0.57102798,
        0.52286364,
        0.49984178,
        0.45529657,
        0.43416309,
        0.3934909,
        0.37433875,
        0.33767914,
        0.32054717,
        0.28793646,
        0.27281637,
        0.24420326,
        0.23104844,
        0.20631272,
        0.19504798,
        0.17402054,
        0.16455221,
        0.147035,
        0.1392606,
        0.12504531,
        0.11886233,
        0.10774676,
        0.10306018,
        0.09486133,
        0.0915894,
        0.0861513,
        0.08422876,
        0.08101295,
    ]
    cres2 = np.array(cres2)

    # Compare the first row in the layer with the answer and compare the
    # last row in the bottom layer with the answer.  This will verify that
    # the results are one-dimensional even though the model is three
    # dimensional
    creslist = [cres1, cres2]
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


@pytest.mark.slow
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
