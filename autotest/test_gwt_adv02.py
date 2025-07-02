"""
Test the advection schemes in the gwt advection package for a one-dimensional
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

cases = ["adv02a", "adv02b", "adv02c", "adv02d"]
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
    nlay, nrow, ncol = 1, 1, 100
    nper = 1
    perlen = [5.0]
    nstp = [200]
    tsmult = [1.0]
    steady = [True]
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

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
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    itri = np.zeros((nrow, ncol), dtype=int)
    itri[:, 1 : ncol - 1] = 1
    verts, iverts = grid_triangulator(itri, delr, delc)
    vertices, cell2d = cvfd_to_cell2d(verts, iverts)
    ncpl = len(cell2d)
    nvert = len(verts)

    c = {0: [[(0, ncpl - 1), 0.0000000]]}
    w = {0: [[(0, 0), 1.0, 1.0]]}

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
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
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
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # time step 200.
    cres1 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999999e-01,
                9.99999998e-01,
                9.99999996e-01,
                9.99999992e-01,
                9.99999986e-01,
                9.99999974e-01,
                9.99999953e-01,
                9.99999916e-01,
                9.99999853e-01,
                9.99999744e-01,
                9.99999560e-01,
                9.99999254e-01,
                9.99998750e-01,
                9.99997930e-01,
                9.99996616e-01,
                9.99994534e-01,
                9.99991280e-01,
                9.99986256e-01,
                9.99978596e-01,
                9.99967064e-01,
                9.99949912e-01,
                9.99924710e-01,
                9.99888122e-01,
                9.99835630e-01,
                9.99761195e-01,
                9.99656858e-01,
                9.99512258e-01,
                9.99314094e-01,
                9.99045508e-01,
                9.98685418e-01,
                9.98207806e-01,
                9.97580979e-01,
                9.96766847e-01,
                9.95720240e-01,
                9.94388311e-01,
                9.92710079e-01,
                9.90616155e-01,
                9.88028708e-01,
                9.84861727e-01,
                9.81021610e-01,
                9.76408132e-01,
                9.70915801e-01,
                9.64435624e-01,
                9.56857251e-01,
                9.48071482e-01,
                9.37973073e-01,
                9.26463768e-01,
                9.13455460e-01,
                8.98873378e-01,
                8.82659167e-01,
                8.64773747e-01,
                8.45199820e-01,
                8.23943904e-01,
                8.01037798e-01,
                7.76539388e-01,
                7.50532734e-01,
                7.23127419e-01,
                6.94457149e-01,
                6.64677639e-01,
                6.33963857e-01,
                6.02506698e-01,
                5.70509218e-01,
                5.38182550e-01,
                5.05741631e-01,
                4.73400912e-01,
                4.41370157e-01,
                4.09850495e-01,
                3.79030822e-01,
                3.49084661e-01,
                3.20167556e-01,
                2.92415050e-01,
                2.65941266e-01,
                2.40838116e-01,
                2.17175095e-01,
                1.94999626e-01,
                1.74337914e-01,
                1.55196221e-01,
                1.37562492e-01,
                1.21408260e-01,
                1.06690730e-01,
                9.33549731e-02,
                8.13361604e-02,
                7.05617550e-02,
                6.09536208e-02,
                5.24299924e-02,
                4.49072719e-02,
                3.83016284e-02,
                3.25303834e-02,
                2.75131759e-02,
                2.31729100e-02,
                1.94364916e-02,
                1.62353686e-02,
                1.35058926e-02,
                1.11895205e-02,
                9.23287950e-03,
                7.58771614e-03,
                6.21075097e-03,
                5.06346009e-03,
                4.11180172e-03,
                3.32590494e-03,
                2.67973538e-03,
                2.15075039e-03,
                1.71955399e-03,
                1.36956006e-03,
                1.08666981e-03,
                8.58968428e-04,
                6.76443820e-04,
                5.30729288e-04,
                4.14870946e-04,
                3.23119730e-04,
                2.50747289e-04,
                1.93884515e-04,
                1.49381177e-04,
                1.14684894e-04,
                8.77376133e-05,
                6.68877051e-05,
                5.08158521e-05,
                3.84730011e-05,
                2.90287464e-05,
                2.18286673e-05,
                1.63592796e-05,
                1.22194132e-05,
                9.09697239e-06,
                6.75017149e-06,
                4.99246589e-06,
                3.68051556e-06,
                2.70462002e-06,
                1.98115573e-06,
                1.44662627e-06,
                1.05300391e-06,
                7.64099584e-07,
                5.52747307e-07,
                3.98630253e-07,
                2.86609767e-07,
                2.05446572e-07,
                1.46826341e-07,
                1.04620334e-07,
                7.43267373e-08,
                5.26502675e-08,
                3.71870946e-08,
                2.61896435e-08,
                1.83917063e-08,
                1.28789035e-08,
                8.99310023e-09,
                6.26214049e-09,
                4.34838583e-09,
                3.01116258e-09,
                2.07945770e-09,
                1.43213687e-09,
                9.83662737e-10,
                6.73819846e-10,
                4.60347508e-10,
                3.13675549e-10,
                2.13175248e-10,
                1.44498191e-10,
                9.76935212e-11,
                6.58802074e-11,
                4.43137016e-11,
                2.97319183e-11,
                1.98983717e-11,
                1.32840240e-11,
                8.84640833e-12,
                4.39186836e-12,
            ]
        ]
    ]
    cres1 = np.array(cres1)

    cres2 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999999e-01,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999998e-01,
                9.99999993e-01,
                9.99999994e-01,
                1.00000000e00,
                1.00000002e00,
                1.00000002e00,
                9.99999999e-01,
                9.99999962e-01,
                9.99999940e-01,
                9.99999973e-01,
                1.00000007e00,
                1.00000018e00,
                1.00000017e00,
                9.99999959e-01,
                9.99999556e-01,
                9.99999217e-01,
                9.99999365e-01,
                1.00000035e00,
                1.00000203e00,
                1.00000340e00,
                1.00000262e00,
                9.99997675e-01,
                9.99987864e-01,
                9.99975605e-01,
                9.99967699e-01,
                9.99974537e-01,
                1.00000558e00,
                1.00005977e00,
                1.00011022e00,
                1.00008407e00,
                9.99839670e-01,
                9.99144904e-01,
                9.97660969e-01,
                9.94936340e-01,
                9.90414434e-01,
                9.83456888e-01,
                9.73381978e-01,
                9.59515273e-01,
                9.41247475e-01,
                9.18092979e-01,
                8.89742279e-01,
                8.56101925e-01,
                8.17317280e-01,
                7.73775496e-01,
                7.26088513e-01,
                6.75058258e-01,
                6.21628023e-01,
                5.66825293e-01,
                5.11701685e-01,
                4.57275459e-01,
                4.04481176e-01,
                3.54129825e-01,
                3.06881347e-01,
                2.63230002e-01,
                2.23501864e-01,
                1.87862728e-01,
                1.56334150e-01,
                1.28815062e-01,
                1.05106442e-01,
                8.49367728e-02,
                6.79864524e-02,
                5.39097837e-02,
                4.23536725e-02,
                3.29726109e-02,
                2.54398934e-02,
                1.94552984e-02,
                1.47496575e-02,
                1.10868391e-02,
                8.26370985e-03,
                6.10861676e-03,
                4.47887907e-03,
                3.25770405e-03,
                2.35085553e-03,
                1.68332282e-03,
                1.19616223e-03,
                8.43620069e-04,
                5.90594884e-04,
                4.10458367e-04,
                2.83227080e-04,
                1.94059595e-04,
                1.32043648e-04,
                8.92335599e-05,
                5.98980167e-05,
                3.99405773e-05,
                2.64592197e-05,
                1.74157687e-05,
                1.13907502e-05,
                7.40365644e-06,
                4.78259593e-06,
                3.07073947e-06,
                1.95984426e-06,
                1.24347214e-06,
                7.84372275e-07,
                4.91943284e-07,
                3.06795247e-07,
                1.90263904e-07,
                1.17346736e-07,
                7.19820791e-08,
                4.39185249e-08,
                2.66545669e-08,
                1.60925890e-08,
                9.66584556e-09,
                5.77619748e-09,
                3.43447585e-09,
                2.03199468e-09,
                1.19634152e-09,
                7.00946025e-10,
                4.08730000e-10,
                2.37211756e-10,
                1.37027700e-10,
                7.87911071e-11,
                4.50989518e-11,
                2.56980023e-11,
                1.45780272e-11,
                8.23354736e-12,
                4.63005332e-12,
                2.59249113e-12,
                1.44544628e-12,
                8.02528982e-13,
                4.43725011e-13,
                2.44332619e-13,
                1.33993251e-13,
                7.31876077e-14,
                3.98166029e-14,
                2.15765206e-14,
                1.16468205e-14,
                6.26267665e-15,
                3.35472405e-15,
                1.79025576e-15,
                9.51813339e-16,
                5.04177652e-16,
                2.66088975e-16,
                1.39925789e-16,
                7.33182482e-17,
                3.82811880e-17,
                1.99174186e-17,
                1.03269036e-17,
                5.33594407e-18,
                2.74771548e-18,
                1.41014305e-18,
                7.21291388e-19,
                3.67691580e-19,
                1.86885340e-19,
                9.45533475e-20,
                4.79404225e-20,
                2.37115056e-20,
                1.27310694e-20,
                2.67369800e-21,
            ]
        ]
    ]
    cres2 = np.array(cres2)

    cres3 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999999e-01,
                9.99999998e-01,
                9.99999997e-01,
                9.99999996e-01,
                9.99999992e-01,
                9.99999988e-01,
                9.99999977e-01,
                9.99999965e-01,
                9.99999935e-01,
                9.99999903e-01,
                9.99999820e-01,
                9.99999731e-01,
                9.99999499e-01,
                9.99999255e-01,
                9.99998614e-01,
                9.99997941e-01,
                9.99996190e-01,
                9.99994359e-01,
                9.99989643e-01,
                9.99984740e-01,
                9.99972298e-01,
                9.99959473e-01,
                9.99927548e-01,
                9.99895008e-01,
                9.99815884e-01,
                9.99736350e-01,
                9.99548166e-01,
                9.99362087e-01,
                9.98935116e-01,
                9.98520712e-01,
                9.97601139e-01,
                9.96726795e-01,
                9.94854670e-01,
                9.93113734e-01,
                9.89523166e-01,
                9.86262124e-01,
                9.79792542e-01,
                9.74060601e-01,
                9.63133239e-01,
                9.53698327e-01,
                9.36427537e-01,
                9.21907768e-01,
                8.96401341e-01,
                8.75537178e-01,
                8.40378779e-01,
                8.12414499e-01,
                7.67224282e-01,
                7.32294666e-01,
                6.78174275e-01,
                6.37542814e-01,
                5.77189801e-01,
                5.33200175e-01,
                4.70564614e-01,
                4.26261319e-01,
                3.65793900e-01,
                3.24305464e-01,
                2.70028089e-01,
                2.33916037e-01,
                1.88631553e-01,
                1.59426651e-01,
                1.24322194e-01,
                1.02384383e-01,
                7.71079810e-02,
                6.18067033e-02,
                4.49072754e-02,
                3.50006334e-02,
                2.45124747e-02,
                1.85605248e-02,
                1.25201360e-02,
                9.20274911e-03,
                5.97547120e-03,
                4.26070056e-03,
                2.66158628e-03,
                1.83980553e-03,
                1.10519401e-03,
                7.40179413e-04,
                4.27404601e-04,
                2.77183690e-04,
                1.53799502e-04,
                9.65363514e-05,
                5.14550934e-05,
                3.12434950e-05,
                1.59926916e-05,
                9.38944096e-06,
                4.61422039e-06,
                2.61810100e-06,
                1.23478054e-06,
                6.76699550e-07,
                3.06141480e-07,
                1.61921388e-07,
                7.02044355e-08,
                3.57908578e-08,
                1.48458395e-08,
                7.27822667e-09,
                2.87772343e-09,
                1.35034049e-09,
                5.04924114e-10,
                2.24462051e-10,
                7.79256525e-11,
                3.19978311e-11,
                9.80115906e-12,
                3.42172760e-12,
                7.38995558e-13,
                1.04998282e-13,
                -6.08063778e-14,
                -7.00202752e-14,
                -4.57797962e-14,
                -2.73273037e-14,
                -1.26821070e-14,
                -5.66553096e-15,
                -1.79362990e-15,
                -4.40132769e-16,
                1.05900605e-16,
                2.30438418e-16,
                2.68347207e-16,
                2.60728717e-16,
                2.32632994e-16,
                1.97839789e-16,
                1.63276013e-16,
                1.32036738e-16,
                1.05203244e-16,
                8.28614219e-17,
                6.46461341e-17,
                5.00217396e-17,
                3.84212334e-17,
                2.93113868e-17,
                2.22197912e-17,
                1.67426198e-17,
                1.25429101e-17,
                9.34448747e-18,
                6.92421024e-18,
                5.10394011e-18,
                3.74295077e-18,
                2.73111210e-18,
                1.98297560e-18,
                1.43276787e-18,
                1.03023499e-18,
                7.37248326e-19,
                5.25070294e-19,
                3.72178066e-19,
                2.62549962e-19,
                1.84329161e-19,
                1.28790859e-19,
                8.95504779e-20,
                6.19612940e-20,
                4.26594479e-20,
                2.92226046e-20,
                1.99155401e-20,
                1.35016962e-20,
                9.10444343e-21,
                4.57287749e-21,
            ]
        ]
    ]
    cres3 = np.array(cres3)

    cres4 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999999e-01,
                9.99999998e-01,
                9.99999997e-01,
                9.99999996e-01,
                9.99999993e-01,
                9.99999988e-01,
                9.99999980e-01,
                9.99999967e-01,
                9.99999945e-01,
                9.99999909e-01,
                9.99999849e-01,
                9.99999748e-01,
                9.99999580e-01,
                9.99999300e-01,
                9.99998835e-01,
                9.99998064e-01,
                9.99996788e-01,
                9.99994679e-01,
                9.99991207e-01,
                9.99985519e-01,
                9.99976245e-01,
                9.99961221e-01,
                9.99937056e-01,
                9.99898507e-01,
                9.99837576e-01,
                9.99742235e-01,
                9.99594683e-01,
                9.99369014e-01,
                9.99028179e-01,
                9.98520185e-01,
                9.97773486e-01,
                9.96691640e-01,
                9.95147449e-01,
                9.92976938e-01,
                9.89973781e-01,
                9.85884941e-01,
                9.80408526e-01,
                9.73194950e-01,
                9.63852496e-01,
                9.51958228e-01,
                9.37074848e-01,
                9.18773513e-01,
                8.96661955e-01,
                8.70416363e-01,
                8.39814727e-01,
                8.04768610e-01,
                7.65349931e-01,
                7.21809331e-01,
                6.74583139e-01,
                6.24286813e-01,
                5.71694587e-01,
                5.17708928e-01,
                4.63475646e-01,
                4.10621606e-01,
                3.60090882e-01,
                3.12561498e-01,
                2.68548170e-01,
                2.28398779e-01,
                1.92300754e-01,
                1.60295787e-01,
                1.32300488e-01,
                1.08130463e-01,
                8.75255237e-02,
                7.01741203e-02,
                5.57355716e-02,
                4.38591512e-02,
                3.41995388e-02,
                2.64285307e-02,
                2.02431998e-02,
                1.53708955e-02,
                1.15715937e-02,
                8.63815025e-03,
                6.39500342e-03,
                4.69581865e-03,
                3.42049930e-03,
                2.47190356e-03,
                1.77252555e-03,
                1.26132360e-03,
                8.90813017e-04,
                6.24488550e-04,
                4.34601338e-04,
                3.00286592e-04,
                2.06019188e-04,
                1.40363342e-04,
                9.49774268e-05,
                6.38341487e-05,
                4.26182859e-05,
                2.82678226e-05,
                1.86287813e-05,
                1.21986984e-05,
                7.93815930e-06,
                5.13384506e-06,
                3.30005212e-06,
                2.10858694e-06,
                1.33934239e-06,
                8.45779171e-07,
                5.31034450e-07,
                3.31530143e-07,
                2.05821848e-07,
                1.27075205e-07,
                7.80302654e-08,
                4.76573231e-08,
                2.89528335e-08,
                1.74975602e-08,
                1.05200557e-08,
                6.29276848e-09,
                3.74521020e-09,
                2.21793995e-09,
                1.30703695e-09,
                7.66512576e-10,
                4.47368886e-10,
                2.59871508e-10,
                1.50249706e-10,
                8.64692928e-11,
                4.95369461e-11,
                2.82485119e-11,
                1.60400425e-11,
                9.06356925e-12,
                5.10160108e-12,
                2.85609391e-12,
                1.59463729e-12,
                8.83001720e-13,
                4.88912295e-13,
                2.66402206e-13,
                1.47206301e-13,
                7.85891341e-14,
                4.65296803e-14,
                3.41047008e-14,
                2.48440107e-14,
                1.80183366e-14,
                1.30166484e-14,
                9.36861656e-15,
                6.71933416e-15,
                4.80311981e-15,
                3.42235516e-15,
                2.43094613e-15,
                1.72149946e-15,
                1.21547008e-15,
                8.55670685e-16,
                6.00634004e-16,
                4.20405055e-16,
                2.93422427e-16,
                2.04219766e-16,
                1.41739986e-16,
                9.81042000e-17,
                6.77163197e-17,
                4.66142906e-17,
                3.20017576e-17,
                2.19112061e-17,
                1.49625226e-17,
                7.71070649e-18,
            ]
        ]
    ]
    cres4 = np.array(cres4)

    creslist = [cres1, cres2, cres3, cres4]

    assert np.allclose(creslist[idx], conc), (
        "simulated concentrations do not match with known solution."
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
