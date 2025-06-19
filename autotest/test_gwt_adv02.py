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

cases = ["adv02a", "adv02b", "adv02c"]
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
                9.99999998e-01,
                9.99999996e-01,
                9.99999990e-01,
                9.99999979e-01,
                9.99999954e-01,
                9.99999902e-01,
                9.99999795e-01,
                9.99999580e-01,
                9.99999152e-01,
                9.99998319e-01,
                9.99996731e-01,
                9.99993757e-01,
                9.99988295e-01,
                9.99978459e-01,
                9.99961085e-01,
                9.99930988e-01,
                9.99879863e-01,
                9.99794694e-01,
                9.99655559e-01,
                9.99432657e-01,
                9.99082457e-01,
                9.98542886e-01,
                9.97727574e-01,
                9.96519340e-01,
                9.94763270e-01,
                9.92259990e-01,
                9.88760005e-01,
                9.83960181e-01,
                9.77503642e-01,
                9.68984370e-01,
                9.57957654e-01,
                9.43957177e-01,
                9.26518878e-01,
                9.05210968e-01,
                8.79668487e-01,
                8.49629868e-01,
                8.14972180e-01,
                7.75741223e-01,
                7.32172621e-01,
                6.84700553e-01,
                6.33951743e-01,
                5.80724281e-01,
                5.25954375e-01,
                4.70801201e-01,
                4.16955335e-01,
                3.65457844e-01,
                3.17016523e-01,
                2.72168709e-01,
                2.31276918e-01,
                1.94535762e-01,
                1.61987773e-01,
                1.33545414e-01,
                1.09016570e-01,
                8.81310767e-02,
                7.05663170e-02,
                5.59703742e-02,
                4.39818117e-02,
                3.42455867e-02,
                2.64250261e-02,
                2.02100890e-02,
                1.53223494e-02,
                1.15172457e-02,
                8.58418638e-03,
                6.34508032e-03,
                4.65180911e-03,
                3.38307547e-03,
                2.44097709e-03,
                1.74756758e-03,
                1.24158770e-03,
                8.75482929e-04,
                6.12769793e-04,
                4.25772639e-04,
                2.93723748e-04,
                2.01201076e-04,
                1.36867183e-04,
                9.24683038e-05,
                6.20521199e-05,
                4.13651949e-05,
                2.73950275e-05,
                1.80264033e-05,
                1.17866001e-05,
                7.65862168e-06,
                4.94578004e-06,
                3.17453177e-06,
                2.02545661e-06,
                1.28469904e-06,
                8.10123071e-07,
                5.07933514e-07,
                3.16667384e-07,
                1.96324242e-07,
                1.21046283e-07,
                7.42280310e-08,
                4.52746149e-08,
                2.74689782e-08,
                1.65791113e-08,
                9.95497454e-09,
                5.94713939e-09,
                3.53502926e-09,
                2.09084715e-09,
                1.23061882e-09,
                7.20809863e-10,
                4.20189085e-10,
                2.43787100e-10,
                1.40785767e-10,
                8.09258984e-11,
                4.63085408e-11,
                2.63768086e-11,
                1.49593924e-11,
                8.44425616e-12,
                4.74686214e-12,
                2.65644329e-12,
                1.47734550e-12,
                8.21530475e-13,
                4.50690158e-13,
                2.48123448e-13,
                1.32512251e-13,
                7.35980830e-14,
                4.48099043e-14,
                3.30421545e-14,
                2.40628013e-14,
                1.74025653e-14,
                1.25292347e-14,
                8.98931442e-15,
                6.42978942e-15,
                4.58562207e-15,
                3.26099157e-15,
                2.31235810e-15,
                1.63499398e-15,
                1.15275274e-15,
                8.10439886e-16,
                5.68169377e-16,
                3.97206490e-16,
                2.76914874e-16,
                1.92520944e-16,
                1.33481735e-16,
                9.22971240e-17,
                6.36483340e-17,
                4.37752100e-17,
                3.00276390e-17,
                2.05435232e-17,
                1.40183891e-17,
                9.29269070e-18,
            ]
        ]
    ]
    cres3 = np.array(cres3)

    creslist = [cres1, cres2, cres3]

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
