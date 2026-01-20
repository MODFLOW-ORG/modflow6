"""
Test the dispersion schemes in the gwt dispersion package for a one-dimensional
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

cases = ["dsp02a", "dsp02b"]
xt3d = [True, False]


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

    # constant heads on left and right so there is no flow
    c = {0: [((0, 0), 0.0000000), ((0, ncpl - 1), 0.0000000)]}

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
    cncs = {0: [[(0, 0), 1.0]]}
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
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
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
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # time step 200.
    cres1 = [
        [
            [
                1.0,
                0.97843231,
                0.97001134,
                0.95317805,
                0.94476996,
                0.9279709,
                0.91958413,
                0.90283612,
                0.89447905,
                0.87779875,
                0.86947965,
                0.85288352,
                0.84461056,
                0.82811481,
                0.81989601,
                0.80351654,
                0.79535977,
                0.77911214,
                0.77102507,
                0.75492445,
                0.74691457,
                0.73097567,
                0.72305022,
                0.70728731,
                0.69945329,
                0.68388011,
                0.67614426,
                0.66077401,
                0.65314277,
                0.63798806,
                0.63046759,
                0.61554043,
                0.60813657,
                0.59344832,
                0.5861666,
                0.57172796,
                0.56457354,
                0.55039454,
                0.54337227,
                0.52946222,
                0.52257658,
                0.50894408,
                0.50219919,
                0.48885212,
                0.48225173,
                0.46919723,
                0.46274472,
                0.44998918,
                0.44368756,
                0.43123664,
                0.42508856,
                0.41294716,
                0.40695489,
                0.39512719,
                0.38929262,
                0.37778204,
                0.37210672,
                0.36091599,
                0.35540108,
                0.3445322,
                0.33917854,
                0.3286328,
                0.32344087,
                0.3132189,
                0.30818885,
                0.2982906,
                0.29342225,
                0.28384706,
                0.2791399,
                0.26988646,
                0.26533971,
                0.25640614,
                0.25201871,
                0.24340255,
                0.23917307,
                0.23087133,
                0.22679818,
                0.21880735,
                0.21488866,
                0.20720475,
                0.20343841,
                0.19605698,
                0.19244068,
                0.18535688,
                0.18188808,
                0.17509667,
                0.17177266,
                0.16526803,
                0.16208593,
                0.15586217,
                0.15281894,
                0.14686982,
                0.14396229,
                0.13828133,
                0.13550621,
                0.13008669,
                0.12744057,
                0.12227559,
                0.11975497,
                0.11483745,
                0.11243875,
                0.10776147,
                0.10548106,
                0.10103668,
                0.09887086,
                0.09465197,
                0.09259703,
                0.08859614,
                0.08664832,
                0.08285793,
                0.08101347,
                0.07742605,
                0.07568121,
                0.07228925,
                0.07064027,
                0.06743631,
                0.06587947,
                0.06285607,
                0.06138769,
                0.05853752,
                0.05715393,
                0.05446974,
                0.05316736,
                0.050642,
                0.04941727,
                0.04704373,
                0.04589318,
                0.04366456,
                0.04258479,
                0.04049435,
                0.03948204,
                0.03752321,
                0.03657509,
                0.03474147,
                0.0338544,
                0.03213974,
                0.03131065,
                0.02970893,
                0.02893484,
                0.02744019,
                0.02671823,
                0.02532501,
                0.02465239,
                0.02335515,
                0.02272921,
                0.02152269,
                0.02094085,
                0.01982002,
                0.01927982,
                0.01823984,
                0.01773891,
                0.01677517,
                0.01631126,
                0.01541935,
                0.0149903,
                0.01416602,
                0.01376979,
                0.01300915,
                0.0126438,
                0.01194301,
                0.01160669,
                0.01096219,
                0.01065317,
                0.01006157,
                0.00977821,
                0.00923636,
                0.00897713,
                0.00848205,
                0.0082455,
                0.00779441,
                0.00757922,
                0.00716954,
                0.00697446,
                0.0066038,
                0.00642767,
                0.00609383,
                0.0059356,
                0.00563654,
                0.00549525,
                0.00522914,
                0.00510391,
                0.00486909,
                0.00475913,
                0.0045541,
                0.00445871,
                0.00428216,
                0.00420072,
                0.00405151,
                0.0039835,
                0.00386064,
                0.0038056,
                0.00370829,
                0.00366586,
                0.00359344,
                0.00356334,
                0.00351533,
                0.00349735,
                0.00346792,
            ]
        ]
    ]
    cres1 = np.array(cres1)

    cres2 = [
        [
            [
                1.0,
                0.9789382,
                0.97051702,
                0.95368454,
                0.94527607,
                0.92847748,
                0.92009017,
                0.90334227,
                0.89498447,
                0.87830394,
                0.86998394,
                0.85338722,
                0.84511319,
                0.82861649,
                0.82039646,
                0.80401567,
                0.79585751,
                0.7796082,
                0.77151959,
                0.75541693,
                0.74740534,
                0.73146407,
                0.72353676,
                0.70777112,
                0.69993509,
                0.68435885,
                0.67662084,
                0.6612472,
                0.65361367,
                0.63845524,
                0.63093235,
                0.61600116,
                0.60859476,
                0.59390218,
                0.58661778,
                0.57217453,
                0.56501733,
                0.55083343,
                0.54380828,
                0.52989307,
                0.52300445,
                0.50936655,
                0.50261858,
                0.48926588,
                0.48266233,
                0.46960198,
                0.46314624,
                0.45038465,
                0.44407974,
                0.43162258,
                0.42547115,
                0.41332335,
                0.40732768,
                0.39549342,
                0.38965542,
                0.37813816,
                0.37245936,
                0.36126184,
                0.35574344,
                0.34486766,
                0.3395105,
                0.32895779,
                0.32376235,
                0.31353335,
                0.30849978,
                0.29859447,
                0.29372261,
                0.28414033,
                0.27942969,
                0.27016915,
                0.26561893,
                0.25667827,
                0.25228741,
                0.24366417,
                0.23943131,
                0.23112253,
                0.22704603,
                0.21904822,
                0.21512624,
                0.20743541,
                0.20366584,
                0.19627758,
                0.19265809,
                0.18556756,
                0.18209565,
                0.17529761,
                0.17197055,
                0.16545942,
                0.16227435,
                0.1560442,
                0.15299808,
                0.14704271,
                0.14413238,
                0.13844531,
                0.13566746,
                0.130242,
                0.12759323,
                0.12242247,
                0.11989929,
                0.11497616,
                0.11257499,
                0.10789226,
                0.10560948,
                0.10115983,
                0.09899173,
                0.09476774,
                0.0927106,
                0.0887048,
                0.08675489,
                0.08295976,
                0.0811133,
                0.07752134,
                0.07577458,
                0.07237826,
                0.07072746,
                0.06751932,
                0.06596074,
                0.06293336,
                0.06146332,
                0.05860935,
                0.0572242,
                0.05453638,
                0.05323252,
                0.05070371,
                0.0494776,
                0.04710077,
                0.04594892,
                0.04371719,
                0.04263619,
                0.04054281,
                0.03952934,
                0.03756773,
                0.03661854,
                0.03478229,
                0.03389421,
                0.0321771,
                0.03134706,
                0.02974302,
                0.02896805,
                0.02747124,
                0.02674846,
                0.02535321,
                0.02467983,
                0.02338069,
                0.02275405,
                0.02154576,
                0.02096328,
                0.0198408,
                0.0193,
                0.0182585,
                0.01775702,
                0.01679186,
                0.01632745,
                0.01543423,
                0.01500472,
                0.01417922,
                0.01378257,
                0.01302081,
                0.01265508,
                0.01195326,
                0.0116166,
                0.01097115,
                0.01066181,
                0.01006936,
                0.00978571,
                0.00924308,
                0.00898359,
                0.00848779,
                0.00825101,
                0.00779927,
                0.00758388,
                0.00717361,
                0.00697834,
                0.00660715,
                0.00643086,
                0.00609654,
                0.00593817,
                0.00563869,
                0.00549728,
                0.0052308,
                0.00510545,
                0.0048703,
                0.00476025,
                0.00455494,
                0.00445946,
                0.00428267,
                0.00420117,
                0.00405175,
                0.00398368,
                0.00386065,
                0.00380557,
                0.00370812,
                0.00366566,
                0.00359314,
                0.00356302,
                0.00351495,
                0.00349695,
                0.00346702,
            ]
        ]
    ]
    cres2 = np.array(cres2)

    creslist = [cres1, cres2]

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
