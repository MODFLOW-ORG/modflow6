"""
Test the advection schemes in the gwt advection package for a one-dimensional
model grid of square cells.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["adv01a_gwtgwt", "adv01b_gwtgwt", "adv01c_gwtgwt"]
scheme = ["upstream", "central", "tvd"]
gdelr = 1.0

# solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def get_gwf_model(sim, gwfname, gwfpath, modelshape, chdspd=None, welspd=None):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = gdelr
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xshift,
        yorigin=yshift,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=laytyp,
        k=hk,
        save_specific_discharge=True,
    )

    # chd files
    if chdspd is not None:
        chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
            gwf,
            # maxbound=len(c),
            stress_period_data=chdspd,
            save_flows=False,
            pname="CHD-1",
        )

    # wel files
    if welspd is not None:
        wel = flopy.mf6.ModflowGwfwel(
            gwf,
            print_input=True,
            print_flows=True,
            # maxbound=len(w),
            stress_period_data=welspd,
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
    gwf.set_model_relative_path(gwfpath)
    return gwf


def get_gwt_model(sim, gwtname, gwtpath, modelshape, scheme, sourcerecarray=None):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
    )
    gwt.name_file.save_flows = True

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=gdelr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xshift,
        yorigin=yshift,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme)

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # obs_data = {
    #     "conc_obs.csv": [
    #         ("(1-1-10)", "CONCENTRATION", (0, 0, 9)),
    #         ("(1-1-50)", "CONCENTRATION", (0, 0, 49)),
    #     ],
    #     "flow_obs.csv": [
    #         ("c10-c11", "FLOW-JA-FACE", (0, 0, 9), (0, 0, 10)),
    #         ("c50-c51", "FLOW-JA-FACE", (0, 0, 49), (0, 0, 50)),
    #         ("c99-c100", "FLOW-JA-FACE", (0, 0, 98), (0, 0, 99)),
    #     ],
    # }
    #
    # obs_package = flopy.mf6.ModflowUtlobs(
    #     gwt,
    #     pname="conc_obs",
    #     filename="{}.obs".format(gwtname),
    #     digits=10,
    #     print_input=True,
    #     continuous=obs_data,
    # )

    gwt.set_model_relative_path(gwtpath)
    return gwt


def build_models(idx, test):
    # temporal discretization
    nper = 1
    perlen = [5.0]
    nstp = [200]
    tsmult = [1.0]
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=ws, version="mf6", exe_name="mf6", sim_ws=ws)
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=tdis_rc, pname="sim.tdis"
    )

    # grid information
    nlay, nrow, ncol = 1, 1, 50

    # Create gwf1 model
    welspd = {0: [[(0, 0, 0), 1.0, 1.0]]}
    chdspd = None  # {0: [[(0, 0, 99), 0.0000000]]}
    gwf1 = get_gwf_model(
        sim,
        "flow1",
        "flow1",
        (nlay, nrow, ncol, 0.0, 0.0),
        chdspd=chdspd,
        welspd=welspd,
    )

    # Create gwf2 model
    welspd = None  # {0: [[(0, 0, 0), 1.0, 1.0]]}
    chdspd = {0: [[(0, 0, ncol - 1), 0.0000000]]}
    gwf2 = get_gwf_model(
        sim,
        "flow2",
        "flow2",
        (nlay, nrow, ncol, 50.0 * gdelr, 0.0),
        chdspd=chdspd,
        welspd=welspd,
    )

    # gwf-gwf with interface model enabled
    gwfgwf_data = [[(0, 0, ncol - 1), (0, 0, 0), 1, 0.5, 0.5, 1.0, 0.0, 1.0]]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwf1.name,
        exgmnameb=gwf2.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="flow1_flow2.gwfgwf",
        dev_interfacemodel_on=True,
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
        filename="flow.ims",
    )
    sim.register_ims_package(imsgwf, [gwf1.name, gwf2.name])

    # Create gwt model
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    gwt1 = get_gwt_model(
        sim,
        "transport1",
        "transport1",
        (nlay, nrow, ncol, 0.0, 0.0),
        scheme[idx],
        sourcerecarray=sourcerecarray,
    )

    # Create gwt model
    sourcerecarray = None
    gwt2 = get_gwt_model(
        sim,
        "transport2",
        "transport2",
        (nlay, nrow, ncol, 50.0 * gdelr, 0.0),
        scheme[idx],
        sourcerecarray=sourcerecarray,
    )

    # Create GWT GWT exchange
    gwt1gwt2 = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1=gwf1.name,
        gwfmodelname2=gwf2.name,
        adv_scheme=scheme[idx],
        nexg=len(gwfgwf_data),
        exgmnamea=gwt1.name,
        exgmnameb=gwt2.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="transport1_transport2.gwtgwt",
    )

    # GWF GWT exchange
    gwfgwt1 = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="flow1",
        exgmnameb="transport1",
        filename="flow1_transport1.gwfgwt",
    )
    gwfgwt2 = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="flow2",
        exgmnameb="transport2",
        filename="flow2_transport2.gwfgwt",
    )

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
        filename="transport.ims",
    )
    sim.register_ims_package(imsgwt, [gwt1.name, gwt2.name])

    return sim, None


def check_output(idx, test):
    gwtname = "transport1"

    fpth = os.path.join(test.workspace, gwtname, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc1 = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    gwtname = "transport2"

    fpth = os.path.join(test.workspace, gwtname, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc2 = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    conc = np.append(conc1, conc2)

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
                9.99999999e-01,
                9.99999997e-01,
                9.99999991e-01,
                9.99999971e-01,
                9.99999914e-01,
                9.99999761e-01,
                9.99999372e-01,
                9.99998435e-01,
                9.99996286e-01,
                9.99991577e-01,
                9.99981712e-01,
                9.99961893e-01,
                9.99923632e-01,
                9.99852532e-01,
                9.99725120e-01,
                9.99504599e-01,
                9.99135431e-01,
                9.98536850e-01,
                9.97595635e-01,
                9.96158712e-01,
                9.94026505e-01,
                9.90948130e-01,
                9.86619748e-01,
                9.80687319e-01,
                9.72754814e-01,
                9.62398489e-01,
                9.49187176e-01,
                9.32707801e-01,
                9.12594513e-01,
                8.88559134e-01,
                8.60420154e-01,
                8.28127324e-01,
                7.91779115e-01,
                7.51630867e-01,
                7.08092322e-01,
                6.61714306e-01,
                6.13165405e-01,
                5.63200494e-01,
                5.12623768e-01,
                4.62249349e-01,
                4.12862664e-01,
                3.65185517e-01,
                3.19847250e-01,
                2.77363614e-01,
                2.38124183e-01,
                2.02388273e-01,
                1.70288648e-01,
                1.41841739e-01,
                1.16962748e-01,
                9.54838854e-02,
                7.71740354e-02,
                6.17583229e-02,
                4.89363652e-02,
                3.83983188e-02,
                2.98381826e-02,
                2.29641338e-02,
                1.75059339e-02,
                1.32196416e-02,
                9.89000005e-03,
                7.33093269e-03,
                5.38459977e-03,
                3.91944360e-03,
                2.82760119e-03,
                2.02199855e-03,
                1.43337156e-03,
                1.00739149e-03,
                7.02013580e-04,
                4.85116958e-04,
                3.32465664e-04,
                2.25991387e-04,
                1.52379541e-04,
                1.01928496e-04,
                6.76460984e-05,
                4.45462926e-05,
                2.91101871e-05,
                1.88792800e-05,
                1.21527525e-05,
                7.76522212e-06,
                4.92565188e-06,
                3.10201677e-06,
                1.93969988e-06,
                1.20440812e-06,
                7.42676511e-07,
                4.54831064e-07,
                2.76669882e-07,
                1.67174989e-07,
                1.00349240e-07,
                5.98446532e-08,
                3.54600737e-08,
            ]
        ]
    ]
    cres1 = np.array(cres1)

    cres2 = [
        [
            [
                9.99996617e-01,
                1.00001184e00,
                1.00000294e00,
                9.99972914e-01,
                9.99992627e-01,
                1.00004237e00,
                1.00002081e00,
                9.99945149e-01,
                9.99952654e-01,
                1.00005669e00,
                1.00008810e00,
                9.99966402e-01,
                9.99865541e-01,
                9.99967791e-01,
                1.00015792e00,
                1.00014755e00,
                9.99895530e-01,
                9.99724106e-01,
                9.99916592e-01,
                1.00029941e00,
                1.00038455e00,
                9.99960678e-01,
                9.99433053e-01,
                9.99453350e-01,
                1.00018163e00,
                1.00097923e00,
                1.00093550e00,
                9.99790199e-01,
                9.98371554e-01,
                9.98054584e-01,
                9.99598363e-01,
                1.00229288e00,
                1.00416575e00,
                1.00323035e00,
                9.98995210e-01,
                9.93234271e-01,
                9.89448228e-01,
                9.91206357e-01,
                1.00016889e00,
                1.01473298e00,
                1.02990960e00,
                1.03846239e00,
                1.03282855e00,
                1.00710727e00,
                9.58480908e-01,
                8.87726436e-01,
                7.98820097e-01,
                6.97900399e-01,
                5.91969549e-01,
                4.87686471e-01,
                3.90487541e-01,
                3.04127133e-01,
                2.30608327e-01,
                1.70400015e-01,
                1.22812141e-01,
                8.64138068e-02,
                5.94120233e-02,
                3.99463958e-02,
                2.62868102e-02,
                1.69426845e-02,
                1.07033555e-02,
                6.63198283e-03,
                4.03300421e-03,
                2.40844447e-03,
                1.41323306e-03,
                8.15254552e-04,
                4.62589305e-04,
                2.58303233e-04,
                1.42001900e-04,
                7.68911977e-05,
                4.10256980e-05,
                2.15775541e-05,
                1.11912143e-05,
                5.72578796e-06,
                2.89083689e-06,
                1.44073067e-06,
                7.09001789e-07,
                3.44624235e-07,
                1.65501321e-07,
                7.85475047e-08,
                3.68512253e-08,
                1.70949923e-08,
                7.84310280e-09,
                3.55966819e-09,
                1.59856594e-09,
                7.10467596e-10,
                3.12565151e-10,
                1.36146377e-10,
                5.87252052e-11,
                2.50886169e-11,
                1.06179506e-11,
                4.45237718e-12,
                1.85013624e-12,
                7.61982955e-13,
                3.11095972e-13,
                1.25908830e-13,
                5.05704707e-14,
                2.00370648e-14,
                8.15003576e-15,
                2.57563506e-15,
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
                9.99999999e-01,
                9.99999997e-01,
                9.99999991e-01,
                9.99999975e-01,
                9.99999926e-01,
                9.99999789e-01,
                9.99999407e-01,
                9.99998374e-01,
                9.99995665e-01,
                9.99988785e-01,
                9.99971918e-01,
                9.99932078e-01,
                9.99841550e-01,
                9.99643930e-01,
                9.99229970e-01,
                9.98398720e-01,
                9.96800070e-01,
                9.93857995e-01,
                9.88681096e-01,
                9.79978744e-01,
                9.66015902e-01,
                9.44652308e-01,
                9.13514114e-01,
                8.70328697e-01,
                8.13410724e-01,
                7.42224214e-01,
                6.57879960e-01,
                5.63390876e-01,
                4.63530320e-01,
                3.64233335e-01,
                2.71628522e-01,
                1.90935412e-01,
                1.25541007e-01,
                7.65316248e-02,
                4.28052252e-02,
                2.16851758e-02,
                9.78976172e-03,
                3.85613094e-03,
                1.28872611e-03,
                3.52070089e-04,
                7.49188445e-05,
                1.17688715e-05,
                1.33952025e-06,
                1.08174095e-07,
                -4.82019087e-08,
                -5.67180537e-08,
                -4.65251289e-08,
                -3.25511455e-08,
                -1.94644548e-08,
                -9.78876693e-09,
                -4.07380361e-09,
                -1.38097809e-09,
                -3.72934181e-10,
                -7.83508455e-11,
                -1.26040926e-11,
                -1.48260453e-12,
                4.10392230e-14,
                2.44993743e-13,
                2.46295025e-13,
                1.90964563e-13,
                1.03476379e-13,
                3.96502895e-14,
                1.04500247e-14,
                2.00830327e-15,
                4.70831032e-16,
                3.38440506e-16,
                2.49848438e-16,
                1.83245111e-16,
                1.32361223e-16,
                9.39406563e-17,
                6.54891851e-17,
                4.48667613e-17,
                3.02333440e-17,
                2.00567815e-17,
                1.31110206e-17,
                8.45177289e-18,
                5.37610069e-18,
                3.37597383e-18,
            ]
        ]
    ]
    cres3 = np.array(cres3)

    creslist = [cres1, cres2, cres3]

    assert np.allclose(creslist[idx], conc), (
        "simulated concentrations do not match with known solution."
    )

    # check budget
    for mname in ["transport1", "transport2"]:
        fpth = os.path.join(test.workspace, mname, f"{mname}.lst")
        for line in open(fpth):
            if line.lstrip().startswith("PERCENT"):
                cumul_balance_error = float(line.split()[3])
                assert abs(cumul_balance_error) < 0.00001, (
                    f"Cumulative balance error = {cumul_balance_error} for {mname}, "
                    "should equal 0.0"
                )

        # get grid data (from GWF)
        gwfname = "flow1" if mname == "transport1" else "flow2"
        fpth = os.path.join(test.workspace, gwfname, f"{gwfname}.dis.grb")
        grb = flopy.mf6.utils.MfGrdFile(fpth)

        # Check on residual, which is stored in diagonal position of
        # flow-ja-face.  Residual should be less than convergence tolerance,
        # or this means the residual term is not added correctly.
        fpth = os.path.join(test.workspace, mname, f"{mname}.cbc")
        cbb = flopy.utils.CellBudgetFile(fpth)
        flow_ja_face = cbb.get_data(text="FLOW-JA-FACE")
        ia = grb._datadict["IA"] - 1
        for fjf in flow_ja_face:
            fjf = fjf.flatten()
            res = fjf[ia[:-1]]
            errmsg = f"min or max flowja residual too large {res.min()} {res.max()}"
            # TODO: this is not implemented yet:
            # assert np.allclose(res, 0.0, atol=1.0e-6), errmsg


@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
