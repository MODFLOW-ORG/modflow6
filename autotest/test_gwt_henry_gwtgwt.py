import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "henry01-gwtgwt-ups",
    "henry01-gwtgwt-cen",
    "henry01-gwtgwt-tvd",
    "henry01-gwtgwt-utvd",
]
advection_scheme = ["UPSTREAM", "CENTRAL", "TVD", "UTVD"]

lx = 2.0
lz = 1.0

nlay = 10
nrow = 1
ncol = 20
ncol_sub = 10
nper = 1
delr = lx / ncol
delc = 1.0
top = 1.0
delz = lz / nlay
botm = list(top - np.arange(delz, nlay * delz + delz, delz))

perlen = [0.5]
nstp = [500]
tsmult = [1.0]
steady = [True]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-9, 0.97

# tolerance for concentrations, relate to hclose and rclose
conc_tol = 1e-06

hc = 864.0


def get_gwf_model(sim, model_shape, model_desc):
    nlay, nrow, ncol = model_shape

    # create gwf model
    gwfname = "gwf_" + model_desc
    gwtname = "gwt_" + model_desc

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)

    xoff = 0.0
    yoff = 0.0
    if model_desc == "right":
        xoff = ncol * delr

    _ = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xoff,
        yorigin=yoff,
    )

    # initial conditions
    _ = flopy.mf6.ModflowGwfic(gwf, strt=1.0)

    # node property flow
    _ = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=hc,
    )

    pd = [(0, 0.7, 0.0, gwtname, "none")]
    _ = flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    def chd_value(k):
        # depth = k * delz + 0.5 * delz
        # hf = top + 0.025 * depth
        hf = top
        return hf

    # chd files for right model
    if not model_desc == "left":
        chdlist1 = []
        for k in range(nlay):
            chdlist1.append([(k, 0, ncol - 1), chd_value(k), 35.0])
        _ = flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist1,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="CHD-1",
            auxiliary="CONCENTRATION",
            filename=f"{gwfname}.chd",
        )

    # WEL for left model
    if not model_desc == "right":
        wellist1 = []
        qwell = 5.7024 / nlay
        for k in range(nlay):
            wellist1.append([(k, 0, 0), qwell, 0.0])
        _ = flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=wellist1,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="WEL-1",
            auxiliary="CONCENTRATION",
            filename=f"{gwfname}.wel",
        )

    # output control
    _ = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def get_gwt_model(sim, model_shape, model_desc, adv_scheme):
    nlay, nrow, ncol = model_shape

    # create gwf model
    gwtname = "gwt_" + model_desc

    # create gwt model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    xoff = 0.0
    yoff = 0.0
    if model_desc == "right":
        xoff = ncol * delr

    _ = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xoff,
        yorigin=yoff,
    )

    # initial conditions
    _ = flopy.mf6.ModflowGwtic(gwt, strt=35.0, filename=f"{gwtname}.ic")

    # advection
    _ = flopy.mf6.ModflowGwtadv(gwt, scheme=adv_scheme, filename=f"{gwtname}.adv")

    # dispersion
    diffc = 0.57024
    _ = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        diffc=diffc,
        # alh=0., alv=0., ath=0., atv=0.,
        filename=f"{gwtname}.dsp",
    )

    # mass storage and transfer
    porosity = 0.35
    _ = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")

    # sources
    if model_desc == "right":
        sourcerecarray = [
            ("CHD-1", "AUX", "CONCENTRATION"),
        ]
    elif model_desc == "left":
        sourcerecarray = [
            ("WEL-1", "AUX", "CONCENTRATION"),
        ]
    elif model_desc == "ref":
        sourcerecarray = [
            ("WEL-1", "AUX", "CONCENTRATION"),
            ("CHD-1", "AUX", "CONCENTRATION"),
        ]

    _ = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm")

    # output control
    _ = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return gwt


def build_models(idx, test):
    name = cases[idx]
    print("RUINNING: ", name, advection_scheme[idx])

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create flow models and GWF-GWF exchange
    gwf_ref = get_gwf_model(sim, (nlay, nrow, ncol), "ref")
    gwf_left = get_gwf_model(sim, (nlay, nrow, ncol_sub), "left")
    gwf_right = get_gwf_model(sim, (nlay, nrow, ncol_sub), "right")

    imsgwf_ref = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
        filename="gwf.ims",
    )
    sim.register_ims_package(imsgwf_ref, [gwf_ref.name, gwf_left.name, gwf_right.name])

    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (ilay, 0, ncol_sub - 1),
            (ilay, 0, 0),
            1,
            0.5 * delr,
            0.5 * delr,
            delc,
            angldegx,
            cdist,
        ]
        for ilay in range(nlay)
    ]
    _ = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwf_left.name,
        exgmnameb=gwf_right.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="leftright.gwfgwf",
    )

    # create transport models and GWT-GWT exchange
    gwt_ref = get_gwt_model(sim, (nlay, nrow, ncol), "ref", advection_scheme[idx])
    gwt_left = get_gwt_model(sim, (nlay, nrow, ncol_sub), "left", advection_scheme[idx])
    gwt_right = get_gwt_model(
        sim, (nlay, nrow, ncol_sub), "right", advection_scheme[idx]
    )

    imsgwt_ref = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
        filename="gwt.ims",
    )
    sim.register_ims_package(imsgwt_ref, [gwt_ref.name, gwt_left.name, gwt_right.name])

    _ = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1=gwf_left.name,
        gwfmodelname2=gwf_right.name,
        adv_scheme=advection_scheme[idx],
        nexg=len(gwfgwf_data),
        exgmnamea=gwt_left.name,
        exgmnameb=gwt_right.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="leftright.gwtgwt",
    )

    # Finally connect flow to transport: GWF GWT exchange
    _ = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_ref.name,
        exgmnameb=gwt_ref.name,
        filename="ref.gwfgwt",
    )

    _ = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_left.name,
        exgmnameb=gwt_left.name,
        filename="left.gwfgwt",
    )

    _ = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_right.name,
        exgmnameb=gwt_right.name,
        filename="right.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, "gwf_ref.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()

    fpth = os.path.join(test.workspace, "gwf_left.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_left = hds.get_data()

    fpth = os.path.join(test.workspace, "gwf_right.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_right = hds.get_data()

    heads_gwfgwf = np.append(heads_left, heads_right, axis=2)

    # compare heads
    maxdiff = np.amax(abs(heads - heads_gwfgwf))
    assert maxdiff < 10 * hclose, (
        f"Max. head diff. {maxdiff} should \
        be within solver tolerance (x10): {10 * hclose}"
    )

    fpth = os.path.join(test.workspace, "gwt_ref.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc_ref = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    fpth = os.path.join(test.workspace, "gwt_left.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc_left = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    fpth = os.path.join(test.workspace, "gwt_right.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc_right = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # merge left and right concentrations for bottom layer:
    conc_gwtgwt = np.append(conc_left, conc_right, axis=2)

    maxdiff = np.amax(abs(conc_gwtgwt - conc_ref))
    assert maxdiff < conc_tol, (
        f"Max. concentration diff. {maxdiff} should \
        be within solver tolerance (x10): {conc_tol}"
    )


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
