"""
General test for the interface model approach.
It compares the result of a single reference model
to the equivalent case where the domain is decomposed
and joined by a GWF-GWF exchange.

       'refmodel'              'leftmodel'     'rightmodel'

   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1    VS    1 1 1 1 1   +   1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1

We assert equality on the head values and the (components of)
specific discharges. All models are part of the same solution
for convenience. Finally, the budget error is checked.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ifmod_buy01"]

# some global convenience...:
# model names
mname_ref = "refmodel"
mname_gwtref = "refgwtmodel"
mname_left = "leftmodel"
mname_gwtleft = "leftgwtmodel"
mname_right = "rightmodel"
mname_gwtright = "rightgwtmodel"

# solver criterion
hclose_check = 1e-9

# model spatial discretization
nlay = 10
ncol = 10
ncol_left = 5
ncol_right = 5
nrow = 10

# cell spacing
delr = 10.0
delc = 10.0
area = delr * delc

# shift (hor. and vert.)
shift_some_x = -20 * delr  # avoids overlap
shift_x = 5 * delr
shift_y = 0.0

# top/bot of the aquifer
tops = [0.0 - ilay * 5.0 for ilay in range(nlay + 1)]

# hydraulic conductivity
k11 = 10.0

# boundary stress period data
h_left = -1.0
h_right = -2.0

# initial head
h_start = -2.0

# head boundaries
left_chd = [[(ilay, irow, 0), h_left] for irow in range(nrow) for ilay in range(nlay)]
right_chd = [
    [(ilay, irow, ncol - 1), h_right] for irow in range(nrow) for ilay in range(nlay)
]
right_chd_split = [
    [(ilay, irow, ncol_right - 1), h_right]
    for irow in range(nrow)
    for ilay in range(nlay)
]

chd_data = left_chd + right_chd
chd_spd = {0: chd_data}
chd_spd_left = {0: left_chd}
chd_spd_right = {0: right_chd_split}

# initial conc
c_strt = 1.1
porosity = 0.30


def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = hclose_check, 1e-3, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
        filename="gwf.ims",
    )

    # the full gwf model as a reference
    add_refmodel(sim)

    # now add two coupled models with the interface model enabled,
    # to be stored in the same solution as the reference model
    add_leftmodel(sim)
    add_rightmodel(sim)
    add_gwfexchange(sim)

    # the transport sector
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="gwt.ims",
    )

    gwt = add_gwtrefmodel(sim)
    gwtleft = add_gwtleftmodel(sim)
    gwtright = add_gwtrightmodel(sim)
    sim.register_ims_package(imsgwt, [gwt.name, gwtleft.name, gwtright.name])

    add_gwtexchange(sim)

    # add GWF GWT exchanges
    gwfgwt_ref = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=mname_ref,
        exgmnameb=mname_gwtref,
        filename="reference.gwfgwt",
    )
    gwfgwt_left = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=mname_left,
        exgmnameb=mname_gwtleft,
        filename="left.gwfgwt",
    )
    gwfgwt_right = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=mname_right,
        exgmnameb=mname_gwtright,
        filename="right.gwfgwt",
    )

    return sim


def add_refmodel(sim):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_ref, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        xorigin=shift_some_x,
        yorigin=0.0,
        top=tops[0],
        botm=tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )

    # buy
    pd = [(0, 0.7, 0.0, mname_gwtref, "CONCENTRATION")]
    buy = flopy.mf6.ModflowGwfbuy(
        gwf,
        packagedata=pd,
        denseref=1000.0,
    )

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_ref}.hds",
        budget_filerecord=f"{mname_ref}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_leftmodel(sim):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_left, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_left,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    pd = [(0, 0.7, 0.0, mname_gwtleft, "CONCENTRATION")]
    buy = flopy.mf6.ModflowGwfbuy(
        gwf,
        packagedata=pd,
        denseref=1000.0,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_left)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_left}.hds",
        budget_filerecord=f"{mname_left}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_rightmodel(sim):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_right, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_right,
        delr=delr,
        delc=delc,
        xorigin=shift_x,
        yorigin=shift_y,
        top=tops[0],
        botm=tops[1:],
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    pd = [(0, 0.7, 0.0, mname_gwtright, "CONCENTRATION")]
    buy = flopy.mf6.ModflowGwfbuy(
        gwf,
        packagedata=pd,
        denseref=1000.0,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_right)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_right}.hds",
        budget_filerecord=f"{mname_right}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_gwfexchange(sim):
    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (ilay, irow, ncol_left - 1),
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
        exgmnamea=mname_left,
        exgmnameb=mname_right,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
    )


def add_gwtrefmodel(sim):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=mname_gwtref)

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=c_strt)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # storage
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{mname_gwtref}.cbc",
        concentration_filerecord=f"{mname_gwtref}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    fmi = flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)

    return gwt


def add_gwtleftmodel(sim):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=mname_gwtleft)

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_left,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=c_strt)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # storage
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{mname_gwtleft}.cbc",
        concentration_filerecord=f"{mname_gwtleft}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    fmi = flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)

    return gwt


def add_gwtrightmodel(sim):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=mname_gwtright)

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_right,
        delr=delr,
        delc=delc,
        xorigin=shift_x,
        yorigin=shift_y,
        top=tops[0],
        botm=tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=c_strt)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # storage
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{mname_gwtright}.cbc",
        concentration_filerecord=f"{mname_gwtright}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    fmi = flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)

    return gwt


def add_gwtexchange(sim):
    angldegx = 0.0
    cdist = delr
    gwtgwt_data = [
        [
            (ilay, irow, ncol_left - 1),
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
    gwtgwt = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        nexg=len(gwtgwt_data),
        gwfmodelname1=mname_left,
        gwfmodelname2=mname_right,
        exgmnamea=mname_gwtleft,
        exgmnameb=mname_gwtright,
        exchangedata=gwtgwt_data,
        auxiliary=["ANGLDEGX", "CDIST"],
    )


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def qxqyqz(fname, nlay, nrow, ncol):
    nodes = nlay * nrow * ncol
    cbb = flopy.utils.CellBudgetFile(fname, precision="double")
    spdis = cbb.get_data(text="DATA-SPDIS")[0]
    qx = np.ones((nodes), dtype=float) * 1.0e30
    qy = np.ones((nodes), dtype=float) * 1.0e30
    qz = np.ones((nodes), dtype=float) * 1.0e30
    n0 = spdis["node"] - 1
    qx[n0] = spdis["qx"]
    qy[n0] = spdis["qy"]
    qz[n0] = spdis["qz"]
    qx = qx.reshape(nlay, nrow, ncol)
    qy = qy.reshape(nlay, nrow, ncol)
    qz = qz.reshape(nlay, nrow, ncol)
    qx = np.ma.masked_equal(qx, 1.0e30)
    qy = np.ma.masked_equal(qy, 1.0e30)
    qz = np.ma.masked_equal(qz, 1.0e30)
    return qx, qy, qz


def check_output(idx, test):
    print("comparing heads and spec. discharge to single model reference...")

    fpth = os.path.join(test.workspace, f"{mname_ref}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_ref}.cbc")
    nlay, nrow, ncol = heads.shape
    qxb, qyb, qzb = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{mname_left}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_left = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_left}.cbc")
    nlay, nrow, ncol = heads_left.shape
    qxb_left, qyb_left, qzb_left = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{mname_right}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_right = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_right}.cbc")
    nlay, nrow, ncol = heads_right.shape
    qxb_right, qyb_right, qzb_right = qxqyqz(fpth, nlay, nrow, ncol)

    heads_2models = np.append(heads_left[0], heads_right[0], axis=1)

    # compare heads
    maxdiff = np.amax(abs(heads - heads_2models))
    assert maxdiff < 10 * hclose_check, (
        f"Max. head diff. {maxdiff} should \
                     be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_x left
    maxdiff = np.amax(abs(qxb[:, :, 0:5] - qxb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (x) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_y left
    maxdiff = np.amax(abs(qyb[:, :, 0:5] - qyb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (y) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_z left
    maxdiff = np.amax(abs(qzb[:, :, 0:5] - qzb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (z) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_x right
    maxdiff = np.amax(abs(qxb[:, :, 5:] - qxb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (x) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_y right
    maxdiff = np.amax(abs(qyb[:, :, 5:] - qyb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (y) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_z right
    maxdiff = np.amax(abs(qzb[:, :, 5:] - qzb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (z) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # check budget error from .lst file
    for mname in [mname_ref, mname_left, mname_right]:
        fpth = os.path.join(test.workspace, f"{mname}.lst")
        for line in open(fpth):
            if line.lstrip().startswith("PERCENT"):
                cumul_balance_error = float(line.split()[3])
                assert abs(cumul_balance_error) < 0.00001, (
                    f"Cumulative balance error = {cumul_balance_error} for {mname}, "
                    "should equal 0.0"
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
