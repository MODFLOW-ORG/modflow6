"""
Tests the starting saturation of a GWT model that simulates transport on a
subset of the GWF model cells.

A one-row flow model with convertible cells starts from a sloping water table
and drains to a constant head at the downgradient end. The first column is
separated from the rest by a horizontal flow barrier with a zero hydraulic
characteristic and held by its own constant head, so no water crosses it. A
full-domain transport model and a transport model that excludes the first
column both simulate linear sorption, which uses the saturation at the start of
the first time step.

Cases:
  - idomsat : concentrations in the reduced transport domain, whose node numbers
              are offset from the flow model node numbers, must match the
              full-domain model.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["idomsat"]

nlay, nrow, ncol = 1, 1, 12
delr, delc = 10.0, 10.0
top, botm = 10.0, [0.0]
hk = 1.0
ss, sy = 1.0e-5, 0.2
porosity = 0.3
bulk_density = 1600.0
distcoef = 1.0e-4
perlen, nstp = 10.0, 10

# sloping water table, so the starting saturation differs between cells
strt = np.linspace(9.0, 4.0, ncol)
hisolated = 5.0
hout = 1.0


def build_flow_model(sim, gwfname):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-10,
        outer_maximum=100,
        inner_maximum=300,
        inner_dvclose=1e-10,
        rcloserecord=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    h0 = strt.copy()
    h0[0] = hisolated
    flopy.mf6.ModflowGwfic(gwf, strt=h0.reshape((nlay, nrow, ncol)))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk, save_specific_discharge=True)
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=ss, sy=sy, transient={0: True}, save_flows=True
    )
    # no water crosses the barrier between the first two columns
    flopy.mf6.ModflowGwfhfb(
        gwf, maxhfb=1, stress_period_data={0: [((0, 0, 0), (0, 0, 1), 0.0)]}
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [((0, 0, 0), hisolated), ((0, 0, ncol - 1), hout)]},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL")],
    )


def build_transport_model(sim, gwtname, idomain):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-10,
        outer_maximum=100,
        inner_maximum=300,
        inner_dvclose=1e-10,
        rcloserecord=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(ims, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=1.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        sorption="linear",
        bulk_density=bulk_density,
        distcoef=distcoef,
    )
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )

    gwfname = "gwf_" + name
    build_flow_model(sim, gwfname)

    idomain_sub = np.ones((nlay, nrow, ncol), dtype=int)
    idomain_sub[0, 0, 0] = 0
    for gwtname, idomain in (
        ("gwtfull_" + name, np.ones((nlay, nrow, ncol), dtype=int)),
        ("gwtsub_" + name, idomain_sub),
    ):
        build_transport_model(sim, gwtname, idomain)
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            filename=f"{gwtname}.gwfgwt",
        )

    return sim, None


def check_output(idx, test):
    name = cases[idx]

    conc = {}
    for key in ("full", "sub"):
        fpth = os.path.join(test.workspace, f"gwt{key}_{name}.ucn")
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc[key] = cobj.get_alldata()

    # the cells drain, so sorbate released into the remaining water raises the
    # concentration and the starting saturation matters
    assert conc["full"][0, 0, 0, 1:].max() > 1.0 + 1e-3, (
        "concentrations did not change in the first time step"
    )

    assert np.allclose(
        conc["sub"][:, 0, 0, 1:], conc["full"][:, 0, 0, 1:], rtol=0.0, atol=1e-8
    ), (
        "concentrations in the reduced transport domain do not match the "
        "full-domain solution"
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
