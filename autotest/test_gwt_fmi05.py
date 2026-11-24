"""
Tests that a GWT model reading flows and heads from file with FMI can
retry a failed time step with ATS.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwtfmi05a", "gwtfmi05b"]
nper = {"gwtfmi05a": 1, "gwtfmi05b": 2}

nlay, nrow, ncol = 1, 1, 10
delr = 100.0
delc = 1.0
top = 1.0
botm = [0.0]
perlen = 10.0
porosity = 0.3


def add_tdis(sim, nper, ats=False):
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(perlen, 1, 1.0)] * nper
    )
    if ats:
        # dt0, dtmin, dtmax, dtadj, dtfailadj
        tdis.ats.initialize(
            maxats=nper,
            perioddata=[(k, perlen, 0.1, perlen, 2.0, 5.0) for k in range(nper)],
            filename=f"{sim.name}.ats",
        )


def add_gwf(sim, name, files=False):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=files)
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
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=files, save_saturation=files)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, ncol - 1), 1.0]])
    flopy.mf6.ModflowGwfwel(
        gwf,
        pname="WEL-1",
        auxiliary=["CONCENTRATION"],
        stress_period_data=[[(0, 0, 0), 0.1, 1.0]],
    )
    if files:
        flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{name}.bud",
            head_filerecord=f"{name}.hds",
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
    return gwf


def add_gwt(sim, name, fmi_pd=None):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name)
    # with one outer iteration, the transport solution converges only if
    # the concentration change is small, so the first full-length step fails
    ims = flopy.mf6.ModflowIms(
        sim,
        outer_maximum=1,
        outer_dvclose=2e-2,
        inner_maximum=100,
        inner_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
        filename=f"{name}.ims",
    )
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(gwt, sources=[("WEL-1", "AUX", "CONCENTRATION")])
    if fmi_pd is not None:
        flopy.mf6.ModflowGwtfmi(gwt, packagedata=fmi_pd)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{name}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )
    return gwt, ims


def build_models(idx, test):
    name = test.name
    np_ = nper[name]
    mf6 = test.targets["mf6"]
    gwf_ws = test.workspace / "gwf"
    gwt_ws = test.workspace / "gwt"
    cpl_ws = test.workspace / "coupled"

    # flow model writing budget and head files
    gwf_sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwf_ws, exe_name=mf6)
    add_tdis(gwf_sim, np_)
    flopy.mf6.ModflowIms(gwf_sim)
    add_gwf(gwf_sim, "gwf", files=True)

    # transport model reading them with FMI
    gwt_sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwt_ws, exe_name=mf6)
    add_tdis(gwt_sim, np_, ats=True)
    pd = [
        ("GWFBUDGET", gwf_ws / "gwf.bud", None),
        ("GWFHEAD", gwf_ws / "gwf.hds", None),
    ]
    add_gwt(gwt_sim, "gwt", fmi_pd=pd)

    # same models coupled with an exchange
    cpl_sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=cpl_ws, exe_name=mf6)
    add_tdis(cpl_sim, np_, ats=True)
    ims_gwf = flopy.mf6.ModflowIms(cpl_sim, filename="gwf.ims")
    add_gwf(cpl_sim, "gwf")
    _, ims_gwt = add_gwt(cpl_sim, "gwt")
    cpl_sim.register_ims_package(ims_gwf, ["gwf"])
    cpl_sim.register_ims_package(ims_gwt, ["gwt"])
    flopy.mf6.ModflowGwfgwt(cpl_sim, exgmnamea="gwf", exgmnameb="gwt")

    return gwf_sim, gwt_sim, cpl_sim


def check_output(idx, test):
    gwt_ws = test.workspace / "gwt"
    cpl_ws = test.workspace / "coupled"

    lst = (gwt_ws / "mfsim.lst").read_text()
    assert "Failed solution for step 1 and period 1 will be retried" in lst, (
        "expected the first time step to be retried"
    )

    ucn = flopy.utils.HeadFile(gwt_ws / "gwt.ucn", text="CONCENTRATION")
    ucn_cpl = flopy.utils.HeadFile(cpl_ws / "gwt.ucn", text="CONCENTRATION")
    assert np.allclose(ucn.get_times(), ucn_cpl.get_times())
    assert np.allclose(ucn.get_alldata(), ucn_cpl.get_alldata())


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()
