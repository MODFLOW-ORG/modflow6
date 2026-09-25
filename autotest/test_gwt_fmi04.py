"""
Tests that a GWT model can read the GWF model's binary grid
file via the FMI package's GWFGRID entry.

Previously GWT treated GWFGRID as the name of an advanced
package and tried to read the grid file as a budget file,
which failed with an invalid method code error.

Reported in https://github.com/MODFLOW-ORG/modflow6/issues/3009.

The second case gives the GWT model an IDOMAIN that differs
from the GWF model's and expects the grid check to fail.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

simname = "gwtfmi04"
cases = [simname, f"{simname}idm"]

nlay, nrow, ncol = 1, 1, 10
delr = delc = 1.0
top = 1.0
botm = [0.0]
nper = 2
tdis_pd = [(5.0, 5, 1.0)] * nper
porosity = 0.1


def get_model_name(name, mdl):
    return f"{name}_{mdl}"


def build_gwf_sim(name, ws, mf6):
    gwf_name = get_model_name(name, "gwf")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_pd)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)
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
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_saturation=True,
    )
    chd_spd = {0: [[(0, 0, 0), 1.0, 1.0], [(0, 0, ncol - 1), 0.0, 0.0]]}
    flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=chd_spd,
        auxiliary=["concentration"],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_name}.bud",
        head_filerecord=f"{gwf_name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_gwt_sim(name, gwf_ws, gwt_ws, mf6):
    gwf_name = get_model_name(name, "gwf")
    gwt_name = get_model_name(name, "gwt")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwt_ws, exe_name=mf6)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_pd)
    flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwt_name, save_flows=True)
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    if "idm" in name:
        idomain[0, 0, ncol // 2] = 0
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
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    pd = [
        ("GWFGRID", gwf_ws / f"{gwf_name}.dis.grb", None),
        ("GWFBUDGET", gwf_ws / f"{gwf_name}.bud", None),
        ("GWFHEAD", gwf_ws / f"{gwf_name}.hds", None),
    ]
    flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwt_name}.bud",
        concentration_filerecord=f"{gwt_name}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    gwf_ws = test.workspace / "gwf"
    gwt_ws = test.workspace / "gwt"
    gwf_sim = build_gwf_sim(test.name, gwf_ws, test.targets["mf6"])
    gwt_sim = build_gwt_sim(test.name, gwf_ws, gwt_ws, test.targets["mf6"])
    return gwf_sim, gwt_sim


def check_output(idx, test):
    if "idm" in test.name:
        buff = test.buffs[1]
        assert any("do not have the same discretization" in l for l in buff)
        return

    gwt_name = get_model_name(test.name, "gwt")
    fpth = test.workspace / "gwt" / f"{gwt_name}.ucn"
    conc = flopy.utils.HeadFile(
        fpth, precision="double", text="CONCENTRATION"
    ).get_alldata()

    # solute should have entered the domain and stay bounded
    assert conc[-1].max() > 0.0
    assert np.all(conc >= 0.0)
    assert np.all(conc <= 1.0 + 1.0e-6)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        xfail=[False, "idm" in name],
    )
    test.run()
