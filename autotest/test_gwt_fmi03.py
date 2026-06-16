"""
Tests that multiple GWT models in the same simulation can reference
the same GWF head and budget files via the FMI package without error.

Previously this caused an error because the second FMI open attempt
found the output files already open on units held by the first FMI.

Reported in https://github.com/MODFLOW-ORG/modflow6/issues/2832.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

simname = "gwtfmi03"
cases = [simname]

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


def build_gwt_model(sim, name, gwf_ws, gwf_name, model_suffix):
    gwt_name = get_model_name(name, model_suffix)
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwt_name, save_flows=True)
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
    pd = [
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
    return gwt


def build_gwt_sim(name, gwf_ws, gwt_ws, mf6):
    gwf_name = get_model_name(name, "gwf")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwt_ws, exe_name=mf6)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_pd)
    ims = flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="BICGSTAB",
        filename="gwt.ims",
    )
    gwt1 = build_gwt_model(sim, name, gwf_ws, gwf_name, "gwt1")
    gwt2 = build_gwt_model(sim, name, gwf_ws, gwf_name, "gwt2")
    sim.register_ims_package(ims, [gwt1.name, gwt2.name])
    return sim


def build_models(idx, test):
    gwf_ws = test.workspace / "gwf"
    gwt_ws = test.workspace / "gwt"
    gwf_sim = build_gwf_sim(test.name, gwf_ws, test.targets["mf6"])
    gwt_sim = build_gwt_sim(test.name, gwf_ws, gwt_ws, test.targets["mf6"])
    return gwf_sim, gwt_sim


def check_output(idx, test):
    gwt_ws = test.workspace / "gwt"
    name = test.name

    def read_conc(gwt_name):
        fpth = gwt_ws / f"{gwt_name}.ucn"
        return flopy.utils.HeadFile(
            fpth, precision="double", text="CONCENTRATION"
        ).get_alldata()

    conc1 = read_conc(get_model_name(name, "gwt1"))
    conc2 = read_conc(get_model_name(name, "gwt2"))

    assert np.allclose(conc1, conc2)
    assert np.all(conc1 >= 0.0)
    assert np.all(conc1 <= 1.0 + 1.0e-6)


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
