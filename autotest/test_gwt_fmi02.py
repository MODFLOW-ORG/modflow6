"""
Tests ability to run a transport model in an independent
simulation, consuming flows with a flow model interface.

Includes 4 cases, first with the transport model's TDIS
exactly matching the flow model's TDIS, then with more
time steps per transport model stress period (allowed),
then two more cases with fewer in the transport model,
which is forbidden.

The last test case is needed because the error trap is
implemented separately for the simulation's last period.
"""

import flopy
import pytest
from framework import TestFramework

simname = "gwtfmi02"
cases = [
    f"{simname}_a",  # matching TDIS
    f"{simname}_b",  # transport nstp > flow nstp = 1 (allowed)
    f"{simname}_c",  # transport nstp < flow nstp (not allowed)
    f"{simname}_d",  # like c, but only in final stress period
]


def get_model_name(name, mdl):
    return f"{name}_{mdl}"


def build_gwf_sim(name, ws, mf6, tdis_pd):
    gwf_name = get_model_name(name, "gwf")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)
    tdis = flopy.mf6.ModflowTdis(sim, nper=len(tdis_pd), perioddata=tdis_pd)
    ims = flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(gwf, nrow=10, ncol=10)
    ic = flopy.mf6.ModflowGwfic(gwf)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, save_saturation=True
    )
    spd = {
        0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
        1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf, pname="CHD-1", stress_period_data=spd, auxiliary=["concentration"]
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_name}.bud",
        head_filerecord=f"{gwf_name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_gwt_sim(name, gwf_ws, gwt_ws, mf6, tdis_pd):
    gwf_name = get_model_name(name, "gwf")
    gwt_name = get_model_name(name, "gwt")
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwt_ws, exe_name=mf6)
    tdis = flopy.mf6.ModflowTdis(sim, nper=len(tdis_pd), perioddata=tdis_pd)
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwt_name, save_flows=True)
    dis = flopy.mf6.ModflowGwtdis(gwt, nrow=10, ncol=10)
    ic = flopy.mf6.ModflowGwtic(gwt)
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.2)
    adv = flopy.mf6.ModflowGwtadv(gwt)
    pd = [("GWFBUDGET", gwf_ws / f"{gwf_name}.bud", None)]
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
    sources = [("CHD-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(gwt, print_flows=True, sources=sources)
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwt_name}.bud",
        concentration_filerecord=f"{gwt_name}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def build_models(idx, test):
    if "_a" in test.name:
        gwf_tdis_pd = [(1.0, 3, 1.0), (1.0, 3, 1.0)]
        gwt_tdis_pd = [(1.0, 3, 1.0), (1.0, 3, 1.0)]
    elif "_b" in test.name:
        gwf_tdis_pd = [(1.0, 1, 1.0), (1.0, 1, 1.0)]
        gwt_tdis_pd = [(1.0, 3, 1.0), (1.0, 3, 1.0)]
    elif "_c" in test.name:
        gwf_tdis_pd = [(1.0, 3, 1.0), (1.0, 3, 1.0)]
        gwt_tdis_pd = [(1.0, 1, 1.0), (1.0, 1, 1.0)]
    elif "_d" in test.name:
        gwf_tdis_pd = [(1.0, 3, 1.0), (1.0, 3, 1.0)]
        gwt_tdis_pd = [(1.0, 3, 1.0), (1.0, 1, 1.0)]
    gwf_sim = build_gwf_sim(
        test.name, test.workspace / "gwf", test.targets["mf6"], tdis_pd=gwf_tdis_pd
    )
    gwt_sim = build_gwt_sim(
        test.name,
        test.workspace / "gwf",
        test.workspace / "gwt",
        test.targets["mf6"],
        tdis_pd=gwt_tdis_pd,
    )
    return gwf_sim, gwt_sim


def check_output(idx, test):
    if "_a" in test.name:
        pass
    elif "_b" in test.name:
        pass
    elif "_c" in test.name:
        buff = test.buffs[1]
        assert any("INCOMPATIBLE WITH TIME" in line for line in buff)
    elif "_d" in test.name:
        buff = test.buffs[1]
        assert any("INCOMPATIBLE WITH TIME" in line for line in buff)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
        xfail=[False, "_c" in name or "_d" in name],
    )
    test.run()
