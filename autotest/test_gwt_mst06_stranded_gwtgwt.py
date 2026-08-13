"""
Tests the MST stranded mass option for a domain split into two models coupled
by a GWT-GWT exchange, which is the arrangement a parallel run uses.

Four cells in a row are divided between a left and a right model. A constant
head in the right model drains and rewets the whole row, so cells of both
models strand mass and return it. The reservoirs belong to a cell rather than
to a connection, so nothing crosses the exchange, and the checks below hold
whether the two models are solved on one process or on two.

Cases:
  - mst06_gwtgwt : the budget of each model closes, both models strand mass,
                   and both reservoirs empty when the row is saturated again.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mst06_gwtgwt"]

top, botm = 40.0, 0.0
delr = delc = 100.0
porosity = 0.3
sy = 0.15
rhob = 1600.0
kd = 1.0e-4
cinit = 100.0
ncol = 2  # columns in each model

# the leading period does not drain; the row is fully saturated again at the end
heads = [40.0, 30.0, 20.0, 10.0, 20.0, 30.0] + [40.0] * 6
nper = len(heads)
perioddata = [(10.0, 1, 1.0)] * nper


def add_solver(sim, filename):
    return flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=filename,
    )


def add_gwf(sim, name, chd):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=1,
        k=100.0,
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=0.0, sy=sy, transient={0: True}, save_flows=True
    )
    if chd:
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={
                k: [[(0, 0, ncol - 1), h, 0.0]] for k, h in enumerate(heads)
            },
            auxiliary="CONCENTRATION",
            pname="CHD-1",
            save_flows=True,
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def add_gwt(sim, name, chd):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        save_flows=True,
        sorption="linear",
        bulk_density=rhob,
        distcoef=kd,
        stranded_mass=True,
        stranded_filerecord=[(f"{name}.strand.bin",)],
    )
    if chd:
        flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{name}.ucn",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    return gwt


def build_models(idx, test):
    sim = flopy.mf6.MFSimulation(
        sim_name=cases[idx], version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=perioddata)

    for name, chd in (("fl", False), ("fr", True)):
        add_gwf(sim, name, chd)
    # the models of an exchange share a solution
    sim.register_ims_package(add_solver(sim, "f.ims"), ["fl", "fr"])

    # the last column of the left model abuts the first of the right
    exchangedata = [
        ((0, 0, ncol - 1), (0, 0, 0), 1, delr / 2.0, delr / 2.0, delc, 0.0, delr)
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=1,
        exchangedata=exchangedata,
        exgmnamea="fl",
        exgmnameb="fr",
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="flfr.gwfgwf",
    )

    for name, chd in (("tl", False), ("tr", True)):
        add_gwt(sim, name, chd)
    sim.register_ims_package(add_solver(sim, "t.ims"), ["tl", "tr"])

    flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1="fl",
        gwfmodelname2="fr",
        nexg=1,
        exchangedata=exchangedata,
        exgmnamea="tl",
        exgmnameb="tr",
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="tltr.gwtgwt",
    )
    for f, t in (("fl", "tl"), ("fr", "tr")):
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=f,
            exgmnameb=t,
            filename=f"{f}{t}.gwfgwt",
        )
    return sim, None


def discrepancies(ws, name):
    """Cumulative percent discrepancies of the budget of one model."""
    out = []
    model_budget = False
    with open(ws / f"{name}.lst") as f:
        for line in f:
            if "BUDGET FOR ENTIRE MODEL" in line:
                model_budget = True
            elif "PERCENT DISCREPANCY" in line and model_budget:
                out.append(float(line.replace("PERCENT DISCREPANCY =", " ").split()[0]))
                model_budget = False
    return np.array(out)


def check_output(idx, test):
    ws = test.workspace
    for name in ("tl", "tr"):
        # the first period does not drain, so its discrepancy is a ratio of
        # two numbers near zero and carries no information
        disc = discrepancies(ws, name)[1:]
        assert len(disc) > 0, f"no model budget was written for {name}"
        assert np.allclose(disc, 0.0, atol=1e-2), (
            f"the budget of {name} does not close, {disc[np.abs(disc) > 1e-2]}"
        )

        f = flopy.utils.HeadFile(ws / f"{name}.strand.bin", text="STRANDED")
        st = np.array([f.get_data(totim=t).sum() for t in f.get_times()])
        assert st.max() > 0.0, f"{name} stranded no mass"
        # the row is saturated again at the end, so the reservoirs are empty
        assert np.isclose(st[-1], 0.0, atol=st.max() * 1e-3), (
            f"the reservoirs of {name} did not empty, {st[-1]} of {st.max()} remains"
        )


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
