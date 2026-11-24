"""
Test head-dependent boundary conductance with the viscosity (VSC)
package.

Steady flow between a constant head at one end of a row of cells and
head-dependent boundaries at the other, with a uniform temperature. VSC
then scales hydraulic conductivity and the boundary conductance by the
same viscosity ratio, so the heads must match the model without VSC and
every flow is the flow without VSC times that ratio. Two stress periods
with two time steps each; period 1 has a boundary at the last cell,
period 2 a boundary at the next to last cell listed before it.

Case 1 (vsc08-ghb): GHB.
Case 2 (vsc08-ghbm): GHB with an AUXMULTCOL conductance multiplier.
Case 3 (vsc08-drnm): DRN with an AUXMULTCOL conductance multiplier.
Case 4 (vsc08-rivm): RIV with an AUXMULTCOL conductance multiplier.
Case 5 (vsc08-ghbts): GHB with conductance from a time series that
    changes at the second time step.
Case 6 (vsc08-ghbmts): GHB with an AUXMULTCOL multiplier from a time
    series that changes at the second time step.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "vsc08-ghb",
    "vsc08-ghbm",
    "vsc08-drnm",
    "vsc08-rivm",
    "vsc08-ghbts",
    "vsc08-ghbmts",
]
package = ["ghb", "ghb", "drn", "riv", "ghb", "ghb"]
multiplier = [None, 0.5, 0.5, 0.5, None, "multts"]
timeseries = [None, None, None, None, "condts", "multts"]

nlay, nrow, ncol = 1, 1, 10
k = 10.0
cond0, cond1 = 5.0, 50.0
mult0, mult1 = 0.5, 0.25
temperature = 35.0
bcells = {0: [(0, 0, ncol - 1)], 1: [(0, 0, ncol - 2), (0, 0, ncol - 1)]}


def build_sim(ws, name, exe, vsc):
    idx = cases.index(name)
    gwfname = f"gwf-{name}"
    gwtname = f"gwt-{name}"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(2.0, 2, 1.0)] * 2)

    # flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    ims = flopy.mf6.ModflowIms(
        sim,
        inner_dvclose=1e-10,
        outer_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=1.0, delc=1.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    flopy.mf6.ModflowGwfic(gwf, strt=1.5)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 2.0, temperature]],
        auxiliary="temperature",
        save_flows=True,
        pname="CHD-1",
    )

    # head-dependent boundaries at the other end
    aux = [temperature]
    auxnames = ["temperature"]
    kwargs = {}
    if multiplier[idx] is not None:
        aux.append(multiplier[idx])
        auxnames.append("mult")
        kwargs["auxmultname"] = "mult"
    cond = "condts" if timeseries[idx] == "condts" else cond0
    pkg = package[idx]
    rivbot = [0.5] if pkg == "riv" else []
    spd = {
        kper: [[cell, 1.0, cond, *rivbot, *aux] for cell in cells]
        for kper, cells in bcells.items()
    }
    pkgtype = {
        "ghb": flopy.mf6.ModflowGwfghb,
        "drn": flopy.mf6.ModflowGwfdrn,
        "riv": flopy.mf6.ModflowGwfriv,
    }[pkg]
    bnd = pkgtype(
        gwf,
        stress_period_data=spd,
        auxiliary=auxnames,
        save_flows=True,
        pname="BND-1",
        **kwargs,
    )
    if timeseries[idx] is not None:
        v0, v1 = (cond0, cond1) if timeseries[idx] == "condts" else (mult0, mult1)
        bnd.ts.initialize(
            filename=f"{gwfname}.{pkg}.ts",
            timeseries=[(0.0, v0), (1.0, v1), (4.0, v1)],
            time_series_namerecord=timeseries[idx],
            interpolation_methodrecord="stepwise",
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    if vsc:
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=8.904e-4,
            thermal_formulation="nonlinear",
            thermal_a2=10.0,
            thermal_a3=248.37,
            thermal_a4=133.16,
            nviscspecies=1,
            packagedata=[(0, 0.0, 20.0, gwtname, "temperature")],
            filename=f"{gwfname}.vsc",
        )

        # transport model holding a uniform temperature
        gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
        imsgwt = flopy.mf6.ModflowIms(
            sim,
            inner_dvclose=1e-10,
            outer_dvclose=1e-10,
            linear_acceleration="BICGSTAB",
            filename=f"{gwtname}.ims",
        )
        sim.register_ims_package(imsgwt, [gwtname])
        flopy.mf6.ModflowGwtdis(
            gwt,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=1.0,
            delc=1.0,
            top=10.0,
            botm=0.0,
        )
        flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
        flopy.mf6.ModflowGwtic(gwt, strt=temperature)
        flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
        flopy.mf6.ModflowGwtssm(
            gwt,
            sources=[
                ("CHD-1", "AUX", "TEMPERATURE"),
                ("BND-1", "AUX", "TEMPERATURE"),
            ],
        )
        flopy.mf6.ModflowGwfgwt(
            sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
        )
    sim.write_simulation(silent=True)
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, cases[idx], "mf6", vsc=True), None


def check_output(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"

    # reference: the model without VSC
    ws = test.workspace / "ref"
    sim = build_sim(ws, name, test.targets["mf6"], vsc=False)
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference simulation failed"

    h = flopy.utils.HeadFile(test.workspace / f"{gwfname}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{gwfname}.hds").get_alldata()
    for i in range(len(h)):
        assert np.allclose(h[i], href[i], atol=1e-8), (
            f"time step {i + 1} heads differ from the model without VSC: "
            f"max diff {np.abs(h[i] - href[i]).max()}"
        )

    # flows: the boundary flows are the flows without VSC times the
    # viscosity ratio, measured from the constant head flow
    text = package[idx].upper()
    cbc = flopy.utils.CellBudgetFile(
        test.workspace / f"{gwfname}.cbc", precision="double"
    )
    cbcref = flopy.utils.CellBudgetFile(ws / f"{gwfname}.cbc", precision="double")
    chd = cbc.get_data(text="CHD")
    chdref = cbcref.get_data(text="CHD")
    q = cbc.get_data(text=text)
    qref = cbcref.get_data(text=text)
    assert len(q) == len(qref) == len(h), f"{text} budget records missing"
    for i in range(len(h)):
        ratio = chd[i]["q"].sum() / chdref[i]["q"].sum()
        assert abs(ratio - 1.0) > 1e-3, "viscosity ratio is one"
        assert np.allclose(q[i]["q"], ratio * qref[i]["q"], atol=1e-8), (
            f"time step {i + 1} {text} flows differ from the flows without VSC "
            f"times the viscosity ratio {ratio}: "
            f"max diff {np.abs(q[i]['q'] - ratio * qref[i]['q']).max()}"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
