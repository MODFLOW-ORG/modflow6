"""
Test time-varying hydraulic conductivity (TVK) with the viscosity (VSC)
package. VSC scales hydraulic conductivity by the viscosity ratio each
time step; a TVK change must be scaled, not replaced by the scaled
original value.

Steady flow between constant heads, with a uniform temperature (one
viscosity ratio for every cell). TVK changes hydraulic conductivity at
the middle cells in period 2, which has three time steps. The period 2
heads and flows must match a model with that conductivity given in NPF.

Case 1 (vsc07-tvk01): one layer, flow along a row; TVK changes K.
Case 2 (vsc07-tvk02): two layers, flow from layer 1 to layer 2; TVK
    changes K33, given as a ratio of K (K33OVERK).
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["vsc07-tvk01", "vsc07-tvk02"]

nlay = [1, 2]
nrow, ncol = 1, 10
k0, k1 = 10.0, 100.0
k33ratio0, k33ratio1 = 0.1, 2.0
middle = range(3, 7)
temperature = 35.0


def build_sim(ws, name, exe, tvk):
    idx = cases.index(name)
    gwfname = f"gwf-{name}"
    gwtname = f"gwt-{name}"
    nl = nlay[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 3, 1.0)])

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
    botm = [-float(k) for k in range(nl)]
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nl, nrow=nrow, ncol=ncol, delr=1.0, delc=1.0, top=1.0, botm=botm
    )

    # period 2 values given in NPF (reference) or by TVK
    k = np.full((nl, nrow, ncol), k0)
    k33 = np.full((nl, nrow, ncol), k33ratio0)
    if idx == 0:
        cells = [(0, 0, j) for j in middle]
        changes = [(cell, "K", k1) for cell in cells]
    else:
        cells = [(kl, 0, j) for kl in range(nl) for j in middle]
        changes = [(cell, "K33", k33ratio1) for cell in cells]
    if not tvk:
        for cell in cells:
            if idx == 0:
                k[cell] = k1
            else:
                k33[cell] = k33ratio1
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=k, k33=k33, k33overk=True, save_flows=True
    )
    if tvk:
        npf.tvk.initialize(filename=f"{gwfname}.tvk", perioddata={1: changes})
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
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
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), 2.0, temperature],
            [(nl - 1, 0, ncol - 1), 1.0, temperature],
        ],
        auxiliary="temperature",
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
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
        gwt, nlay=nl, nrow=nrow, ncol=ncol, delr=1.0, delc=1.0, top=1.0, botm=botm
    )
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtic(gwt, strt=temperature)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "TEMPERATURE")])
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    sim.write_simulation(silent=True)
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, cases[idx], "mf6", tvk=True), None


def check_output(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"

    # reference: the period 2 values given in NPF; period 1 differs, so
    # only period 2 (time steps 2 to 4) is compared
    ws = test.workspace / "ref"
    sim = build_sim(ws, name, test.targets["mf6"], tvk=False)
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference simulation failed"

    h = flopy.utils.HeadFile(test.workspace / f"{gwfname}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{gwfname}.hds").get_alldata()
    assert np.abs(h[1] - h[0]).max() > 1e-3, "TVK change has no effect on heads"
    assert np.allclose(h[1:], href[1:], atol=1e-8), (
        "period 2 heads differ from NPF reference: "
        f"max diff {np.abs(h[1:] - href[1:]).max()}"
    )

    cbc = flopy.utils.CellBudgetFile(
        test.workspace / f"{gwfname}.cbc", precision="double"
    )
    cbcref = flopy.utils.CellBudgetFile(ws / f"{gwfname}.cbc", precision="double")
    q = cbc.get_data(text="FLOW-JA-FACE")
    qref = cbcref.get_data(text="FLOW-JA-FACE")
    assert len(q) == len(qref) == 4, "FLOW-JA-FACE records missing"
    for i in range(1, 4):
        assert np.allclose(q[i], qref[i], atol=1e-8), (
            f"time step {i + 1} flows differ from NPF reference: "
            f"max diff {np.abs(q[i] - qref[i]).max()}"
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
