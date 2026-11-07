"""
Test that a horizontal flow barrier (HFB) stays in effect when hydraulic
conductivity changes during a simulation.

Steady flow between constant heads at the two ends of a row of cells,
with a barrier between the two middle cells. Three stress periods; period
2 has two time steps.

Case 1 (hfb-tvk01): TVK changes K at the cell next to the barrier in
    period 2; the HFB PERIOD block is given again in period 3. The heads
    from period 2 on must match a model with that K given in NPF.
Case 2 (hfb-tvk02): As case 1, with the Newton-Raphson formulation and
    the HFB PERIOD block given again in period 2, with the TVK change.
Case 3 (hfb-vsc03): VSC with a uniform temperature (one viscosity ratio
    for every cell, scaling K and the barrier alike), with the
    Newton-Raphson formulation. The heads must match the model without
    VSC.
Case 4 (hfb-vsc04): As case 3, without the Newton-Raphson formulation.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["hfb-tvk01", "hfb-tvk02", "hfb-vsc03", "hfb-vsc04"]
newton = [False, True, True, False]
tvk_case = [True, True, False, False]

nlay, nrow, ncol = 1, 1, 10
k0, k1 = 10.0, 20.0
hydchr = 1.0
cell_n, cell_m = (0, 0, 4), (0, 0, 5)
temperature = 35.0
perioddata = [(1.0, 1, 1.0), (1.0, 2, 1.0), (1.0, 1, 1.0)]
# periods (zero based) with an HFB PERIOD block
hfb_periods = [(0, 2), (0, 1), (0,), (0,)]


def build_sim(ws, name, exe, tvk=False, vsc=False):
    idx = cases.index(name)
    gwfname = f"gwf-{name}"
    gwtname = f"gwt-{name}"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=perioddata)

    # flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions=newton[idx])
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
    k = np.full((nlay, nrow, ncol), k0)
    if tvk_case[idx] and not tvk:
        # reference: the period 2 K given in NPF
        k[cell_n] = k1
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    if tvk:
        npf.tvk.initialize(
            filename=f"{gwfname}.tvk", perioddata={1: [(cell_n, "K", k1)]}
        )
    flopy.mf6.ModflowGwfic(gwf, strt=1.5)
    flopy.mf6.ModflowGwfhfb(
        gwf,
        maxhfb=1,
        stress_period_data={
            kper: [(cell_n, cell_m, hydchr)] for kper in hfb_periods[idx]
        },
        pname="HFB-1",
    )
    chd = [[(0, 0, 0), 2.0], [(0, 0, ncol - 1), 1.0]]
    if vsc:
        chd = [c + [temperature] for c in chd]
        flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=chd, auxiliary="temperature", pname="CHD-1"
        )
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
    else:
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd, pname="CHD-1")
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{gwfname}.hds", saverecord=[("HEAD", "ALL")]
    )

    if vsc:
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
        flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "TEMPERATURE")])
        flopy.mf6.ModflowGwfgwt(
            sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
        )
    sim.write_simulation(silent=True)
    return sim


def build_models(idx, test):
    name = cases[idx]
    tvk = tvk_case[idx]
    return build_sim(test.workspace, name, "mf6", tvk=tvk, vsc=not tvk), None


def check_output(idx, test):
    name = cases[idx]
    gwfname = f"gwf-{name}"

    # reference: period 2 K given in NPF (TVK cases), or no VSC (VSC cases)
    ws = test.workspace / "ref"
    sim = build_sim(ws, name, test.targets["mf6"])
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference simulation failed"

    h = flopy.utils.HeadFile(test.workspace / f"{gwfname}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{gwfname}.hds").get_alldata()
    # the barrier holds most of the head drop
    drop = href[-1].flat[cell_n[2]] - href[-1].flat[cell_m[2]]
    assert drop > 0.5, f"reference barrier head drop too small: {drop}"
    # time steps compared: from period 2 (TVK cases) or all (VSC cases)
    steps = range(1, len(h)) if tvk_case[idx] else range(len(h))
    for i in steps:
        assert np.allclose(h[i], href[i], atol=1e-8), (
            f"time step {i + 1} heads differ from reference: "
            f"max diff {np.abs(h[i] - href[i]).max()}"
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
