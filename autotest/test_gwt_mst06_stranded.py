"""
Tests the MST stranded mass option, which holds solute out of the mobile domain
when a cell drains and returns it when the cell rewets.

A single thick cell is drained and rewetted by a constant head in the cell next
to it, so the saturation of the tested cell varies over a wide range.

Cases:
  - mst06_noop  : stranding active without sorption has nothing to strand and
                  must reproduce the results of a model with the option off.
  - mst06_sorb  : sorbed mass only; the mass held out of the mobile domain over
                  a drainage step is checked against the analytical amount.
  - mst06_cycle : the cell is drained and rewetted over a full cycle, so all of
                  the stranded mass must return and the reservoir must empty.
  - mst06_osc   : repeated drain and rewet cycles must not drift.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mst06_noop", "mst06_sorb", "mst06_cycle", "mst06_osc"]

top, botm = 40.0, 0.0
perlen = 10.0
delr = delc = 100.0
area = delr * delc
vcell = area * (top - botm)
porosity = 0.3
rhob = 1600.0
kd = 1.0e-4
cinit = 100.0

# head in the constant head cell for each stress period
# the leading period does not drain, which lets the saturation of the
# previous time step initialize before any mass is stranded
drain = [40.0, 35.0, 30.0, 25.0, 20.0, 15.0]
# the last periods let the cell equilibrate back to full saturation
cycle = drain + [20.0, 25.0, 30.0, 35.0] + [40.0] * 6
# identical drain and rewet cycles, each ending fully saturated
osc = ([30.0, 20.0, 30.0] + [40.0] * 5) * 4

heads = {
    "mst06_noop": drain,
    "mst06_sorb": drain,
    "mst06_cycle": cycle,
    "mst06_osc": osc,
}


def get_model(name, ws, hds, strand, sorb):
    nper = len(hds)
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(perlen, 1, 1.0)] * nper
    )

    gwfname, gwtname = "gwf-" + name, "gwt-" + name
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwfname}.ims",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, newtonoptions="NEWTON"
    )
    sim.register_ims_package(imsgwf, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=True, save_saturation=True, icelltype=1, k=100.0
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=porosity, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={k: [[(0, 0, 1), h, 0.0]] for k, h in enumerate(hds)},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=2, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)

    kwargs = {"porosity": porosity, "save_flows": True}
    if sorb:
        kwargs.update(sorption="linear", bulk_density=rhob, distcoef=kd)
    if strand:
        kwargs.update(stranded_mass=True)
        if sorb:
            kwargs.update(stranded_filerecord=[(f"{gwtname}.strand.bin",)])
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)

    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sorb = name != "mst06_noop"
    sim = get_model(name, test.workspace, heads[name], True, sorb)
    mc = None
    if name == "mst06_noop":
        # the same model with the stranded mass option off
        mc = get_model(
            name, os.path.join(test.workspace, "mf6"), heads[name], False, sorb
        )
    return sim, mc


def discrepancies(ws, name):
    """Cumulative percent discrepancies of the model budget.

    Only the cumulative value of the budget of the whole model is used. The rate
    discrepancy of a time step in which the flows have gone to zero, and the
    discrepancy of the budget of the stranded mass before any mass is stranded,
    are ratios of numbers near machine precision and carry no information.
    """
    fpth = os.path.join(ws, f"gwt-{name}.lst")
    assert os.path.isfile(fpth)
    out = []
    model_budget = False
    with open(fpth) as f:
        for line in f:
            if "BUDGET FOR ENTIRE MODEL" in line:
                model_budget = True
            elif "PERCENT DISCREPANCY" in line and model_budget:
                out.append(float(line.replace("PERCENT DISCREPANCY =", " ").split()[0]))
                model_budget = False
    return np.array(out)


def series(ws, name, ext, text):
    """Values for the tested cell over time from a binary output file."""
    fpth = os.path.join(ws, f"gwt-{name}.{ext}")
    assert os.path.isfile(fpth)
    f = flopy.utils.HeadFile(fpth, text=text)
    return np.array([f.get_data(totim=t).flatten()[0] for t in f.get_times()])


def check_output(idx, test):
    name = cases[idx]
    ws = test.workspace

    # the model budget must close in every case; decay of stranded mass never
    # passes through the concentration equation, so an entry there would show up
    # the discrepancy of the first period is a ratio of two numbers near zero,
    # because that period does not drain the cell, and carries no information
    disc = discrepancies(ws, name)[1:]
    assert np.allclose(disc, 0.0, atol=1e-2), (
        f"model budget does not close, discrepancies {disc[np.abs(disc) > 1e-2]}"
    )

    if name == "mst06_noop":
        # the specific yield equals the porosity, so no water is retained
        # against drainage, and without sorption there is nothing to strand
        base = series(os.path.join(ws, "mf6"), name, "ucn", "CONCENTRATION")
        conc = series(ws, name, "ucn", "CONCENTRATION")
        assert np.allclose(conc, base, rtol=1e-6), (
            f"stranding with nothing to strand changed the result: {conc} vs {base}"
        )
        return

    stranded = series(ws, name, "strand.bin", "STRANDED")

    if name == "mst06_sorb":
        # the mass held out of the mobile domain over a drainage step is the
        # sorbed mass of the part of the cell that drained, which follows the
        # change in the saturation of the cell
        conc = series(ws, name, "ucn", "CONCENTRATION")
        hds = flopy.utils.HeadFile(os.path.join(ws, f"gwf-{name}.hds"))
        h = np.array([hds.get_data(totim=t).flatten()[0] for t in hds.get_times()])
        sat = np.clip((h - botm) / (top - botm), 0.0, 1.0)
        sat_prev = np.concatenate(([1.0], sat[:-1]))
        expected = np.cumsum((sat_prev - sat) * vcell * rhob * kd * conc)
        assert np.allclose(stranded, expected, rtol=1e-3), (
            f"stranded mass {stranded} does not match the sorbed mass of the "
            f"part of the cell that drained {expected}"
        )

    elif name == "mst06_cycle":
        # the cell drains and then rewets to its original saturation, so all of
        # the mass that was held out of the mobile domain must return
        assert stranded.max() > 0.0, "no mass was stranded"
        assert np.isclose(stranded[-1], 0.0, atol=stranded.max() * 1e-4), (
            f"the reservoir did not empty on rewetting, {stranded[-1]} remains "
            f"of a maximum of {stranded.max()}"
        )

    elif name == "mst06_osc":
        # the reservoir must empty at the end of every cycle, so that repeated
        # drainage and rewetting cannot ratchet mass out of the mobile domain.
        # The amount stranded within a cycle grows as the cell concentrates,
        # which is a change in concentration rather than a loss of mass; the
        # budget check above is what tests conservation.
        assert stranded.max() > 0.0, "no mass was stranded"
        ends = stranded[7::8]
        assert np.allclose(ends, 0.0, atol=stranded.max() * 1e-4), (
            f"the reservoir did not empty at the end of every cycle, {ends}"
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
