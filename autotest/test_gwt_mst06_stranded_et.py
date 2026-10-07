"""
Tests the MST stranded mass option for solute left behind by a sink that
carries no solute, here evapotranspiration given a mixed (AUXMIXED)
concentration of zero.

A thin surficial layer is drained by evapotranspiration whose extinction depth
reaches below the layer, and is rewetted by raising the head in the layer
below. Without the option the solute that evapotranspiration leaves behind is
concentrated into the water that remains, and the concentration of a nearly
dry cell grows without bound.

Cases:
  - standard : the share of that solute in the drained part of a cell is held,
               so the concentration stays bounded, the budget closes, and the
               held mass returns in the cells that resaturate.
  - newton   : the same with the Newton formulation, which keeps nearly dry
               cells active, so the held mass returns where the layer
               resaturates; with the standard formulation drained cells are
               inactive and keep what they hold.
  - nooption : with the Newton formulation and without the option, the
               transport solution fails (xfail).
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

ncol = 5
top, botm = 10.0, [8.0, 0.0]
porosity = 0.3
cinit = 100.0


def get_model(ws, newton, strand):
    sim = flopy.mf6.MFSimulation(sim_name="et", sim_ws=ws, exe_name="mf6")
    # evapotranspiration drains the surficial layer, then a higher head in the
    # layer below rewets it
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(100.0, 100, 1.0)] * 2)

    def add_ims(fname):
        return flopy.mf6.ModflowIms(
            sim,
            complexity="moderate",
            outer_dvclose=1e-8,
            inner_dvclose=1e-10,
            outer_maximum=200,
            linear_acceleration="bicgstab",
            filename=fname,
        )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="gwf",
        save_flows=True,
        newtonoptions="NEWTON" if newton else None,
    )
    sim.register_ims_package(add_ims("gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=2, nrow=1, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, k33=0.1)
    flopy.mf6.ModflowGwfic(gwf, strt=9.5)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=porosity, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={
            k: [[(1, 0, 0), h, cinit], [(1, 0, ncol - 1), h, cinit]]
            for k, h in enumerate((7.0, 10.5))
        },
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )
    # the extinction depth, 5 m, reaches below the 2 m surficial layer
    flopy.mf6.ModflowGwfevt(
        gwf,
        stress_period_data=[[(0, 0, j), top, 0.005, 5.0, 0.0] for j in range(ncol)],
        auxiliary="CONCENTRATION",
        pname="EVT-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    sim.register_ims_package(add_ims("gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=2, nrow=1, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    kwargs = {"porosity": porosity}
    if strand:
        kwargs.update(stranded_mass=True, stranded_filerecord=[("gwt.strand.bin",)])
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[
            ["CHD-1", "AUX", "CONCENTRATION"],
            ["EVT-1", "AUXMIXED", "CONCENTRATION"],
        ],
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def check_output(ws, newton):
    text = (ws / "gwt.lst").read_text()
    disc = re.findall(
        r"BUDGET FOR ENTIRE MODEL.*?PERCENT DISCREPANCY\s*=\s*([-+0-9.Ee]+)",
        text,
        re.S,
    )
    disc = np.abs(np.array(disc, dtype=float))
    assert disc.size > 0 and np.all(disc < 1e-2), f"budget does not close {disc}"

    ucn = flopy.utils.HeadFile(ws / "gwt.ucn", text="CONCENTRATION")
    conc = np.array([ucn.get_data(totim=t)[0, 0] for t in ucn.get_times()])
    active = conc > -1e20
    # the concentration of the surficial layer stays bounded; without the
    # option a nearly dry cell rises without limit
    assert conc[active].max() < 3.0 * cinit, f"peak {conc[active].max()}"

    sf = flopy.utils.HeadFile(ws / "gwt.strand.bin", text="STRANDED")
    strand = np.array([sf.get_data(totim=t)[0, 0] for t in sf.get_times()])
    hds = flopy.utils.HeadFile(ws / "gwf.hds")
    head = hds.get_data(totim=hds.get_times()[-1])[0, 0]
    sat = np.clip((head - botm[0]) / (top - botm[0]), 0.0, 1.0)
    assert strand[99].max() > 0.0, "no solute was held while the layer drained"
    # with the Newton formulation the higher head of the second period
    # resaturates the cells at the ends of the row, which return all of the
    # solute they held; evapotranspiration keeps the cells between them
    # drained, and they keep holding solute. With the standard formulation the
    # drained cells are inactive and, without rewetting, keep what they hold.
    wet = sat >= 1.0
    if newton:
        assert wet.any() and (~wet).any(), f"saturation at the end {sat}"
    assert np.allclose(strand[-1][wet], 0.0, atol=strand[99].max() * 1e-6), (
        f"held solute did not return {strand[-1][wet]}"
    )
    assert np.all(strand[-1][~wet] > 0.0), f"held solute {strand[-1][~wet]}"


@pytest.mark.parametrize("newton", [False, True], ids=["standard", "newton"])
def test_stranded_et(function_tmpdir, targets, newton):
    TestFramework(
        name="et",
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: get_model(t.workspace, newton, True),
        check=lambda t: check_output(t.workspace, newton),
    ).run()


def test_nooption(function_tmpdir, targets):
    def check(test):
        listing = (test.workspace / "mfsim.lst").read_text()
        assert "did not converge" in listing, "expected a convergence failure"

    TestFramework(
        name="et",
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: get_model(t.workspace, True, False),
        check=check,
        xfail=True,
    ).run()
