"""
Test that the pseudo-transient continuation (PTC) factor does not depend on
the order of the models in a solution.

The simulation has two solutions with the same pair of unconnected
steady-state Newton models listed in opposite order. One model starts at its
solution (zero residual) and the other does not, so each model must use its
own residuals for the initial PTC time step to be the same in both solutions.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ptc02"]

# spatial discretization data
nlay, nrow, ncol = 1, 1, 100
delr, delc = 50.0, 1.0
top = 25.0
botm = 0.0

# hydraulic properties
hk = 50.0

# solver options
nouter, ninner = 500, 300
hclose, rclose = 1e-9, 1e-6

# constant heads at both ends; the near model starts at this head with no
# recharge, so its initial residual is zero
hchd = 20.0
rech = 0.001

# solution name, models in solution order
solutions = {
    "sln1": ["near1", "far1"],
    "sln2": ["far2", "near2"],
}


def add_model(sim, mname, near):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=hchd if near else botm)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk)
    if not near:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=rech)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), hchd], [(0, 0, ncol - 1), hchd]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{mname}.hds", saverecord=[("HEAD", "LAST")]
    )


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    for sname, mnames in solutions.items():
        ims = flopy.mf6.ModflowIms(
            sim,
            filename=f"{sname}.ims",
            print_option="SUMMARY",
            dev_ptc_output_filerecord=f"{sname}.ptc",
            outer_dvclose=hclose,
            outer_maximum=nouter,
            under_relaxation="NONE",
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=rclose,
            linear_acceleration="BICGSTAB",
        )
        for mname in mnames:
            add_model(sim, mname, near=mname.startswith("near"))
        sim.register_ims_package(ims, mnames)
    return sim, None


def read_ptcdel(path):
    """Return the PTC time step of the first outer iteration."""
    with open(path) as f:
        for line in f:
            tokens = line.split()
            if tokens and tokens[0] == "1":
                return float(tokens[1])
    raise ValueError(f"no first outer iteration in {path}")


def check_output(idx, test):
    ptcdel = [read_ptcdel(os.path.join(test.workspace, f"{s}.ptc")) for s in solutions]
    assert np.isclose(ptcdel[0], ptcdel[1], rtol=1e-12, atol=0.0), (
        f"initial PTC time step depends on model order: {ptcdel}"
    )

    for a, b in (("near1", "near2"), ("far1", "far2")):
        ha = flopy.utils.HeadFile(os.path.join(test.workspace, f"{a}.hds")).get_data()
        hb = flopy.utils.HeadFile(os.path.join(test.workspace, f"{b}.hds")).get_data()
        assert np.allclose(ha, hb, atol=10 * hclose), f"heads differ for {a} and {b}"


@pytest.mark.developmode
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
