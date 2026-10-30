"""
Test that the upstream-weighted conductance of a GWF-GWF exchange does not
depend on whether the exchange or the interface model forms it.

The simulation has two solutions with the same pair of Newton models joined
by a GWF-GWF exchange; the second exchange uses the interface model. The cells
on either side of the exchange start at equal heads with different
saturations, so the upstream cell is undefined at the first outer iteration.
The outer iteration histories of the two solutions must be the same.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ifmod_upstream_tie"]

# spatial discretization data; the right model is thinner, so its saturation
# differs from the left model at the same head
nlay, nrow, ncol = 1, 1, 10
delr, delc = 100.0, 1.0
top = 100.0
botm = {"left": 0.0, "right": 40.0}

# hydraulic properties
hk = 10.0

# solver options
nouter, ninner = 100, 300
hclose, rclose = 1e-9, 1e-6

# starting head everywhere and at the left constant-head cell
strt = 50.0
rech = 0.01

# solution name, models in solution order, interface model flag
solutions = {
    "sln1": (["left1", "right1"], False),
    "sln2": (["left2", "right2"], True),
}


def add_model(sim, mname, side):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        xorigin=0.0 if side == "left" else ncol * delr,
        top=top,
        botm=botm[side],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk)
    if side == "left":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), strt]])
    else:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=rech)
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{mname}.hds", saverecord=[("HEAD", "LAST")]
    )


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    for sname, (mnames, ifmod) in solutions.items():
        ims = flopy.mf6.ModflowIms(
            sim,
            filename=f"{sname}.ims",
            print_option="SUMMARY",
            csv_outer_output_filerecord=f"{sname}.outer.csv",
            outer_dvclose=hclose,
            outer_maximum=nouter,
            under_relaxation="NONE",
            inner_maximum=ninner,
            inner_dvclose=hclose,
            rcloserecord=rclose,
            linear_acceleration="BICGSTAB",
        )
        add_model(sim, mnames[0], "left")
        add_model(sim, mnames[1], "right")
        sim.register_ims_package(ims, mnames)
        flopy.mf6.ModflowGwfgwf(
            sim,
            exgtype="GWF6-GWF6",
            exgmnamea=mnames[0],
            exgmnameb=mnames[1],
            filename=f"{sname}.gwfgwf",
            newton=True,
            nexg=1,
            exchangedata=[[(0, 0, ncol - 1), (0, 0, 0), 1, delr / 2, delr / 2, delc]],
            dev_interfacemodel_on=ifmod,
        )
    return sim, None


def check_output(idx, test):
    dvmax = [
        np.genfromtxt(
            os.path.join(test.workspace, f"{s}.outer.csv"), delimiter=",", names=True
        )["solution_outer_dvmax"]
        for s in solutions
    ]
    assert dvmax[0].size == dvmax[1].size, (
        f"outer iterations differ: {dvmax[0].size} and {dvmax[1].size}"
    )
    assert np.allclose(dvmax[0], dvmax[1], rtol=1e-6, atol=0.0), (
        f"outer iteration histories differ: {dvmax[0]} and {dvmax[1]}"
    )

    for a, b in (("left1", "left2"), ("right1", "right2")):
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
