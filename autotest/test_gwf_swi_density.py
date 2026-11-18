"""
Test the RELATIVE_DENSITY_DIFFERENCE input option of the SWI Package.

The interface elevation follows from the two heads as

    zeta = -alphaf * hf + alphas * hs,    alphaf = 1 / delta,  alphas = alphaf + 1

where delta is the relative density difference, (rhos - rhof) / rhof. With the
saltwater head at its default of zero this reduces to zeta = -hf / delta, which
is satisfied cell by cell regardless of the flow solution and so makes a direct check
that the option reached the formulation.

Cases:
  a  default (option omitted)          -> delta = 0.025, alphaf = 40
  b  delta specified explicitly as 0.025 -> must reproduce case a exactly
  c  delta = 0.05                      -> alphaf = 20, a different solution
  d  delta = 0.0                       -> input error (expected failure)
  e  two-fluid, the two models disagree -> exchange error (expected failure)
"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "swi-density-a",
    "swi-density-b",
    "swi-density-c",
    "swi-density-d",
    "swi-density-e",
]

# None means the option is omitted entirely
delta = [None, 0.025, 0.05, 0.0, None]
XFAIL = {3, 4}

Lx = 10000.0
delr, delc = 100.0, 1.0
ncol = int(Lx / delr) + 1
nlay, nrow = 1, 1
top, botm = 0.0, -400.0
recharge = 0.001
k = 10.0
hbnd = 0.0


def build_mismatched_two_fluid(test):
    """Two coupled models that disagree on the relative density difference.

    The exchange must reject this: the interface elevation is computed from
    both heads using coefficients derived from delta, so two values would place
    the interface differently in each model with no outward sign of trouble.
    """
    sim = flopy.mf6.MFSimulation(
        sim_name="mymodel", sim_ws=test.workspace, exe_name="mf6"
    )
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, print_option="summary", linear_acceleration="bicgstab")
    for is_salt, d in ((False, 0.025), (True, 0.030)):
        nm = "saltwater" if is_salt else "freshwater"
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=nm, save_flows=True, newtonoptions="newton"
        )
        flopy.mf6.ModflowGwfdis(
            gwf, nlay=1, nrow=1, ncol=5, delr=delr, delc=delc, top=top, botm=botm
        )
        flopy.mf6.ModflowGwfic(gwf, strt=0.0)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
        flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=0.2, transient=True)
        flopy.mf6.ModflowGwfswi(gwf, relative_density_difference=d)
        flopy.mf6.ModflowGwfoc(
            gwf, head_filerecord=nm + ".hds", saverecord=[("HEAD", "ALL")]
        )
    flopy.mf6.ModflowSwiswi(
        sim, exgtype="SWI6-SWI6", exgmnamea="freshwater", exgmnameb="saltwater"
    )
    return sim, None


def build_models(idx, test):
    if idx == 4:
        return build_mismatched_two_fluid(test)
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=500,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, alternative_cell_averaging="amt-hmk", icelltype=0, k=k)
    swi_kwargs = dict(zeta_filerecord=name + ".zta")
    if delta[idx] is not None:
        swi_kwargs["relative_density_difference"] = delta[idx]
    flopy.mf6.ModflowGwfswi(gwf, **swi_kwargs)
    cghb = 10.0 * delr * delc / 10.0
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=[[0, 0, 0, hbnd, cghb], [0, 0, ncol - 1, hbnd, cghb]],
    )
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim, None


def _head_and_zeta(test):
    ws = pl.Path(test.workspace)
    gwf = test.sims[0].gwf[0]
    head = gwf.output.head().get_data().flatten()
    zeta = flopy.utils.HeadFile(ws / "mymodel.zta", text="zeta").get_data().flatten()
    return head, zeta


def check_output(idx, test):
    if idx in XFAIL:
        return

    head, zeta = _head_and_zeta(test)
    d = 0.025 if delta[idx] is None else delta[idx]

    # zeta = -hf / delta, cell by cell, independent of the flow solution
    expected = -head / d
    assert np.allclose(zeta, expected, rtol=1.0e-9, atol=1.0e-9), (
        f"zeta does not equal -head/delta for delta={d}: "
        f"max |diff| = {np.abs(zeta - expected).max():.3e}"
    )
    print(
        f"delta={d}  alphaf={1.0 / d:g}  max|zeta - (-h/delta)| = "
        f"{np.abs(zeta - expected).max():.3e}"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=idx in XFAIL,
    )
    test.run()
