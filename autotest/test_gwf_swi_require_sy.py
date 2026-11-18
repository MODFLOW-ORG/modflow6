"""
Validation test for the seawater intrusion (SWI) package input.

The SWI Package uses specific yield (SY) as the drainable porosity for the
interface-movement storage.  A moving interface changes both the freshwater
and saltwater storage even where the aquifer is confined (iconvert=0), so SY
is required whenever the STO Package is active with SWI.  STO zero-fills SY
when it is not provided and does not require it for confined cells, so without
this check a confined transient SWI run would silently compute zero interface
storage.

This test builds a confined, transient, single-fluid SWI model whose STO
Package supplies a zero specific yield (equivalent, internally, to a
hand-written STO input that omits SY, since STO zero-fills an unspecified SY)
and confirms that the run terminates with the expected error message (the
model is expected to fail; xfail=True). Note that the flopy STO interface
defaults SY to a nonzero value when the argument is omitted, so SY is set to
zero explicitly here to reach the all-zero code path.
"""

import flopy
import pytest
from framework import TestFramework

cases = ["swi-nosy"]


def build_models(idx, test):
    Lx = 10000.0  # meters
    delr, delc = 100.0, 1.0
    ncol = int(Lx / delr)
    nlay, nrow = 1, 1
    top, botm = 0.0, -400.0
    k = 10.0
    recharge = 0.001
    h0 = 0.0

    ws = test.workspace
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")

    # transient so that storage (and thus SY) is required
    perioddata = [(100000.0, 10, 1.0)]
    flopy.mf6.ModflowTdis(sim, perioddata=perioddata)
    flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab")

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)  # confined
    # STO is active and transient, but SY is not meaningfully specified (zero)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1.0e-5, sy=0.0)
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    cghb = 1.0 * delr * delc / 10.0
    flopy.mf6.ModflowGwfghb(
        gwf, stress_period_data=[[0, 0, 0, h0, cghb], [0, 0, ncol - 1, h0, cghb]]
    )
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    return sim, None


def check_output(idx, test):
    with open(test.workspace / "mfsim.lst", "r") as f:
        lines = f.readlines()
    error_count = sum(1 for line in lines if "specific yield (SY) be specified" in line)
    assert error_count == 1, (
        "expected exactly one SWI 'SY required' error in mfsim.lst, "
        f"found {error_count}"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=True,
    )
    test.run()
