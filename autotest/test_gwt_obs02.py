"""
Tests that continuous observations with three-digit exponents are written to
a comma-separated value file with the E character.

Two cells that exchange no water keep their initial concentrations of 1e-120
and 1e130.

Cases:
  - digits : DIGITS 12, which wrote 0.100000000000-119 rather than
             0.100000000000E-119 before the fix.
  - default: DIGITS not specified.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["obs02_digits", "obs02_default"]
digits = [12, None]
cinit = [1e-120, 1e130]


def build_models(idx, test):
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim)

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    sim.register_ims_package(flopy.mf6.ModflowIms(sim, filename="gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 1.0], [(0, 0, 1), 1.0]]
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    sim.register_ims_package(flopy.mf6.ModflowIms(sim, filename="gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(gwt, nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowUtlobs(
        gwt,
        digits=digits[idx],
        continuous={
            "conc.csv": [
                ("small", "concentration", (0, 0, 0)),
                ("large", "concentration", (0, 0, 1)),
            ]
        },
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def check_output(idx, test):
    line = (test.workspace / "conc.csv").read_text().splitlines()[1]
    # every value has to be read as a number, so none may drop the E
    conc = np.array([float(v) for v in line.split(",")[1:]])
    assert np.allclose(conc, cinit, rtol=1e-10, atol=0.0), f"{line}"


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
