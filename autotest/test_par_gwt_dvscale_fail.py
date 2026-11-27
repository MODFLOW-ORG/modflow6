"""
Test that the dependent_variable_scaling option must be specified for all GWT
models in the same IMS solution when the models are on different processes.
Reuses the simulation in test_gwt_dvscale.py, split into two models.

Cases:
  - par_dvsfail : the GWT model on one process scales and the GWT model on the
                  other does not (xfail).
"""

import numpy as np
import pytest
from flopy.mf6.utils import Mf6Splitter
from framework import TestFramework

cases = ["par_dvsfail"]
errmsg = "dependent_variable_scaling must be specified for all models"


def build_models(idx, test):
    from test_gwt_dvscale import build_model as build

    ws = test.workspace
    sim = build(idx, ws, dvscale=True)

    gwf = sim.get_model()
    nrow, ncol = gwf.dis.nrow.array, gwf.dis.ncol.array
    split_array = np.ones((nrow, ncol), dtype=int)
    split_array[:, int(ncol / 2) :] = 2
    new_sim = Mf6Splitter(sim).split_multi_model(split_array)
    new_sim.set_sim_path(ws)

    # turn off scaling for one of the split GWT models
    gwtnames = [name for name in new_sim.model_names if name.startswith("gwt")]
    assert len(gwtnames) == 2, f"expected two GWT models, found {gwtnames}"
    gwt = new_sim.get_model(gwtnames[-1])
    gwt.name_file.dependent_variable_scaling = False

    return new_sim, None


def check_output(test):
    found = any(errmsg in f.read_text() for f in test.workspace.glob("*.lst"))
    assert found, f"expected '{errmsg}' in a simulation listing file"


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=check_output,
        xfail=True,
        parallel=True,
        ncpus=2,
    )
    test.run()
