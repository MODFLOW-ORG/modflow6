"""
Test that the memory manager's readonly/output variable attributes are
interpreted/enforced through the BMI.

- get_input_var_names()/get_output_var_names() are correctly filtered
- a readonly variable cannot be written via set_value()
"""

import flopy
import numpy as np
import pytest
from modflow_devtools.markers import requires_pkg

name = "bmiattrs"
nlay, nrow, ncol = 1, 1, 5
delr = delc = 1.0
top = 1.0
botm = [0.0]
k = 1.0
strt = 10.0
chd_head = 10.0
wellq = -1.0


@pytest.fixture
def simple_sim(tmp_path):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(tmp_path))
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-8, inner_dvclose=1e-8)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), chd_head]])
    flopy.mf6.ModflowGwfwel(
        gwf, pname="WEL_0", stress_period_data=[[(0, 0, ncol - 1), wellq]]
    )
    sim.write_simulation()
    return sim


@requires_pkg("xmipy")
def test_var_attributes(simple_sim, targets):
    from xmipy import XmiWrapper
    from xmipy.errors import XMIError

    sim = simple_sim
    mf6 = XmiWrapper(lib_path=targets["libmf6"], working_directory=sim.sim_path)
    mf6.initialize()

    simvals_addr = mf6.get_var_address("SIMVALS", name, "WEL_0")
    x_addr = mf6.get_var_address("X", name)

    input_vars = mf6.get_input_var_names()
    output_vars = mf6.get_output_var_names()

    # SIMVALS is readonly output
    assert simvals_addr not in input_vars
    assert simvals_addr in output_vars

    # X is writable (input) and output
    assert x_addr in input_vars
    assert x_addr in output_vars

    simvals = mf6.get_value_ptr(simvals_addr)
    with pytest.raises(XMIError, match="read-only"):
        mf6.set_value(simvals_addr, np.zeros_like(simvals))

    x = mf6.get_value_ptr(x_addr)
    new_x = np.full_like(x, 5.0)
    mf6.set_value(x_addr, new_x)
    np.testing.assert_array_equal(mf6.get_value_ptr(x_addr), new_x)

    mf6.finalize()
