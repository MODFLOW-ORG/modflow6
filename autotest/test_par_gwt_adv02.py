"""
This test reuses the simulation data and config in
test_gwt_adv01_gwtgwt.py and runs it in parallel mode.
"""

import pytest
from framework import TestFramework

cases = ["par_adv01a_gwtgwt", "par_adv01b_gwtgwt", "par_adv01c_gwtgwt"]


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    from test_gwt_adv02_gwtgwt import build_models, check_output, number_of_models

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=number_of_models,
    )
    test.run()
