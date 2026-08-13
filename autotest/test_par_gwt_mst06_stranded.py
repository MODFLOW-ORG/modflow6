"""
This test reuses the simulation data and config in
test_gwt_mst06_stranded_gwtgwt.py and runs it in parallel mode.

The stranded mass reservoirs belong to a cell rather than to a connection, so
nothing about them crosses the exchange between the two models. This confirms
that a domain which strands and returns mass gives the same answer whether its
models are solved on one process or on two.
"""

import pytest
from framework import TestFramework

cases = ["par_mst06_stranded"]


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    from test_gwt_mst06_stranded_gwtgwt import build_models, check_output

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=2,
    )
    test.run()
