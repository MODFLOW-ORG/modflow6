"""
flopy4-run version of test_par_gwf01: instead of MODFLOW 6's native parallel
run path, flopy4 re-writes the simulation as ASCII input and runs it itself,
under MPI, via its own `Simulation.run(ncpus=...)`. The flopy4 run is checked
against the same expected head values as the base test.
"""

import pytest
from framework import TestFramework, require_flopy4, write_fp4_ascii
from test_par_gwf01 import build_models, cases, check_output
from test_par_gwf_splitter01 import build_models as splitter_build_models
from test_par_gwf_splitter01 import cases as splitter_cases
from test_par_gwf_splitter01 import check_output as splitter_check_output


@pytest.mark.netcdf
@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model_fp4(idx, name, function_tmpdir, targets):
    require_flopy4()
    from flopy4.mf6 import Simulation

    def check(test):
        check_output(idx, test)

        workspace = test.workspace / "fp4_parallel"
        write_fp4_ascii(test.workspace, workspace)
        sim = Simulation.load(workspace / "mfsim.nam")
        sim.run(exe=test.targets["mf6"], ncpus=2)

        test.workspace = workspace
        check_output(idx, test)

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=check,
        compare=None,
    )
    test.run()


@pytest.mark.netcdf
@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(splitter_cases))
def test_mf6model_fp4_splitter(idx, name, function_tmpdir, targets):
    """Same as test_mf6model_fp4, but for a real Mf6Splitter/HPC-partition
    topology (4 domains) instead of a hand-built 2-model exchange -- also
    confirms flopy4 round-trips the HPC partition file (sim.hpc)."""
    require_flopy4()
    from flopy4.mf6 import Simulation

    def check(test):
        splitter_check_output(idx, test)

        workspace = test.workspace / "fp4_parallel"
        write_fp4_ascii(test.workspace, workspace)
        sim = Simulation.load(workspace / "mfsim.nam")
        sim.run(exe=test.targets["mf6"], ncpus=4)

        test.workspace = workspace
        splitter_check_output(idx, test)

    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: splitter_build_models(idx, t),
        check=check,
        compare=None,
    )
    test.run()
