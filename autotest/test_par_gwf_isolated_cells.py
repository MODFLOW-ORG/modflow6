"""
Tests the warning for groups of cells with an undetermined dependent variable
in a parallel simulation, where a group can cross processes.

The two-model cases of test_gwf_isolated_cells.py run with one model on each
of two processes.

Cases:
  - par_iso_x_island    : the group of 21 cells crosses the processes and has
                          no boundary; each process warns about its own cells.
  - par_iso_x_connected : the second model is connected to the constant head
                          only through the first, on the other process; no
                          warning.
  - par_iso_x_island_ci : as par_iso_x_island with -ci, so every time step is
                          checked; each process still warns once.
"""

import pytest
from framework import TestFramework
from test_gwf_isolated_cells import WARNING, warning_text
from test_gwf_isolated_cells import build_models as build_serial
from test_gwf_isolated_cells import cases as serial_cases

cases = ["par_iso_x_island", "par_iso_x_connected", "par_iso_x_island_ci"]


def build_models(idx, test):
    name = cases[idx].removeprefix("par_").removesuffix("_ci")
    return build_serial(serial_cases.index(name), test)


def check_output(idx, test):
    name = cases[idx]
    text = warning_text(test.workspace)
    if name.startswith("par_iso_x_island"):
        # -- each process warns once about its own cells
        assert text.count(WARNING) == 2, "each process should warn once"
        assert "6 cell(s) in solution SLN_1 (model(s) LEFT)" in text
        assert "15 cell(s) in solution SLN_1 (model(s) RIGHT)" in text
        assert "LEFT (1,1,4) INPUT" in text
        assert "RIGHT (1,1,1) INPUT" in text
    else:
        assert WARNING not in text, f"unexpected warning for {name}"


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        cargs=[["-ci"]] if name.endswith("_ci") else None,
        parallel=True,
        ncpus=2,
    )
    test.run()
