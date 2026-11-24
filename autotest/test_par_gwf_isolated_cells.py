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
  - par_iso_x_solo      : as par_iso_x_connected with a third model in a
                          second solution on the first process only, so the
                          second process has no model in that solution; runs
                          without hanging and no warning.
  - par_iso_x_solo_val  : as par_iso_x_solo with -m validate.
"""

import flopy
import pytest
from framework import TestFramework
from test_gwf_isolated_cells import WARNING, add_gwf, warning_text
from test_gwf_isolated_cells import build_models as build_serial
from test_gwf_isolated_cells import cases as serial_cases

cases = [
    "par_iso_x_island",
    "par_iso_x_connected",
    "par_iso_x_island_ci",
    "par_iso_x_solo",
    "par_iso_x_solo_val",
]
CARGS = {"ci": [["-ci"]], "val": [["-m", "validate"]]}


def build_models(idx, test):
    name = cases[idx]
    solo = name.startswith("par_iso_x_solo")
    serial = "iso_x_connected" if solo else name.removeprefix("par_")
    sim, _ = build_serial(serial_cases.index(serial.removesuffix("_ci")), test)
    if solo:
        # -- a second solution whose only model is on the first process
        add_gwf(sim, "solo", 5, 1, chd=True)
        ims = flopy.mf6.ModflowIms(
            sim, linear_acceleration="bicgstab", filename="solo.ims"
        )
        sim.register_ims_package(ims, ["solo"])
        flopy.mf6.ModflowUtlhpc(
            sim, partitions=[("left", 0), ("right", 1), ("solo", 0)]
        )
    return sim, None


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
        cargs=CARGS.get(name.rsplit("_", 1)[-1]),
        parallel=True,
        ncpus=2,
    )
    test.run()
