"""
Tests the warning for groups of connected cells with no storage,
dependent-variable boundary, or constant value, whose dependent variable is
not determined.

Each steady-state model is 1 layer by 3 rows by 10 columns with a constant
head in the first column. Setting IDOMAIN to 0 in column 5 cuts off columns 6
to 10. The first assembly is checked by default, the first outer iteration of
every time step with the -ci command-line option, every outer iteration with
-ci ALL, and the starting heads of every stress period with -m validate.

Cases:
  - iso_island      : the cut-off cells have no boundary, so one warning
                      counts them and a table lists them as cut off by the
                      input.
  - iso_ghb         : the cut-off cells have a GHB; no warning.
  - iso_later       : the cut-off cells have a GHB in stress period 1 only, so
                      the check of the first assembly finds nothing; no
                      warning.
  - iso_later_ci    : as iso_later with -ci, so the cells are found in stress
                      period 2.
  - iso_later_all   : as iso_later with -ci ALL; the cells are found in stress
                      period 2.
  - iso_later_val   : as iso_later with -m validate; the cells are found in
                      stress period 2 without solving.
  - iso_later_gwt   : as iso_later_val with a GWT model coupled in memory and a
                      GHB source in SSM; the GWT solution is not assembled.
  - iso_x_island    : two 5-column models joined by an exchange; IDOMAIN cuts
                      off the last two columns of the first model, which with
                      the second model have no boundary, so the warning
                      counts all 21 cells.
  - iso_x_connected : as iso_x_island without the cut, so the second model is
                      connected to the constant head through the exchange; no
                      warning.
"""

import re
from pathlib import Path

import flopy
import pytest
from framework import TestFramework

cases = [
    "iso_island",
    "iso_ghb",
    "iso_later",
    "iso_later_ci",
    "iso_later_all",
    "iso_later_val",
    "iso_later_gwt",
    "iso_x_island",
    "iso_x_connected",
]
nrow = 3
strt = 1.0
WARNING = "had no storage, dependent-variable boundary, or constant value"
VALIDATE = [["-m", "validate"]]
CARGS = {"ci": [["-ci"]], "all": [["-ci", "ALL"]], "val": VALIDATE, "gwt": VALIDATE}
SKIPPED = 'mode="validation" -- Skipping matrix assembly and solution.'
ASSEMBLED = 'mode="validation" -- Assembled at the starting heads'


def add_gwf(sim, name, ncol, idomain, chd=False, ghb=None, aux=False):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=1.0,
        botm=0.0,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    if chd:
        flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=[((0, i, 0), strt) for i in range(nrow)]
        )
    if ghb is not None:
        cells = [
            ((0, i, ncol - 1), strt, 1.0) + ((1.0,) if aux else ()) for i in range(nrow)
        ]
        flopy.mf6.ModflowGwfghb(
            gwf,
            pname="ghb",
            auxiliary=["concentration"] if aux else None,
            stress_period_data={kper: cells if on else [] for kper, on in ghb},
        )
    return gwf


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    nper = 2 if name.startswith("iso_later") else 1
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab")
    if name.startswith("iso_x"):
        ncol = 5
        idomain = [[1] * ncol for _ in range(nrow)]
        if name == "iso_x_island":
            for row in idomain:
                row[2] = 0
        left = add_gwf(sim, "left", ncol, [idomain], chd=True)
        right = add_gwf(sim, "right", ncol, 1)
        flopy.mf6.ModflowGwfgwf(
            sim,
            exgmnamea="left",
            exgmnameb="right",
            nexg=nrow,
            exchangedata=[
                ((0, i, ncol - 1), (0, i, 0), 1, 0.5, 0.5, 1.0) for i in range(nrow)
            ],
        )
        sim.register_ims_package(ims, [left.name, right.name])
    else:
        ncol = 10
        idomain = [[1] * ncol for _ in range(nrow)]
        for row in idomain:
            row[4] = 0
        ghb = None
        if name == "iso_ghb":
            ghb = [(0, True)]
        elif name.startswith("iso_later"):
            ghb = [(0, True), (1, False)]
        gwt = name.endswith("_gwt")
        add_gwf(sim, name, ncol, [idomain], chd=True, ghb=ghb, aux=gwt)
        if gwt:
            add_gwt(sim, name, ncol, [idomain])
    return sim, None


def add_gwt(sim, gwfname, ncol, idomain):
    name = f"t_{gwfname}"
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab", filename="gwt.ims")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name)
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=1,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=1.0,
        botm=0.0,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["ghb", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=name)
    sim.register_ims_package(ims, [name])


def warning_text(workspace):
    """Return the listing warnings with the line wrapping removed."""
    text = " ".join(
        Path(f).read_text() for f in sorted(Path(workspace).glob("mfsim*.lst"))
    )
    return re.sub(r"\s+", " ", text)


def check_output(idx, test):
    name = cases[idx]
    text = warning_text(test.workspace)
    if name in ("iso_ghb", "iso_later", "iso_x_connected"):
        assert WARNING not in text, f"unexpected warning for {name}"
        return
    # -- one warning for the simulation, and a table row for each cell
    assert text.count(WARNING) == 1, "expected one warning"
    model = name.upper()
    if name == "iso_island":
        assert f"15 cell(s) in solution SLN_1 (model(s) {model})" in text
        assert "15 not dry," in text
        assert f"{model} (1,1,6) INPUT NOT RECHECKED 1 1 1 1 1 1" in text
    elif name.startswith("iso_later"):
        # -- found in stress period 2, after the GHB is removed
        assert "15 not dry," in text
        assert f"{model} (1,1,6) INPUT NOT RECHECKED 2 1 2 1 1 1" in text
    if name == "iso_later_gwt":
        # -- the GWF solution is assembled and the GWT solution is not
        ws = Path(test.workspace)
        gwf_text = (ws / f"{name}.lst").read_text()
        gwt_text = (ws / f"t_{name}.lst").read_text()
        assert gwf_text.count(ASSEMBLED) == 2, "GWF not assembled in each period"
        assert ASSEMBLED not in gwt_text, "GWT solution was assembled"
        assert gwt_text.count(SKIPPED) == 2, "GWT not skipped in each period"
    elif name == "iso_x_island":
        assert "21 cell(s) in solution SLN_1 (model(s) LEFT, RIGHT)" in text
        assert "LEFT (1,1,4) INPUT NOT RECHECKED" in text
        assert "RIGHT (1,1,1) INPUT NOT RECHECKED" in text


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        cargs=CARGS.get(name.rsplit("_", 1)[-1]),
    )
    test.run()
