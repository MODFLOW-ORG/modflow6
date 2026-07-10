"""
Test that boundary packages with MAXBOUND=0 and no period data run successfully
with a warning, and that MAXBOUND=0 with period data still produces an error.

Cases:
  *_warn  - MAXBOUND=0, no period data: should succeed with a warning
  *_err   - MAXBOUND=0, with period data: should fail
"""

import pathlib as pl
import re

import flopy
import pytest
from framework import TestFramework

cases = [
    "bwel_warn",
    "bwel_err",
    "bghb_warn",
    "bghb_err",
    "bdrn_warn",
    "bdrn_err",
    "bchd_warn",
    "bchd_err",
]

xfail = [
    False,  # bwel_warn - MAXBOUND=0, no data: should succeed
    True,  # bwel_err  - MAXBOUND=0 with data: should fail
    False,  # bghb_warn - MAXBOUND=0, no data: should succeed
    True,  # bghb_err  - MAXBOUND=0 with data: should fail
    False,  # bdrn_warn - MAXBOUND=0, no data: should succeed
    True,  # bdrn_err  - MAXBOUND=0 with data: should fail
    False,  # bchd_warn - MAXBOUND=0, no data: should succeed
    True,  # bchd_err  - MAXBOUND=0 with data: should fail
]

# Simple 1-layer, 1-row, 3-column grid
nlay, nrow, ncol = 1, 1, 3
nper = 2
top = 10.0
botm = [0.0]
strt = 5.0
hk = 1.0


def _make_sim(name, ws):
    """Build base GWF simulation with DIS, IC, NPF, CHD, OC."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=str(ws),
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)],
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ims = flopy.mf6.ModflowIms(sim, complexity="simple")
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=10.0,
        delc=10.0,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk)
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("HEAD", "ALL")])
    return sim, gwf


def _patch_maxbound(ws, pkg_file):
    """Change MAXBOUND to 0 in the given package file."""
    path = pl.Path(ws) / pkg_file
    txt = path.read_text()
    txt = re.sub(r"(MAXBOUND\s+)\d+", r"\g<1>0", txt)
    path.write_text(txt)


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    is_err = "err" in name

    if "chd" in name:
        # CHD tests: don't add a base CHD (the package under test IS the CHD)
        sim, gwf = _make_sim(name, ws)
        if not is_err:
            flopy.mf6.ModflowGwfchd(gwf, maxbound=0)
        else:
            # Write with valid maxbound so flopy doesn't auto-compute to 0
            flopy.mf6.ModflowGwfchd(
                gwf,
                maxbound=1,
                stress_period_data={0: [((0, 0, 0), strt)]},
            )
    elif "wel" in name:
        sim, gwf = _make_sim(name, ws)
        # Add a CHD to keep head from going unbounded
        flopy.mf6.ModflowGwfchd(
            gwf, maxbound=1, stress_period_data={0: [((0, 0, 0), strt)]}
        )
        if not is_err:
            flopy.mf6.ModflowGwfwel(gwf, maxbound=0)
        else:
            flopy.mf6.ModflowGwfwel(
                gwf,
                maxbound=1,
                stress_period_data={0: [((0, 0, 1), -1.0)]},
            )
    elif "ghb" in name:
        sim, gwf = _make_sim(name, ws)
        flopy.mf6.ModflowGwfchd(
            gwf, maxbound=1, stress_period_data={0: [((0, 0, 0), strt)]}
        )
        if not is_err:
            flopy.mf6.ModflowGwfghb(gwf, maxbound=0)
        else:
            flopy.mf6.ModflowGwfghb(
                gwf,
                maxbound=1,
                stress_period_data={0: [((0, 0, 2), strt, 1.0)]},
            )
    elif "drn" in name:
        sim, gwf = _make_sim(name, ws)
        flopy.mf6.ModflowGwfchd(
            gwf, maxbound=1, stress_period_data={0: [((0, 0, 0), strt)]}
        )
        if not is_err:
            flopy.mf6.ModflowGwfdrn(gwf, maxbound=0)
        else:
            flopy.mf6.ModflowGwfdrn(
                gwf,
                maxbound=1,
                stress_period_data={0: [((0, 0, 2), strt - 1.0, 1.0)]},
            )

    if is_err:
        # Write files now so we can patch MAXBOUND before the framework re-writes
        sim.write_simulation()
        if "wel" in name:
            pkg_file = f"{name}.wel"
        elif "ghb" in name:
            pkg_file = f"{name}.ghb"
        elif "drn" in name:
            pkg_file = f"{name}.drn"
        else:  # chd
            pkg_file = f"{name}.chd"
        _patch_maxbound(ws, pkg_file)

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    if "warn" not in name:
        return

    # Warnings are written to mfsim.lst under "WARNING REPORT:"
    lst = pl.Path(test.workspace) / "mfsim.lst"
    lines = lst.read_text().splitlines()

    in_warning_report = False
    found_warning = False
    for line in lines:
        if "WARNING REPORT:" in line:
            in_warning_report = True
        if in_warning_report and "MAXBOUND of zero" in line:
            found_warning = True
            break

    assert found_warning, (
        f"Expected MAXBOUND=0 warning not found in mfsim.lst for {name}"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    is_err = "err" in name
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=xfail[idx],
        compare=None,
        overwrite=not is_err,  # don't overwrite pre-patched files for error cases
    )
    test.run()
