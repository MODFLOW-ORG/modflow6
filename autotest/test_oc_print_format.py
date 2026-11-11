"""
Check which forms of the PRINT_FORMAT record MODFLOW 6 accepts, in GWF's
Output Control and GWT's IST Package. The two read the record through
separate code.

The syntax is COLUMNS <columns> WIDTH <width> DIGITS <digits> [<format>].
The input loader accepts any leading part of the record, and COLUMNS 2 or
COLUMNS 2 WIDTH 16 used to run with the format silently ignored. They are
now an error. Malformed forms that the loader rejects are covered too, so
this test notices if that changes.
"""

import re
import subprocess

import flopy
import pytest

SYNTAX = "PRINT_FORMAT must include COLUMNS, WIDTH and DIGITS"

# (record, expected error or None if it should run)
cases = {
    "full": ("COLUMNS 2 WIDTH 16 DIGITS 9 EXPONENTIAL", None),
    "no-format": ("COLUMNS 2 WIDTH 16 DIGITS 9", None),
    "no-digits": ("COLUMNS 2 WIDTH 16", SYNTAX),
    "columns-only": ("COLUMNS 2", SYNTAX),
    "format-for-digits": ("COLUMNS 2 WIDTH 16 FIXED", 'tag "DIGITS"'),
    "no-width": ("COLUMNS 2 DIGITS 9", 'tag "WIDTH"'),
}


def build_sim(ws, exe):
    name = "pf"
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=exe, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    ims_gwf = flopy.mf6.ModflowIms(sim, filename="gwf.ims")
    sim.register_ims_package(ims_gwf, [gwf.name])
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=3, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 1.0], [(0, 0, 2), 0.5]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        headprintrecord=[("COLUMNS", 10, "WIDTH", 11, "DIGITS", 4, "GENERAL")],
        printrecord=[("HEAD", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    ims_gwt = flopy.mf6.ModflowIms(sim, filename="gwt.ims")
    sim.register_ims_package(ims_gwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, nlay=1, nrow=1, ncol=3, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwtic(gwt, strt=1.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.2)
    flopy.mf6.ModflowGwtist(
        gwt,
        cimprintrecord=[("COLUMNS", 10, "WIDTH", 11, "DIGITS", 4, "GENERAL")],
        porosity=0.1,
        volfrac=0.1,
        zetaim=0.01,
        cim=0.5,
    )
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtoc(gwt, printrecord=[("CONCENTRATION", "ALL")])
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea=gwf.name, exgmnameb=gwt.name)
    return sim


def set_print_format(path, record):
    """Replace the PRINT_FORMAT record flopy wrote, since flopy checks
    the record against its own DFNs and won't write the partial forms."""
    text = path.read_text()
    text, n = re.subn(r"PRINT_FORMAT .*", f"PRINT_FORMAT {record}", text)
    assert n == 1, f"no PRINT_FORMAT record in {path}"
    path.write_text(text)


@pytest.mark.parametrize("package", ["oc", "ist"])
@pytest.mark.parametrize("case", list(cases))
def test_oc_print_format(case, package, function_tmpdir, targets):
    record, error = cases[case]
    ws = function_tmpdir
    sim = build_sim(ws, str(targets["mf6"]))
    sim.write_simulation(silent=True)
    path = ws / ("gwf.oc" if package == "oc" else "gwt.ist")
    set_print_format(path, record)

    proc = subprocess.run([targets["mf6"]], cwd=ws, capture_output=True, text=True)
    out = " ".join((proc.stdout + proc.stderr).split())
    if error is None:
        assert "Normal termination" in out, out
        # the format is applied: two columns, 9 digits
        lst = (ws / ("gwf.lst" if package == "oc" else "gwt.lst")).read_text()
        table = "HEAD IN LAYER" if package == "oc" else "CIM IN LAYER"
        assert table in lst
        assert re.search(r"\d\.\d{9}", lst.split(table)[-1])
    else:
        assert "Normal termination" not in out, out
        assert error in out, out
