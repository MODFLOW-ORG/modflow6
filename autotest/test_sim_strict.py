"""Test STRICT mode, which raises errors on unknown input tags."""

import subprocess

import flopy


def run_mf6(exe, ws):
    proc = subprocess.Popen(
        [exe], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, _ = proc.communicate()
    buff = result.decode("utf-8").splitlines() if result else []
    return proc.returncode, buff


def build_sim(ws, exe, strict=False):
    name = "simstrict"
    simkwargs = {"strict": True} if strict else {}
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=exe, sim_ws=ws, **simkwargs
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, top=0.0, botm=[-1.0])
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), 0.0], [(0, 2, 2), 1.0]]}
    )
    sim.write_simulation()
    return name


def write_bad_option(ws, name):
    npf_file = ws / f"{name}.npf"
    with open(npf_file, "w") as f:
        f.write("BEGIN options\n")
        f.write("  UNKNOWN\n")
        f.write("END options\n\n")
        f.write("BEGIN griddata\n")
        f.write("  ICELLTYPE\n")
        f.write("    CONSTANT 0\n")
        f.write("  K\n")
        f.write("    CONSTANT 1.0\n")
        f.write("END griddata\n")


def test_unknown_option(function_tmpdir, targets):
    """By default, an unrecognized variable is a warning."""
    mf6 = targets["mf6"]
    name = build_sim(function_tmpdir, mf6, strict=False)
    write_bad_option(function_tmpdir, name)

    returncode, buff = run_mf6(mf6, function_tmpdir)
    assert returncode == 0, "mf6 failed unexpectedly:\n" + "\n".join(buff)

    lst = (function_tmpdir / "mfsim.lst").read_text()
    assert "UNKNOWN" in lst
    assert "ignored" in lst.lower()
    assert "Normal termination" in lst


def test_unknown_option_strict(function_tmpdir, targets):
    """With STRICT, an unrecognized variable is an error."""
    mf6 = targets["mf6"]
    name = build_sim(function_tmpdir, mf6, strict=True)
    write_bad_option(function_tmpdir, name)

    returncode, _ = run_mf6(mf6, function_tmpdir)
    assert returncode != 0
