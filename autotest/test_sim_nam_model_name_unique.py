"""
Test to confirm that MODFLOW 6 does not explicitly enforce that model
names (mname, doc/mf6io/mf6ivar/dfn/sim-nam.dfn) are unique across the
"models" block of the simulation name file, even though uniqueness is an
implicit, load-bearing requirement.

check_model_name() (src/SimulationCreate.f90) only validates that a model
name is <=16 characters and has no embedded spaces -- there is no
duplicate-name check. But model names are used to build memory-manager
paths (e.g. this%memoryPath = create_mem_path(modelname) in
src/Model/GroundWaterFlow/gwf.f90) and are resolved via ifind()
(src/Utilities/ArrayHandlers.f90), which returns only the first matching
element. So giving two models the same name isn't rejected outright, but
it corrupts the memory manager (colliding memory addresses trigger
"Already existing variable" warnings, see PtrHashTableType%add in
src/Utilities/PtrHashTable.f90) and breaks solution assignment (the
second, identically-named model is never found by ifind(), so it's
silently left without a solution).

flopy itself refuses to let two models share a name -- MFSimulation
registers models in a dict keyed by name, so the second one just
overwrites the first internally. To get a genuinely duplicate-named
"models" block past flopy, two distinctly-named models are built and
written normally, and the generated mfsim.nam is then hand-edited to
register the second model under the first model's name (same trick used
in autotest/test_gwf_errors.py to exercise error paths flopy won't
generate on its own).
"""

import subprocess

import flopy


def build_two_gwf_models(ws, exe):
    sim = flopy.mf6.MFSimulation(
        sim_name="dupname", exe_name=exe, version="mf6", sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

    def add_model(name, chd_left, chd_right):
        gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
        flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=5)
        flopy.mf6.ModflowGwfic(gwf, strt=0.0)
        flopy.mf6.ModflowGwfnpf(gwf, save_flows=True)
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: [[(0, 0, 0), chd_left], [(0, 0, 4), chd_right]]},
        )
        flopy.mf6.ModflowGwfoc(
            gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
        )
        ims = flopy.mf6.ModflowIms(sim, filename=f"{name}.ims")
        sim.register_solution_package(ims, [gwf.name])
        return gwf

    add_model("model1", 1.0, 0.0)
    add_model("model2", 10.0, 0.0)
    sim.write_simulation()
    return sim


def duplicate_model2_name(ws):
    """Register model2 under model1's name in mfsim.nam.

    Leaves model2's own input files (model2.nam, model2.ims, ...) alone --
    only the name it's registered under in the simulation name file's
    "models" and "solutiongroup" blocks is changed, so the only thing
    being tested is what happens when two models share an mname.
    """
    namfile = ws / "mfsim.nam"
    text = namfile.read_text()
    text = text.replace("model2.nam  model2", "model2.nam  model1")
    text = text.replace("model2.ims  model2", "model2.ims  model1")
    namfile.write_text(text)


def run_mf6(argv, ws):
    proc = subprocess.Popen(
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, _ = proc.communicate()
    buff = result.decode("utf-8").splitlines() if result is not None else []
    return proc.returncode, buff


def test_model_name_uniqueness_is_implicit(function_tmpdir, targets):
    mf6 = targets["mf6"]

    build_two_gwf_models(function_tmpdir, mf6)
    duplicate_model2_name(function_tmpdir)

    returncode, buff = run_mf6([mf6], str(function_tmpdir))
    text = "\n".join(buff)

    assert returncode != 0, "mf6 unexpectedly succeeded with two models sharing a name"

    # there's no dedicated "duplicate model name" check, so mf6 doesn't
    # fail with a clean, purpose-built error message
    assert "duplicate" not in text.lower()

    # every memory-manager entry for the second model instead collides
    # with the first model's identically-named entries
    assert "Already existing variable being added to the HashTable" in text
    assert "MODEL1" in text

    # and because name-based lookups (ifind) resolve to the first match
    # only, the second (duplicate-named) model is silently dropped from
    # its solution, which is reported as a distinct, unrelated-looking
    # error
    assert "Model was not assigned to a solution: MODEL1" in text
