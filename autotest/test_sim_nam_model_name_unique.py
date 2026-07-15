"""
Test that MODFLOW 6 rejects a simulation whose name file "models" block
(doc/mf6io/mf6ivar/dfn/sim-nam.dfn) assigns the same mname to more than
one model.

check_model_name() (src/SimulationCreate.f90) validates mname against the
model names already seen earlier in the models block, in addition to its
existing length (<=16 characters) and no-embedded-spaces checks, and
fails with an informative "Model name is not unique" error.

That check runs in models_create(), which is called after the IDM
subsystem has already loaded each model's own input files into the
memory manager, keyed in part by mname (src/Utilities/Idm/IdmLoad.f90,
load_models()). So a duplicate mname still causes a batch of "Already
existing variable being added to the HashTable" memory-manager warnings
(src/Utilities/PtrHashTable.f90) before the simulation aborts with the
clear, purpose-built error.

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


def test_duplicate_model_name_is_rejected(function_tmpdir, targets):
    mf6 = targets["mf6"]

    build_two_gwf_models(function_tmpdir, mf6)
    duplicate_model2_name(function_tmpdir)

    returncode, buff = run_mf6([mf6], str(function_tmpdir))
    text = "\n".join(buff)

    assert returncode != 0, "mf6 unexpectedly succeeded with two models sharing a name"

    # memory-manager entries for the second model collide with the first
    # model's identically-named entries, logged as warnings, before the
    # simulation aborts
    assert "Already existing variable being added to the HashTable" in text
    assert "MODEL1" in text

    # check_model_name() catches the duplicate and reports a clear,
    # purpose-built error rather than failing later with something more
    # cryptic (e.g. a solution-assignment error)
    assert "Invalid model name: MODEL1" in text
    assert "Model name is not unique" in text
