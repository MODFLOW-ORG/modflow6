"""
Validates the variable-modifiability classes documented in the mf6api guide
(doc/mf6api/) by driving a model through the MODFLOW 6 API and modifying
memory-manager variables at run time.

The guide claims, for a GWF model under the standard (non-XT3D) formulation:
  - NPF/CONDSAT is a "Derived variable": setting it changes intercell
    conductance, so results change. It is computed once during Allocate and Read
    and then used directly, so it can be set a single time and persists.
  - NPF/K11 is an "Inert input": it is only used to compute CONDSAT during
    Allocate and Read, so setting it afterward has NO effect on results.
  - STO/SS and STO/SY are "Live input": they are read every formulate step, so
    setting them (once) changes results and persists.
  - Boundary stress variables such as GHB/COND are "Period-reset": they are
    re-read from the input context at the start of each stress period, so a
    value set once before the time loop is overwritten and has no effect, while
    a value set every time step (after prepare_time_step) takes effect.

The model is a one-dimensional convertible (water-table) row that drains toward
a constant head, with a general-head boundary, so conductance, storage, and the
boundary all measurably affect the transient heads.
"""

import numpy as np
from modflow_devtools.markers import requires_pkg

NAME = "apimod"
MODELNAME = NAME.upper()

# --- grid / properties ---
nlay, nrow, ncol = 1, 1, 10
delr = delc = 10.0
top = 10.0
botm = 0.0
strt = 10.0
chd_head = 1.0  # low head at the right end -> row drains
ghb_col = 4
ghb_head = 5.0
ghb_cond = 1.5
k11 = 1.0
ss = 1e-4
sy = 0.2

# --- time: drain over 10 steps ---
nper, nstp, perlen = 1, 10, 10.0


def _build_model(ws):
    import flopy

    sim = flopy.mf6.MFSimulation(
        sim_name=NAME, version="mf6", exe_name="mf6", sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(perlen, nstp, 1.0)]
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-9,
        outer_maximum=100,
        inner_dvclose=1e-10,
        inner_maximum=100,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=NAME, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=[botm],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k11, save_flows=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=ss, sy=sy, transient={0: True})
    flopy.mf6.ModflowGwfghb(
        gwf,
        pname="ghb",
        stress_period_data={0: [[(0, 0, ghb_col), ghb_head, ghb_cond]]},
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, ncol - 1), chd_head]]})
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{NAME}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    sim.write_simulation(silent=True)


def _run_api(libmf6, ws, address=None, factor=1.0, when="init", callback=None):
    """Run the model through the API and return the final head array.

    If address (a tuple passed to get_var_address) is given, the variable is
    scaled by factor either once after initialize() (when="init", i.e. after
    Allocate and Read) or every time step after prepare_time_step (when="step").
    If callback is given, it is called with the ModflowApi object every time step
    after prepare_solve (used here to relocate a boundary).
    """
    from modflowapi import ModflowApi

    mf6 = ModflowApi(str(libmf6), working_directory=str(ws))
    mf6.initialize()

    def _apply():
        ptr = mf6.get_value_ptr(mf6.get_var_address(*address))
        ptr[:] = ptr[:] * factor

    if address is not None and when == "init":
        _apply()

    max_iter = int(mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))[0])
    head_tag = mf6.get_var_address("X", MODELNAME)

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)
        mf6.prepare_solve()
        if address is not None and when == "step":
            _apply()
        if callback is not None:
            callback(mf6)
        kiter = 0
        while kiter < max_iter:
            if mf6.solve():
                break
            kiter += 1
        else:
            mf6.finalize()
            raise RuntimeError("API model failed to converge")
        mf6.finalize_solve()
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

    head = mf6.get_value(head_tag).copy()
    mf6.finalize()
    return head


@requires_pkg("modflowapi")
def test_api_variable_modifiability(function_tmpdir, targets):
    libmf6 = targets["libmf6"]
    ws = function_tmpdir
    _build_model(ws)

    baseline = _run_api(libmf6, ws)

    # NPF/CONDSAT -- Derived variable: setting it once must change the heads.
    h_condsat = _run_api(libmf6, ws, ("CONDSAT", MODELNAME, "NPF"), 0.5)
    assert not np.allclose(baseline, h_condsat, atol=1e-4), (
        "Setting NPF/CONDSAT had no effect, but it is documented as a Derived "
        "variable that controls intercell conductance."
    )

    # NPF/K11 -- Inert input: setting it after Allocate and Read must NOT change
    # the heads (CONDSAT was already computed from it).
    h_k11 = _run_api(libmf6, ws, ("K11", MODELNAME, "NPF"), 0.5)
    assert np.allclose(baseline, h_k11, atol=1e-9), (
        "Setting NPF/K11 changed the results, but it is documented as an Inert "
        "input under the standard formulation (set CONDSAT instead)."
    )

    # STO/SY and STO/SS -- Live input: setting them once must change the heads.
    h_sy = _run_api(libmf6, ws, ("SY", MODELNAME, "STO"), 0.5)
    assert not np.allclose(baseline, h_sy, atol=1e-4), (
        "Setting STO/SY had no effect, but it is documented as a Live input."
    )
    h_ss = _run_api(libmf6, ws, ("SS", MODELNAME, "STO"), 100.0)
    assert not np.allclose(baseline, h_ss, atol=1e-4), (
        "Setting STO/SS had no effect, but it is documented as a Live input."
    )

    # GHB/COND -- Period-reset: a value set once before the loop is overwritten
    # at the first stress period's Read and Prepare (no effect), but a value set
    # every time step takes effect.
    h_cond_once = _run_api(libmf6, ws, ("COND", MODELNAME, "GHB"), 10.0, "init")
    assert np.allclose(baseline, h_cond_once, atol=1e-9), (
        "Setting GHB/COND once before the time loop changed the results, but it "
        "is documented as Period-reset (overwritten at Read and Prepare)."
    )
    h_cond_step = _run_api(libmf6, ws, ("COND", MODELNAME, "GHB"), 10.0, "step")
    assert not np.allclose(baseline, h_cond_step, atol=1e-4), (
        "Setting GHB/COND every time step had no effect, but it is documented "
        "as a Period-reset variable that can be set after prepare_time_step."
    )

    # GHB/NODELIST -- a boundary can be relocated by setting its (one-based,
    # reduced) node number every time step. Move the GHB from its original cell
    # (col ghb_col) to col 0, which must change the heads.
    def relocate(mf6):
        nodelist = mf6.get_value_ptr(mf6.get_var_address("NODELIST", MODELNAME, "GHB"))
        nodelist[0] = 1  # one-based reduced node number of (0, 0, 0)

    h_moved = _run_api(libmf6, ws, callback=relocate)
    assert not np.allclose(baseline, h_moved, atol=1e-4), (
        "Relocating the GHB by setting NODELIST had no effect, but boundary "
        "locations are documented as modifiable through the API."
    )


HFB_NAME = "apihfb"


def _build_hfb_model(ws):
    import flopy

    sim = flopy.mf6.MFSimulation(
        sim_name=HFB_NAME, version="mf6", exe_name="mf6", sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(10.0, 5, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=HFB_NAME, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=10.0, delc=10.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-4, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), 9.0], [(0, 0, 4), 1.0]]}
    )
    # a strong barrier (low hydraulic characteristic) between col 1 and col 2
    flopy.mf6.ModflowGwfhfb(
        gwf,
        pname="hfb",
        maxhfb=1,
        stress_period_data=[[(0, 0, 1), (0, 0, 2), 0.01]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{HFB_NAME}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation(silent=True)


def _run_hfb(libmf6, ws, modify=None):
    """Run the HFB model, optionally modifying HYDCHR or CONDSAT each step."""
    from modflowapi import ModflowApi

    name = HFB_NAME.upper()
    mf6 = ModflowApi(str(libmf6), working_directory=str(ws))
    mf6.initialize()
    max_iter = int(mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))[0])

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)
        mf6.prepare_solve()
        if modify == "hydchr":
            mf6.get_value_ptr(mf6.get_var_address("HYDCHR", name, "HFB"))[:] = 100.0
        elif modify == "condsat":
            mf6.get_value_ptr(mf6.get_var_address("CONDSAT", name, "NPF"))[:] = 1e6
        kiter = 0
        while kiter < max_iter:
            if mf6.solve():
                break
            kiter += 1
        else:
            mf6.finalize()
            raise RuntimeError("HFB API model failed to converge")
        mf6.finalize_solve()
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

    head = mf6.get_value(mf6.get_var_address("X", name)).copy()
    mf6.finalize()
    return head


@requires_pkg("modflowapi")
def test_hfb_api_modifiability(function_tmpdir, targets):
    libmf6 = targets["libmf6"]
    ws = function_tmpdir
    _build_hfb_model(ws)

    baseline = _run_hfb(libmf6, ws)

    # HFB/HYDCHR -- Inert: HFB re-reads it and bakes the barrier into CONDSAT at
    # Read and Prepare, so setting it through the API has no effect (non-XT3D).
    h_hydchr = _run_hfb(libmf6, ws, "hydchr")
    assert np.allclose(baseline, h_hydchr, atol=1e-9), (
        "Setting HFB/HYDCHR changed the results, but it is documented as inert "
        "to the API (the barrier is applied to CONDSAT at Read and Prepare)."
    )

    # The barrier is changed by modifying NPF/CONDSAT, which does take effect.
    h_condsat = _run_hfb(libmf6, ws, "condsat")
    assert not np.allclose(baseline, h_condsat, atol=1e-4), (
        "Modifying NPF/CONDSAT did not change the barrier, but it is documented "
        "as the run-time knob for the HFB Package."
    )
