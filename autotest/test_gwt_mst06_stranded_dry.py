"""
Tests the MST stranded mass option for a cell that drains completely, is
deactivated for transport, and is later rewet.

The saturation of an inactive cell is recorded, so that the step in which it
becomes active again is seen as rewetting rather than as further drainage.

Cases:
  - nodecay : the reservoir holds its mass unchanged while the cell is dry and
              returns it when the cell is rewet.
  - decay   : the same, except that the reservoir decays while the cell is dry,
              at the rate given for its phase.

Each case runs with the standard formulation, in which the drained cell is
inactive for transport and rewet by the flow model, and with the Newton
formulation, in which it stays active at a saturation of zero.
"""

import flopy
import numpy as np
import pytest

nlay, nrow, ncol = 3, 1, 2
delr = delc = 100.0
top = 30.0
botm = [20.0, 10.0, 0.0]
porosity = 0.3
sy = 0.15
rhob = 1600.0
kd = 1.0e-4
cinit = 100.0
lam = 0.01
perlen = 10.0

# the water table falls below the bottom of the top layer and then returns
heads = [30.0, 22.0, 14.0, 5.0, 5.0, 14.0, 22.0, 30.0, 30.0]
nper = len(heads)
DRY = -1.0e29


def build_model(ws, exe, decay, newton):
    sim = flopy.mf6.MFSimulation(sim_name="d", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, 1, 1.0)] * nper)
    solver = {
        "complexity": "complex",
        "outer_dvclose": 1e-8,
        "inner_dvclose": 1e-9,
        "linear_acceleration": "bicgstab",
        "outer_maximum": 200,
        "inner_maximum": 300,
    }

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="f",
        save_flows=True,
        newtonoptions="NEWTON" if newton else None,
    )
    # rewetting is not used with the Newton formulation, which keeps the cell
    # active when it drains
    rewet = {}
    if not newton:
        rewet = {
            "rewet_record": [("WETFCT", 1.0, "IWETIT", 1, "IHDWET", 0)],
            "wetdry": 1.0,
        }
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, filename="f.ims", **solver), ["f"]
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=1,
        k=50.0,
        **rewet,
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=0.0, sy=sy, transient={0: True}, save_flows=True
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={
            k: [[(nlay - 1, 0, 1), h, 0.0]] for k, h in enumerate(heads)
        },
        auxiliary="CONCENTRATION",
        pname="CHD-1",
        save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="f.hds",
        budget_filerecord="f.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="t", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, filename="t.ims", **solver), ["t"]
    )
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtadv(gwt)
    kwargs = {
        "porosity": porosity,
        "save_flows": True,
        "sorption": "linear",
        "bulk_density": rhob,
        "distcoef": kd,
        "stranded_mass": True,
        "stranded_filerecord": [("t.strand.bin",)],
    }
    if decay:
        kwargs.update(first_order_decay=True, decay=lam, decay_sorbed=lam)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="t.ucn",
        budget_filerecord="t.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="f", exgmnameb="t")
    sim.write_simulation(silent=True)
    return sim


@pytest.mark.parametrize("newton", [False, True], ids=["standard", "newton"])
@pytest.mark.parametrize("decay", [False, True])
def test_stranded_dry_cell(function_tmpdir, targets, decay, newton):
    ws = function_tmpdir
    sim = build_model(ws, str(targets["mf6"]), decay, newton)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"simulation failed\n{buff}"

    hf = flopy.utils.HeadFile(ws / "f.hds")
    sf = flopy.utils.HeadFile(ws / "t.strand.bin", text="STRANDED")
    times = sf.get_times()
    # the tested cell is the top layer of the column without the boundary
    head = np.array([hf.get_data(totim=t)[0, 0, 0] for t in times])
    st = np.array([sf.get_data(totim=t)[0, 0, 0] for t in times])

    # the drained cell is flagged dry with the standard formulation, and has a
    # head at or below its bottom with the Newton formulation
    dry = head < DRY if not newton else head <= botm[0]
    assert dry.any(), "the cell never went dry, so this tests nothing"
    assert not dry[-1], "the cell never came back, so the return is not tested"

    idry = np.flatnonzero(dry)
    first, last = idry[0], idry[-1]
    assert st[first - 1] > 0.0, "nothing was stranded before the cell went dry"

    # with the standard formulation the cell is inactive from the step in
    # which it goes dry; with the Newton formulation it stays active, so that
    # step strands what remained in it and the reservoir is held from its end
    if newton:
        ref, steps = first, idry[1:]
    else:
        ref, steps = first - 1, idry
    assert len(steps) > 0, "the cell was not dry long enough to test"

    # the reservoir is held through the dry period, decaying only if the
    # solute decays
    held = st[steps]
    if decay:
        expected = st[ref] * np.exp(-lam * perlen * np.arange(1, len(steps) + 1))
        assert np.allclose(held, expected, rtol=1e-6), (
            f"stranded mass did not decay at the given rate while the cell was "
            f"dry: {held} against {expected}"
        )
    else:
        assert np.allclose(held, st[ref], rtol=1e-9), (
            f"stranded mass changed while the cell was dry: {held}"
        )

    # the step that reactivates the cell must return mass, not strand more
    assert st[last + 1] < held[-1], (
        "the cell was rewet but its stranded mass grew, so the step was read "
        f"as drainage: {held[-1]} then {st[last + 1]}"
    )
    # the mass returns in proportion to the part of the drained interval that
    # resaturates, relative to the part drained when the cell stopped draining,
    # so what remains at the end is the share of that part still unsaturated,
    # less decay. The saturation is the one the flow model used.
    cbc = flopy.utils.CellBudgetFile(ws / "f.cbc", precision="double")
    sat = np.array([r["sat"][r["node"] == 1][0] for r in cbc.get_data(text="DATA-SAT")])
    drained = 1.0 - sat[ref]
    expected = st[last] * max(drained - sat[-1], 0.0) / drained
    if decay:
        expected *= np.exp(-lam * perlen * (len(st) - 1 - last))
    assert np.isclose(st[-1], expected, rtol=1e-6, atol=st.max() * 1e-9), (
        f"stranded mass after the cell was rewet {st[-1]} expected {expected}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
