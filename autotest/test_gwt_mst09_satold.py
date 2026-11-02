"""
Tests that the MST storage terms carry solute across time steps in a single
cell whose water table moves, using the saturation of the previous time step.

Specific yield equals porosity and specific storage is zero, so the water
volume of the transport model matches the storage of the flow model, and the
concentration at the end of each time step follows from a mass balance of the
cell over the time step.

Cases:
  - drain      : a well drains half the cell and refills it with clean water,
                 with and without linear sorption and with one or four time
                 steps per period; no solute is lost.
  - fmi_file   : the drain case with flows read from the files of a separate
                 flow simulation, where the first time step reconstructs the
                 starting saturation from the specific yield storage flow.
  - ats_retry  : the drain case, starting at a saturation of 0.9, with a
                 transport solution that fails on long time steps, so that ATS
                 repeats them with shorter lengths.
  - inject     : a cell is filled with water at Cin, given a mass M of pure
                 solute, and partly drained (Provost, written communication,
                 2026); the concentration does not depend on whether solute
                 enters with the water or as pure solute, and is unchanged
                 when Cin equals the starting concentration and M is zero
                 without sorption.
  - confined   : a confined cell with specific storage takes in and gives up
                 water at the starting concentration, with linear sorption;
                 the concentration does not change, since the cell stays
                 saturated.
"""

import flopy
import numpy as np
import pytest

top, botm = 10.0, 0.0
delr = delc = 10.0
vcell = delr * delc * (top - botm)
porosity = 0.1
rhob, kd = 1500.0, 1.0e-4
cinit = 100.0
# pump 25 m3/d for 2 d to a saturation of 0.5, then inject clean water at
# 25 m3/d for 2 d; each entry is (length, well rate, injected concentration,
# mass source rate)
drain_periods = [(2.0, -25.0, 0.0, 0.0), (2.0, 25.0, 0.0, 0.0)]


def retardation(sorb):
    return porosity + (rhob * kd if sorb else 0.0)


def add_ims(sim, fname, **kwargs):
    # the moderate settings damp the change in storage coefficient at the top
    # of the cell
    options = {
        "complexity": "moderate",
        "outer_dvclose": 1e-10,
        "inner_dvclose": 1e-12,
        "linear_acceleration": "bicgstab",
    }
    options.update(kwargs)
    return flopy.mf6.ModflowIms(sim, filename=fname, **options)


def add_gwf(sim, periods, strt, ss=0.0):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", save_flows=True, newtonoptions="NEWTON"
    )
    sim.register_ims_package(add_ims(sim, "gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=0 if ss > 0.0 else 1,
        k=10.0,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=0 if ss > 0.0 else 1, ss=ss, sy=porosity, transient={0: True}
    )
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={
            k: [[(0, 0, 0), q, cin]] for k, (_, q, cin, _) in enumerate(periods)
        },
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def add_gwt(sim, periods, sorb, fmi_files=False, ims_kwargs=None):
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    ims = add_ims(sim, "gwt.ims", **(ims_kwargs or {}))
    sim.register_ims_package(ims, ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    kwargs = {"porosity": porosity, "save_flows": True}
    if sorb:
        kwargs.update(sorption="linear", bulk_density=rhob, distcoef=kd)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    if any(m != 0.0 for *_, m in periods):
        flopy.mf6.ModflowGwtsrc(
            gwt,
            stress_period_data={
                k: [[(0, 0, 0), m]] for k, (*_, m) in enumerate(periods)
            },
        )
    if fmi_files:
        flopy.mf6.ModflowGwtfmi(
            gwt,
            packagedata=[
                ("GWFHEAD", "../flow/gwf.hds"),
                ("GWFBUDGET", "../flow/gwf.cbc"),
            ],
        )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    return gwt


def build(ws, exe, periods, nstp, sorb, strt=top, ss=0.0, ims_kwargs=None):
    """Coupled flow and transport simulation of a single cell."""
    sim = flopy.mf6.MFSimulation(sim_name="cell", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(perlen, nstp, 1.0) for perlen, *_ in periods],
    )
    add_gwf(sim, periods, strt, ss=ss)
    add_gwt(sim, periods, sorb, ims_kwargs=ims_kwargs)
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def run(sim):
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"simulation failed\n{buff}"


def simulated(ws, gwf_ws=None):
    """Times, saturations, and concentrations at the end of each time step."""
    gwf_ws = ws if gwf_ws is None else gwf_ws
    cbc = flopy.utils.CellBudgetFile(gwf_ws / "gwf.cbc", precision="double")
    ucn = flopy.utils.HeadFile(ws / "gwt.ucn", text="CONCENTRATION")
    times = np.array(ucn.get_times())
    sat = np.array([cbc.get_data(text="DATA-SAT", totim=t)[0]["sat"][0] for t in times])
    conc = np.array([ucn.get_data(totim=t).ravel()[0] for t in times])
    return times, sat, conc


def expected(times, sat, periods, sorb, sat0):
    """Concentration from a mass balance of the cell over each time step.

    The well and the mass source act over the period that holds the end of
    the time step, and water leaving the cell carries the concentration at the
    end of the time step.
    """
    ends = np.cumsum([perlen for perlen, *_ in periods])
    rv = retardation(sorb) * vcell
    conc = []
    t0, s0, c0 = 0.0, sat0, cinit
    for t, s in zip(times, sat):
        _, q, cin, m = periods[int(np.searchsorted(ends, t - 1e-9))]
        dt = t - t0
        qin, qout = max(q, 0.0), max(-q, 0.0)
        c = (rv * s0 * c0 + dt * (qin * cin + m)) / (rv * s + dt * qout)
        conc.append(c)
        t0, s0, c0 = t, s, c
    return np.array(conc)


@pytest.mark.parametrize("nstp", [1, 4])
@pytest.mark.parametrize("sorb", [False, True], ids=["nosorb", "sorb"])
def test_drain(function_tmpdir, targets, sorb, nstp):
    ws = function_tmpdir
    run(build(ws, str(targets["mf6"]), drain_periods, nstp, sorb))
    times, sat, conc = simulated(ws)
    assert np.isclose(sat[nstp - 1], 0.5), f"saturation {sat}"
    assert np.isclose(sat[-1], 1.0), f"saturation {sat}"
    ref = expected(times, sat, drain_periods, sorb, 1.0)
    assert np.allclose(conc, ref, rtol=1e-6), f"concentration {conc} expected {ref}"
    # clean water dilutes what is left after pumping, so the mass in the cell
    # is unchanged by the refill
    mass = retardation(sorb) * vcell * sat * conc
    assert np.isclose(mass[-1], mass[nstp - 1], rtol=1e-6), f"mass {mass}"


@pytest.mark.parametrize("sorb", [False, True], ids=["nosorb", "sorb"])
def test_fmi_file(function_tmpdir, targets, sorb):
    exe = str(targets["mf6"])
    flow_ws = function_tmpdir / "flow"
    sim = flopy.mf6.MFSimulation(sim_name="flow", sim_ws=flow_ws, exe_name=exe)
    tdis = [(perlen, 1, 1.0) for perlen, *_ in drain_periods]
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=len(drain_periods), perioddata=tdis
    )
    add_gwf(sim, drain_periods, top)
    run(sim)

    ws = function_tmpdir / "transport"
    sim = flopy.mf6.MFSimulation(sim_name="transport", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=len(drain_periods), perioddata=tdis
    )
    add_gwt(sim, drain_periods, sorb, fmi_files=True)
    run(sim)

    times, sat, conc = simulated(ws, flow_ws)
    ref = expected(times, sat, drain_periods, sorb, 1.0)
    assert np.allclose(conc, ref, rtol=1e-6), f"concentration {conc} expected {ref}"


def test_ats_retry(function_tmpdir, targets):
    ws = function_tmpdir
    # one outer iteration converges only if the concentration changes by less
    # than 5 g/m3, so long time steps fail and are repeated a quarter as long;
    # a last period without stresses follows, so that a repeated time step is
    # not the last of the simulation
    periods = drain_periods + [(1.0, 0.0, 0.0, 0.0)]
    sim = build(
        ws,
        str(targets["mf6"]),
        periods,
        1,
        True,
        strt=9.0,
        ims_kwargs={"complexity": "simple", "outer_maximum": 1, "outer_dvclose": 5.0},
    )
    sim.tdis.ats.initialize(
        maxats=2,
        perioddata=[(k, 2.0, 1.0e-4, 2.0, 1.0, 4.0) for k in range(2)],
        filename="cell.ats",
    )
    run(sim)
    assert "Failed solution for step" in (ws / "mfsim.lst").read_text(), (
        "no time step was repeated"
    )
    times, sat, conc = simulated(ws)
    assert np.isclose(times[-1], 5.0), f"the simulation ended at {times[-1]}"
    ref = expected(times, sat, periods, True, 0.9)
    assert np.allclose(conc, ref, rtol=1e-6), f"concentration {conc} expected {ref}"


# the cell starts at a saturation of 0.5, is filled to 0.7 by 20 m3 of water
# at Cin in 1 d, and is drained to 0.4 by pumping 30 m3 in 1 d while a mass M
# of pure solute is added; each case is (Cin, M)
inject_cases = {
    "mixed": (50.0, 2000.0),
    "same": (cinit, 0.0),
    "pure": (0.0, 2000.0 + 20.0 * 50.0),
}


@pytest.mark.parametrize("sorb", [False, True], ids=["nosorb", "sorb"])
@pytest.mark.parametrize("case", list(inject_cases))
def test_inject(function_tmpdir, targets, case, sorb):
    ws = function_tmpdir
    cin, mass = inject_cases[case]
    periods = [(1.0, 20.0, cin, 0.0), (1.0, -30.0, 0.0, mass)]
    run(build(ws, str(targets["mf6"]), periods, 1, sorb, strt=5.0))
    times, sat, conc = simulated(ws)
    assert np.allclose(sat, [0.7, 0.4]), f"saturation {sat}"
    ref = expected(times, sat, periods, sorb, 0.5)
    assert np.allclose(conc, ref, rtol=1e-6), f"concentration {conc} expected {ref}"
    if not sorb:
        # Provost's equation 12, from the water in the cell after the fill
        vs0, vin = porosity * 0.5 * vcell, 20.0
        fs = vin / (vs0 + vin)
        c2 = cinit + fs * (cin - cinit) + mass / (vs0 + vin)
        assert np.isclose(conc[-1], c2), f"concentration {conc[-1]} expected {c2}"
    if case == "same" and not sorb:
        assert np.allclose(conc, cinit), f"concentration {conc}"
    elif case == "pure":
        # the same solute as the mixed case, entering only as pure solute
        cin, mass = inject_cases["mixed"]
        mixed = [(1.0, 20.0, cin, 0.0), (1.0, -30.0, 0.0, mass)]
        c2 = expected(times, sat, mixed, sorb, 0.5)[-1]
        assert np.isclose(conc[-1], c2), f"concentration {conc[-1]} expected {c2}"


def test_confined(function_tmpdir, targets):
    ws = function_tmpdir
    # with a specific storage of 1e-3 1/m, injecting 1 m3 raises the head by 1 m
    # and pumping 2 m3 lowers it by 2 m
    periods = [(1.0, 1.0, cinit, 0.0), (1.0, -2.0, 0.0, 0.0)]
    run(build(ws, str(targets["mf6"]), periods, 1, True, strt=20.0, ss=1.0e-3))
    times, sat, conc = simulated(ws)
    assert np.allclose(sat, 1.0), f"saturation {sat}"
    assert np.allclose(conc, cinit, rtol=1e-9), f"concentration {conc}"
