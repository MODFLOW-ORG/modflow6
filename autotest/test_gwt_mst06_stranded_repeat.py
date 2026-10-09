"""
Tests that the MST stranded mass option gives the same result when a time step
that failed to converge is solved again.

The adaptive time step (ATS) utility solves such a time step again with a
shorter length, so the stranded mass, like the concentration, has to start
again from the values the previous time step ended with.

Cases:
  - ats_retry : a cell is pumped to half saturation and refilled with clean
                water to 0.9, with a transport solution allowed one outer
                iteration, so that ATS repeats the first time step of the
                refill with shorter lengths; the stranded mass and the final
                concentration match the analytical values, which do not depend
                on the time step. A last stress period without stresses
                follows, so that the repeated time step is not the last of the
                simulation.
"""

import flopy
import numpy as np

top, botm = 10.0, 0.0
delr = delc = 10.0
vcell = delr * delc * (top - botm)
porosity = 0.1
rhob, kd = 1500.0, 1.0e-4
cinit = 100.0


def series(ws, fname, text):
    f = flopy.utils.HeadFile(ws / fname, text=text, precision="double")
    return np.array([f.get_data(totim=t).ravel() for t in f.get_times()])


def build_ats(ws, exe):
    # pump 25 m3/d for 2 d, inject 25 m3/d of clean water for 1.6 d, and leave
    # the cell for 1 d
    sim = flopy.mf6.MFSimulation(sim_name="ats", sim_ws=ws, exe_name=exe)
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=3,
        perioddata=[(2.0, 1, 1.0), (1.6, 1, 1.0), (1.0, 1, 1.0)],
    )
    # each period starts with a single step; a step that fails is repeated
    # with a quarter of its length, which is then kept for the period
    tdis.ats.initialize(
        maxats=2,
        perioddata=[(k, 2.0, 1.0e-4, 2.0, 1.0, 4.0) for k in range(2)],
        filename="ats.ats",
    )

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", save_flows=True, newtonoptions="NEWTON"
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-10,
        inner_dvclose=1e-12,
        linear_acceleration="bicgstab",
        filename="gwf.ims",
    )
    sim.register_ims_package(ims, ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0, save_saturation=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=porosity, transient={0: True})
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={
            0: [[(0, 0, 0), -25.0, 0.0]],
            1: [[(0, 0, 0), 25.0, 0.0]],
            2: [[(0, 0, 0), 0.0, 0.0]],
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

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    # one outer iteration converges only if the concentration changes by
    # less than 5 g/m3, so the long time steps of the refill fail
    ims = flopy.mf6.ModflowIms(
        sim,
        outer_maximum=1,
        outer_dvclose=5.0,
        inner_dvclose=1e-12,
        linear_acceleration="bicgstab",
        filename="gwt.ims",
    )
    sim.register_ims_package(ims, ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        save_flows=True,
        sorption="linear",
        bulk_density=rhob,
        distcoef=kd,
        stranded_mass=True,
        stranded_filerecord=[("gwt.strand.bin",)],
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def test_ats_retry(function_tmpdir, targets):
    ws = function_tmpdir
    sim = build_ats(ws, str(targets["mf6"]))
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"simulation failed\n{buff}"

    listing = (ws / "mfsim.lst").read_text()
    assert "Failed solution for step" in listing, "no time step was repeated"

    f = flopy.utils.HeadFile(ws / "gwt.strand.bin", text="STRANDED")
    times = np.array(f.get_times())
    stranded = np.array([f.get_data(totim=t).ravel()[0] for t in times])
    conc = series(ws, "gwt.ucn", "CONCENTRATION")[:, 0]
    # the half of the cell that drains strands its sorbed mass, 7,500 g, and
    # the clean water rewets 0.4 of the 0.5 that drained, returning four
    # fifths of it; the other 18,500 g are in the water and on the solids of
    # the cell at a saturation of 0.9
    pumped = times <= 2.0
    assert np.isclose(stranded[pumped][-1], 7500.0), f"stranded {stranded}"
    assert np.isclose(times[-1], 4.6), f"the simulation ended at {times[-1]}"
    assert np.isclose(stranded[-1], 1500.0), f"stranded {stranded}"
    retard = porosity + rhob * kd
    expected = 18500.0 / (retard * 0.9 * vcell)
    assert np.isclose(conc[-1], expected), (
        f"concentration {conc[-1]} expected {expected}"
    )
