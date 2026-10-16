"""
Tests the MST stranded mass option when the flow model is run as a separate
simulation and the transport model reads its flows through the FMI Package.

The retained water volume is derived from the STO-SY flow term, which a
coupled flow model always supplies but a budget file carries only if the flow
model saved it.

Cases:
  - matches_coupled : results are identical to the same models run in one
                      simulation through a GWF-GWT exchange.
  - requires_stosy  : a budget file without STO-SY turns stranded mass off with
                      a note, and results match LEGACY_STORAGE.
  - confined_flow   : a confined flow model writes no STO-SY because it has no
                      specific yield, and is handled the same way.
"""

import flopy
import numpy as np
import pytest

top, botm = 40.0, 0.0
delr = delc = 100.0
porosity = 0.3
sy = 0.15
rhob = 1600.0
kd = 1.0e-4
cinit = 100.0

# the leading period does not drain, which lets the saturation of the
# previous time step initialize before any mass is stranded
heads = [40.0, 35.0, 30.0, 25.0, 20.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
nper = len(heads)
perioddata = [(10.0, 1, 1.0)] * nper


def add_ims(sim, filename):
    return flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=filename,
    )


def add_gwf(sim, name, save_sto, confined=False):
    # SAVE_FLOWS on the model writes the storage terms whatever the STO
    # package asks for, so suppressing STO-SY means turning it off there too;
    # NPF and CHD keep theirs, so the budget file still has everything else
    # the transport model needs
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=save_sto, newtonoptions="NEWTON"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=0 if confined else 1,
        k=100.0,
    )
    # a confined cell has no specific yield, so STO writes STO-SS but no STO-SY
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=0 if confined else 1,
        ss=1.0e-5 if confined else 0.0,
        sy=sy,
        transient={0: True},
        save_flows=save_sto,
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={k: [[(0, 0, 1), h, 0.0]] for k, h in enumerate(heads)},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
        save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def add_gwt(sim, name, fmi_packagedata=None, legacy=False, strand_output=True):
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=2, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        save_flows=True,
        sorption="linear",
        bulk_density=rhob,
        distcoef=kd,
        legacy_storage=legacy,
        stranded_filerecord=[(f"{name}.strand.bin",)] if strand_output else None,
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    if fmi_packagedata is not None:
        flopy.mf6.ModflowGwtfmi(gwt, packagedata=fmi_packagedata)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{name}.ucn",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    return gwt


def run_coupled(ws, exe):
    """Flow and transport in one simulation through a GWF-GWT exchange."""
    sim = flopy.mf6.MFSimulation(sim_name="c", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
    add_gwf(sim, "f", True)
    sim.register_ims_package(add_ims(sim, "f.ims"), ["f"])
    add_gwt(sim, "t")
    sim.register_ims_package(add_ims(sim, "t.ims"), ["t"])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="f", exgmnameb="t")
    sim.write_simulation(silent=True)
    return sim.run_simulation(silent=True)


def run_flow(ws, exe, save_sto, confined=False):
    sim = flopy.mf6.MFSimulation(sim_name="f", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
    add_gwf(sim, "f", save_sto, confined)
    add_ims(sim, "f.ims")
    sim.write_simulation(silent=True)
    return sim.run_simulation(silent=True)


def run_transport(ws, exe, flow_ws, legacy=False, strand_output=True):
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
    add_gwt(
        sim,
        "t",
        [
            ("GWFHEAD", str(flow_ws / "f.hds"), None),
            ("GWFBUDGET", str(flow_ws / "f.cbc"), None),
        ],
        legacy=legacy,
        strand_output=strand_output,
    )
    add_ims(sim, "t.ims")
    sim.write_simulation(silent=True)
    return sim.run_simulation(silent=True)


def assert_stosy_fallback(tran, exe, flow, buff):
    """Stranded mass is turned off with a note, and the run matches the legacy one."""
    # the message is wrapped in the listing file, so compare without the breaks
    listing = " ".join((tran / "mfsim.lst").read_text().split())
    assert "does not contain the STO-SY flow term" in listing, (
        f"the note was not written\n{buff}"
    )
    legacy = tran.parent / "legacy"
    legacy.mkdir()
    success, buff = run_transport(legacy, exe, flow, legacy=True, strand_output=False)
    assert success, f"legacy transport simulation failed\n{buff}"
    a = series(tran, "t.ucn", "CONCENTRATION")
    b = series(legacy, "t.ucn", "CONCENTRATION")
    assert np.array_equal(a, b), (
        f"results differ from LEGACY_STORAGE by {np.abs(a - b).max()}"
    )


def series(ws, fname, text):
    f = flopy.utils.HeadFile(ws / fname, text=text)
    return np.array([f.get_data(totim=t).flatten() for t in f.get_times()])


def test_stranded_fmi_matches_coupled(function_tmpdir, targets):
    exe = str(targets["mf6"])
    coupled = function_tmpdir / "coupled"
    flow = function_tmpdir / "flow"
    tran = function_tmpdir / "transport"
    for d in (coupled, flow, tran):
        d.mkdir(parents=True)

    success, buff = run_coupled(coupled, exe)
    assert success, f"coupled simulation failed\n{buff}"
    success, buff = run_flow(flow, exe, True)
    assert success, f"flow simulation failed\n{buff}"
    success, buff = run_transport(tran, exe, flow)
    assert success, f"transport simulation failed\n{buff}"

    # the flow model supplies the same saturations and storage flows either
    # way, so the two must agree to the last bit
    for fname, text in (("t.ucn", "CONCENTRATION"), ("t.strand.bin", "STRANDED")):
        a = series(coupled, fname, text)
        b = series(tran, fname, text)
        assert np.array_equal(a, b), (
            f"{text} differs between the coupled and the separate simulation, "
            f"maximum difference {np.abs(a - b).max()}"
        )
    # the option did something, so the comparison above is not vacuous
    assert series(tran, "t.strand.bin", "STRANDED").max() > 0.0


def test_stranded_fmi_requires_stosy(function_tmpdir, targets):
    exe = str(targets["mf6"])
    flow = function_tmpdir / "flow"
    tran = function_tmpdir / "transport"
    for d in (flow, tran):
        d.mkdir(parents=True)

    # the flow model does not save its storage flows, so the budget file has
    # no STO-SY and the retained water volume cannot be determined
    success, buff = run_flow(flow, exe, False)
    assert success, f"flow simulation failed\n{buff}"
    names = flopy.utils.CellBudgetFile(
        flow / "f.cbc", precision="double"
    ).get_unique_record_names(decode=True)
    assert "STO-SY" not in [n.strip() for n in names], (
        "the flow model saved STO-SY, so this case does not test what it claims"
    )

    success, buff = run_transport(tran, exe, flow, strand_output=False)
    assert success, f"transport simulation failed\n{buff}"
    assert_stosy_fallback(tran, exe, flow, buff)


def test_stranded_fmi_confined_flow(function_tmpdir, targets):
    exe = str(targets["mf6"])
    flow = function_tmpdir / "flow"
    tran = function_tmpdir / "transport"
    for d in (flow, tran):
        d.mkdir(parents=True)

    # a confined flow model saves its storage flows but has no specific yield
    success, buff = run_flow(flow, exe, True, confined=True)
    assert success, f"flow simulation failed\n{buff}"
    names = [
        n.strip()
        for n in flopy.utils.CellBudgetFile(
            flow / "f.cbc", precision="double"
        ).get_unique_record_names(decode=True)
    ]
    assert "STO-SS" in names, "the flow model saved no storage flows at all"
    assert "STO-SY" not in names, (
        "the confined flow model wrote STO-SY, so this case does not test "
        "what it claims"
    )

    # no cell drains, but a budget file cannot show that, so it is handled as
    # any budget file without STO-SY
    success, buff = run_transport(tran, exe, flow, strand_output=False)
    assert success, f"transport simulation failed\n{buff}"
    assert_stosy_fallback(tran, exe, flow, buff)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
