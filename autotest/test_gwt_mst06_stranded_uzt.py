"""
Tests the report that the MST stranded mass option makes when the Unsaturated
Zone Transport (UZT) Package is active in the same model.

Both represent the part of a cell that the water table has drained, so a cell
in both carries its solute twice. UZF covers only part of a model in many
simulations, and stranded mass is meaningful in the rest of it, so the
combination is reported rather than refused.

Cases:
  - stranded : STRANDED_MASS with UZT runs and writes the warning.
  - plain    : the same model without STRANDED_MASS runs and does not.
"""

import flopy
import pytest

nlay, nrow, ncol = 5, 1, 1
delr = delc = 1.0
delv = 2.0
top = 0.0
botm = [top - (k + 1) * delv for k in range(nlay)]
strt = -8.0
seconds_to_days = 60.0 * 60.0 * 24.0
hk = 4.0e-6 * seconds_to_days
thts, thtr = 0.4, 0.2
brooks_corey_epsilon = 3.5


def build_model(ws, exe, stranded):
    sim = flopy.mf6.MFSimulation(sim_name="u", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(10.0, 10, 1.0)])

    gwf = flopy.mf6.ModflowGwf(sim, modelname="f", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(
            sim, complexity="simple", linear_acceleration="bicgstab", filename="f.ims"
        ),
        ["f"],
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=1,
        k=hk,
    )
    flopy.mf6.ModflowGwfsto(
        gwf, save_flows=True, iconvert=1, ss=0.0, sy=0.4, transient={0: True}
    )
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data={0: [[(nlay - 1, 0, 0), strt, hk / (0.5 * delv), 0.0]]},
        pname="GHB-1",
        auxiliary="CONCENTRATION",
    )
    pkdat = [
        [
            k,
            (k, 0, 0),
            1 if k == 0 else 0,
            k + 1 if k < nlay - 1 else -1,
            0.05,
            hk,
            thtr,
            thts,
            thtr,
            brooks_corey_epsilon,
            f"uzf{k + 1}",
        ]
        for k in range(nlay)
    ]
    flopy.mf6.ModflowGwfuzf(
        gwf,
        nuzfcells=len(pkdat),
        ntrailwaves=15,
        nwavesets=40,
        boundnames=True,
        save_flows=True,
        packagedata=pkdat,
        perioddata={0: [[0, 0.5 * hk, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]},
        pname="UZF-1",
        budget_filerecord="f.uzf.bud",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="f.hds",
        budget_filerecord="f.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="t", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(
            sim, complexity="simple", linear_acceleration="bicgstab", filename="t.ims"
        ),
        ["t"],
    )
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=100.0)
    flopy.mf6.ModflowGwtadv(gwt)
    kwargs = {
        "porosity": 0.2,
        "save_flows": True,
        "sorption": "linear",
        "bulk_density": 1600.0,
        "distcoef": 1.0e-4,
    }
    if stranded:
        kwargs["stranded_mass"] = True
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["GHB-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtuzt(
        gwt,
        flow_package_name="UZF-1",
        boundnames=True,
        packagedata=[(k, 0.0, f"uzt{k + 1}") for k in range(nlay)],
        uztperioddata=[(0, "INFILTRATION", 0.0)],
        pname="UZT-1",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="t.ucn",
        budget_filerecord="t.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="f", exgmnameb="t")
    sim.write_simulation(silent=True)
    return sim


def warned(ws):
    # the message is wrapped in the listing file, so compare without the breaks
    listing = " ".join((ws / "mfsim.lst").read_text().split())
    return "UZT Package is active in the same model" in listing


@pytest.mark.parametrize("stranded", [True, False])
def test_stranded_with_uzt(function_tmpdir, targets, stranded):
    ws = function_tmpdir / ("stranded" if stranded else "plain")
    ws.mkdir(parents=True)
    sim = build_model(ws, str(targets["mf6"]), stranded)
    success, buff = sim.run_simulation(silent=True)

    # the combination is reported, not refused, so the run must still finish
    assert success, f"simulation failed\n{buff}"
    assert warned(ws) is stranded, (
        f"the report was {'missing' if stranded else 'made'} for a model "
        f"{'with' if stranded else 'without'} STRANDED_MASS"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
