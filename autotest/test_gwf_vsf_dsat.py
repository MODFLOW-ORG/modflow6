"""
Test the NEWTON_DSAT_MAX option of the VSF package, which limits the change
of effective saturation per outer iteration in the Newton formulation.

Infiltration into a dry Haverkamp column (Celia et al., 1990) is solved with
the limit switched off, with the default (0.2) and with a tight limit (0.05).
The limit only changes the path of the outer iterations, not the equations,
so all three give the same heads; the tight limit takes more iterations.
"""

import flopy
import numpy as np
from gwf_test_utils import get_vsf_soil_data

# (name, NEWTON_DSAT_MAX), None leaves the option out to get the default
cases = [("dsat-off", 0.0), ("dsat-dflt", None), ("dsat-005", 0.05)]

column_height = 42.0
nlay = 42
dt = 30.0
perlen = 360.0
dvclose = 1.0e-8


def build_model(ws, name, dsat_max, exe):
    delz = column_height / nlay
    top = delz
    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
    hp_upper, hp_lower = -20.7, -61.5
    h_upper = hp_upper + botm[0] + 0.5 * delz
    h_lower = hp_lower + botm[-1] + 0.5 * delz
    strt = np.array([botm[i] + 0.5 * delz + hp_lower for i in range(nlay)])

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim, time_units="SECONDS", perioddata=[(perlen, int(perlen / dt), 1.0)]
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=dvclose,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=300,
        inner_dvclose=dvclose,
        rcloserecord=1.0e-8,
        linear_acceleration="BICGSTAB",
        csv_outer_output_filerecord=f"{name}.outer.csv",
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=1, ncol=1, delr=10.0, delc=10.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=0.00944, k33=0.00944)
    flopy.mf6.ModflowGwfsto(gwf, ss=0.0, sy=0.0, iconvert=0, transient={0: True})
    soil = get_vsf_soil_data("Celia1990-eq10-Haverkamp")
    opts = {} if dsat_max is None else {"newton_dsat_max": dsat_max}
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        soil_model="Haverkamp",
        porosity=soil["porosity"],
        satres=soil["satres"],
        alphahvk=soil["alpha"],
        nhvk=soil["n"],
        betahvk=soil["beta"],
        khvk=soil["k"],
        **opts,
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), h_upper], [(nlay - 1, 0, 0), h_lower]]}
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim


def test_vsf_newton_dsat_max(function_tmpdir, targets):
    heads, outer = {}, {}
    for name, dsat_max in cases:
        ws = function_tmpdir / name
        sim = build_model(ws, name, dsat_max, targets["mf6"])
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True)
        assert success, f"{name} did not run successfully"

        lst = (ws / f"{name}.lst").read_text()
        if dsat_max == 0.0:
            assert "saturation change limit switched off" in lst
        else:
            expected = 0.2 if dsat_max is None else dsat_max
            assert "saturation change per outer iteration limited to" in lst
            line = next(s for s in lst.splitlines() if "limited to" in s)
            assert np.isclose(float(line.split()[-1]), expected)

        heads[name] = flopy.utils.HeadFile(ws / f"{name}.hds").get_data().ravel()
        csv = np.genfromtxt(ws / f"{name}.outer.csv", delimiter=",", names=True)
        outer[name] = csv.size

    # same equations, so the same answer, whatever the limit
    for name in ("dsat-dflt", "dsat-005"):
        dh = np.abs(heads[name] - heads["dsat-off"]).max()
        assert dh < 1.0e-5, f"{name} heads differ from unlimited run by {dh}"

    # the tight limit must have been active: it shortens steps across the front
    assert outer["dsat-005"] > outer["dsat-off"], outer
