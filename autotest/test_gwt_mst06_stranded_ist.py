"""
Tests MST stranded mass in a model with an immobile domain (IST Package), using
the single cell of test_gwt_mst06_stranded_cell.py.

The immobile domain holds its solute and sorbate over the whole cell whatever
the saturation, and exchanges solute with the mobile domain in proportion to the
saturation, so a cell that drains strands mobile mass alone. The specific yield
is half of the mobile porosity, so both the solute of the retained water and
the sorbate are stranded. The mass of the cell, mobile, immobile, and stranded,
is computed from the concentrations and the stranded mass written to file.

Cases:
  - conserved : a well drains half of the cell and refills it with clean water,
                with and without linear sorption; stranded mass is simulated,
                both budgets close, and the mass of the cell changes only by
                the mass the well removes. With LEGACY_STORAGE the same check
                fails, by up to 2,300 g of 25,000 g.
"""

import os

import flopy
import numpy as np
import pytest
from test_gwt_mst06_stranded_cell import (
    assert_budget_closes,
    cinit,
    get_model,
    kd,
    porosity,
    rhob,
    top,
    vcell,
)

# immobile domain: a fifth of the cell, with a porosity of 0.1, exchanging
# solute with the mobile domain at a rate coefficient of 0.1 per day
volfrac = 0.2
porosity_im = 0.1
zetaim = 0.1
thetam = (1.0 - volfrac) * porosity
thetaim = volfrac * porosity_im
sy = 0.5 * thetam
# 20 m3 drains half of the cell, 4 m3 of pore space per meter of head
periods = [(2.0, [(-10.0, 0.0)]), (2.0, [(10.0, 0.0)])]


def build_model(ws, sorb, legacy, nstp):
    sim = get_model(
        ws,
        periods,
        nstp,
        sy,
        sorb,
        strand=not legacy,
    )
    gwt = sim.get_model("gwt")
    kwargs = {
        "porosity": porosity_im,
        "volfrac": volfrac,
        "zetaim": zetaim,
        "cim": cinit,
        "cim_filerecord": "gwt.cim",
    }
    if sorb:
        kwargs.update(sorption="linear", bulk_density=rhob, distcoef=kd)
    flopy.mf6.ModflowGwtist(gwt, **kwargs)
    return sim


def cell_mass(ws, sorb, legacy):
    """Mass of the cell, mobile, immobile, and stranded, at each time."""
    hf = flopy.utils.HeadFile(os.path.join(ws, "gwf.hds"))
    times = np.array(hf.get_times())
    head = np.array([hf.get_data(totim=t).ravel()[0] for t in times])
    sat = np.clip(head / top, 0.0, 1.0)
    cf = flopy.utils.HeadFile(os.path.join(ws, "gwt.ucn"), text="CONCENTRATION")
    conc = np.array([cf.get_data(totim=t).ravel()[0] for t in times])
    imf = flopy.utils.HeadFile(os.path.join(ws, "gwt.cim"), text="CIM")
    cim = np.array([imf.get_data(totim=t).ravel()[0] for t in times])
    stranded = np.zeros_like(times)
    if not legacy:
        sf = flopy.utils.HeadFile(os.path.join(ws, "gwt.strand.bin"), text="STRANDED")
        stranded = np.array([sf.get_data(totim=t).ravel()[0] for t in times])
    rkm = (1.0 - volfrac) * rhob * kd if sorb else 0.0
    rkim = volfrac * rhob * kd if sorb else 0.0
    mobile = (thetam + rkm) * sat * vcell * conc
    immobile = (thetaim + rkim) * vcell * cim
    mass0 = (thetam + rkm + thetaim + rkim) * vcell * cinit
    return times, mass0, mobile + immobile + stranded


def well_mass(ws):
    """Times of the printed budgets and the mass the wells added by then."""
    lb = flopy.utils.Mf6ListBudget(
        os.path.join(ws, "gwt.lst"), budgetkey="MASS BUDGET FOR ENTIRE MODEL"
    )
    cum = lb.get_dataframes(diff=False, start_datetime=None)[1]
    return cum.index.to_numpy(), (cum["WEL_IN"] - cum["WEL_OUT"]).to_numpy()


@pytest.mark.parametrize("nstp", [1, 10])
@pytest.mark.parametrize("sorb", [False, True])
def test_conserved(function_tmpdir, targets, sorb, nstp):
    exe = str(targets["mf6"])
    error = {}
    for legacy in (False, True):
        ws = function_tmpdir / ("legacy" if legacy else "stranded")
        sim = build_model(ws, sorb, legacy, nstp)
        sim.exe_name = exe
        sim.write_simulation(silent=True)
        success, buff = sim.run_simulation(silent=True)
        assert success, f"simulation failed\n{buff}"
        listing = " ".join((ws / "gwt.lst").read_text().split())
        if not legacy:
            assert "STRANDED MASS IS ACTIVE" in listing, "stranded mass is off"
            assert_budget_closes(ws)
        times, mass0, mass = cell_mass(ws, sorb, legacy)
        # the listing file prints the budget at the end of each period
        tbud, wel = well_mass(ws)
        mass = mass[np.isin(times, tbud)]
        error[legacy] = mass - (mass0 + wel)
    tol = 1e-6 * mass0
    assert np.all(np.abs(error[False]) < tol), (
        f"the mass of the cell is not conserved with stranded mass {error[False]}"
    )
    assert np.abs(error[True]).max() > 1.0, (
        f"LEGACY_STORAGE conserves the mass of the cell {error[True]}, so the "
        "check does not discriminate"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
