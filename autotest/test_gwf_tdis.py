"""Test TDIS package."""

import flopy
import numpy as np
import pytest
from modflow_devtools.markers import requires_pkg


@pytest.fixture
def simple_sim(tmp_path):
    """Create a simple and bogus GWF simulation without TDIS package."""
    sim = flopy.mf6.MFSimulation(sim_ws=str(tmp_path))
    _ = flopy.mf6.ModflowTdis(sim)  # placeholder
    _ = flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim)
    _ = flopy.mf6.ModflowGwfdis(gwf)
    _ = flopy.mf6.ModflowGwfic(gwf)
    _ = flopy.mf6.ModflowGwfnpf(gwf)
    sim.write_simulation()

    try:
        sim.remove_package("TDIS")
    except AttributeError:
        pass

    return sim


@requires_pkg("xmipy")
@pytest.mark.parametrize("tsmult", [1.0, 1.2])
def test_tdis_tsmult(tsmult, simple_sim, targets):
    """Check totim values to ensure they avoid accumulation errors."""
    from xmipy import XmiWrapper

    sim = simple_sim

    # Add TDIS package using time variables
    nper = 4
    nstp = 3
    perlen = 7.0
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)] * nper,
    )
    tdis.write()

    # Run within libmf6
    mf6 = XmiWrapper(lib_path=targets["libmf6"], working_directory=sim.sim_path)

    mf6.initialize()
    dt_list = []
    totim = mf6.get_current_time()
    all_times = [totim]
    sp_times = [totim]
    for kper in range(nper):
        for kstep in range(nstp):
            mf6.update()
            delt = mf6.get_time_step()
            totim = mf6.get_current_time()
            all_times.append(totim)
            dt_list.append(delt)
        sp_times.append(totim)
    mf6.finalize()

    # Stress period times should be exact
    np.testing.assert_equal(sp_times, np.arange(nper + 1) * perlen)

    if tsmult == 1.0:
        # Note that delt may not always be exactly unique, but should be close
        assert max(dt_list) - min(dt_list) < 1e-14
        assert pytest.approx(dt_list[0]) == perlen / nstp

        # Check all times
        np.testing.assert_allclose(
            all_times,
            np.linspace(0.0, perlen * nper, nper * nstp + 1),
        )
    else:
        # Recreate geometric progression to check delt
        dt0 = perlen * (1.0 - tsmult) / (1.0 - tsmult**nstp)
        expected_dt = np.tile(dt0 * tsmult ** np.arange(nstp), nper)
        np.testing.assert_allclose(dt_list, expected_dt)

        # Check all times
        np.testing.assert_allclose(
            all_times,
            np.cumsum(np.insert(expected_dt, 0, 0.0)),
        )
