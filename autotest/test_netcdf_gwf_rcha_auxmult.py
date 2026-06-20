"""
NetCDF structured export test for RCHA READASARRAYS with AUXMULTNAME and
time-varying recharge across multiple stress periods.

This test verifies that auxiliary multiplier arrays are correctly exported
to NetCDF and that when read back via the NETCDF keyword, the multiplier
persists across stress periods (same bug as the non-netcdf case).
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

xa = pytest.importorskip("xarray")
nc = pytest.importorskip("netCDF4")

nlay, nrow, ncol = 1, 5, 5
nper = 3
name = "rcha_auxmult_nc"

# recharge varies per period; multiplier only in period 1
rch_vals = [0.001, 0.002, 0.002]
mult_val = 0.5


def build_models(idx, test):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=str(test.workspace),
        exe_name=test.targets["mf6"],
    )

    tdis_rc = [(1.0, 1, 1.0)] * nper
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(sim)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=100.0,
        top=10.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)

    # fixed heads on right column
    chd_spd = [[(0, j, ncol - 1), 10.0] for j in range(nrow)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_spd})

    # recharge with auxmultname - multiplier only in period 1
    rch_spd = {i: np.full((nrow, ncol), v) for i, v in enumerate(rch_vals)}
    rch_mult = np.full((nrow, ncol), mult_val)
    rch_aux = {0: [rch_mult]}

    flopy.mf6.ModflowGwfrcha(
        gwf,
        readasarrays=True,
        recharge=rch_spd,
        auxiliary="RCH_MULT",
        auxmultname="RCH_MULT",
        aux=rch_aux,
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # netcdf export settings
    gwf.dis.export_array_netcdf = True
    gwf.ic.export_array_netcdf = True
    gwf.npf.export_array_netcdf = True
    gwf.rcha_0.export_array_netcdf = True

    gwf.name_file.nc_structured_filerecord = f"{name}.input.nc"
    flopy.mf6.ModflowUtlncf(gwf.dis, filename=f"{name}.dis.ncf")

    return sim, None


def check_output(idx, test):
    # --- verify validate-mode NetCDF input file ---
    with xa.open_dataset(test.workspace / f"{name}.input.nc") as xds:
        # check recharge variable
        assert "rcha_0_recharge" in xds, "rcha_0_recharge missing from NetCDF"
        rch_var = xds["rcha_0_recharge"]
        assert rch_var.shape[0] == nper, (
            f"Expected {nper} time slices for recharge, got {rch_var.shape[0]}"
        )
        for i, v in enumerate(rch_vals):
            assert np.allclose(rch_var[i].values, v), (
                f"Period {i + 1} recharge mismatch"
            )

        # check auxiliary multiplier variable
        assert "rcha_0_rch_mult" in xds, "rcha_0_rch_mult missing from NetCDF"
        aux_var = xds["rcha_0_rch_mult"]
        # multiplier provided in period 1 should be exported
        assert np.allclose(aux_var[0].values, mult_val), (
            f"Period 1 RCH_MULT mismatch: expected {mult_val}"
        )

    # --- re-run reading from NetCDF (Pass 2) ---
    with open(test.workspace / f"{name}.nam", "w") as f:
        f.write("BEGIN options\n")
        f.write("  SAVE_FLOWS\n")
        f.write(f"  NETCDF  FILEIN {name}.input.nc\n")
        f.write("END options\n\n")
        f.write("BEGIN packages\n")
        f.write(f"  DIS6  {name}.dis  dis\n")
        f.write(f"  IC6  {name}.ic  ic\n")
        f.write(f"  NPF6  {name}.npf  npf\n")
        f.write(f"  CHD6  {name}.chd  chd_0\n")
        f.write(f"  RCH6  {name}.rcha  rcha_0\n")
        f.write(f"  OC6  {name}.oc  oc\n")
        f.write("END packages\n")

    with open(test.workspace / f"{name}.dis", "w") as f:
        f.write("BEGIN options\n")
        f.write(f"  NCF6  FILEIN  {name}.dis.ncf\n")
        f.write("END options\n\n")
        f.write("BEGIN dimensions\n")
        f.write(f"  NLAY  {nlay}\n")
        f.write(f"  NROW  {nrow}\n")
        f.write(f"  NCOL  {ncol}\n")
        f.write("END dimensions\n\n")
        f.write("BEGIN griddata\n")
        f.write("  delr NETCDF\n")
        f.write("  delc NETCDF\n")
        f.write("  top NETCDF\n")
        f.write("  botm NETCDF\n")
        f.write("END griddata\n\n")

    with open(test.workspace / f"{name}.ic", "w") as f:
        f.write("BEGIN options\n")
        f.write("END options\n\n")
        f.write("BEGIN griddata\n")
        f.write("  strt NETCDF\n")
        f.write("END griddata\n")

    with open(test.workspace / f"{name}.npf", "w") as f:
        f.write("BEGIN options\n")
        f.write("END options\n\n")
        f.write("BEGIN griddata\n")
        f.write("  icelltype NETCDF\n")
        f.write("  k NETCDF\n")
        f.write("END griddata\n")

    with open(test.workspace / f"{name}.rcha", "w") as f:
        f.write("BEGIN options\n")
        f.write("  READASARRAYS\n")
        f.write("  auxiliary  RCH_MULT\n")
        f.write("  AUXMULTNAME  RCH_MULT\n")
        f.write("END options\n\n")
        # period 1: both recharge and aux
        f.write("BEGIN period 1\n")
        f.write("  recharge NETCDF\n")
        f.write("  RCH_MULT NETCDF\n")
        f.write("END period 1\n\n")
        # period 2: only recharge (aux should persist)
        f.write("BEGIN period 2\n")
        f.write("  recharge NETCDF\n")
        f.write("END period 2\n\n")
        # period 3: nothing (both should persist from period 2)

    success, buff = flopy.run_model(
        test.targets["mf6"],
        test.workspace / "mfsim.nam",
        model_ws=test.workspace,
        report=True,
    )
    assert success, "Pass 2 run with NETCDF recharge failed"

    # verify recharge budget from pass 2
    cbc = flopy.utils.CellBudgetFile(
        str(test.workspace / f"{name}.cbc"), precision="double"
    )
    rch_budgets = cbc.get_data(text="RCH")

    expected_q = [
        rch_vals[0] * mult_val * 100.0 * 100.0,
        rch_vals[1] * mult_val * 100.0 * 100.0,
        rch_vals[2] * mult_val * 100.0 * 100.0,
    ]

    for iper, (rch_data, exp_flux) in enumerate(zip(rch_budgets, expected_q)):
        q = rch_data["q"]
        actual_max = q.max()
        msg = f"Period {iper + 1}: expected rch flux {exp_flux}, got {actual_max}"
        assert np.isclose(actual_max, exp_flux, rtol=1e-6), msg
        if iper > 0:
            assert actual_max > 0.0, (
                f"Period {iper + 1}: recharge is zero - "
                f"AUXMULTNAME not persisting across periods (NETCDF path)"
            )


@pytest.mark.netcdf
@pytest.mark.developmode
def test_mf6model(function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(0, t),
        check=lambda t: check_output(0, t),
        cargs=["--mode=validate"],
    )
    test.run()
