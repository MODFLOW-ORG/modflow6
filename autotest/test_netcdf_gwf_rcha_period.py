"""
NetCDF structured export test for RCHA READASARRAYS with time-varying
per-period recharge arrays.

The existing rch01 and rch03 netcdf tests use a constant recharge value that
is the same across all stress periods. This test exercises the time-varying
case by assigning a distinct recharge array to each of three stress periods,
covering the full NCPL timeseries write path in the structured NetCDF export.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

xa = pytest.importorskip("xarray")
nc = pytest.importorskip("netCDF4")

# minimal DIS model; small enough to run fast, large enough to be realistic
nlay, nrow, ncol = 1, 5, 5
# three periods with distinct recharge values so all three slices must be written
rch_vals = [0.0001, 0.0002, 0.0003]
nper = len(rch_vals)
name = "rcha_period"


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

    # fixed heads on left and right columns; recharge varies by period
    chd_spd = {
        i: [(0, j, 0, 10.0) for j in range(nrow)]
        + [(0, j, ncol - 1, 9.0) for j in range(nrow)]
        for i in range(nper)
    }
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    rch_spd = {i: np.full((nrow, ncol), v) for i, v in enumerate(rch_vals)}
    flopy.mf6.ModflowGwfrcha(gwf, readasarrays=True, recharge=rch_spd)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )

    # mark packages for NetCDF input export (written during validate mode)
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
        assert "rcha_0_recharge" in xds, "rcha_0_recharge missing from NetCDF input"
        rch_var = xds["rcha_0_recharge"]
        assert rch_var.shape[0] == nper, (
            f"Expected {nper} time slices, got {rch_var.shape[0]}"
        )
        for i, v in enumerate(rch_vals):
            assert np.allclose(rch_var[i].values, v), (
                f"Period {i + 1} recharge mismatch: expected {v}, "
                f"got {rch_var[i].values}"
            )

    # --- re-run with NETCDF recharge for all periods (Pass 2) ---
    with open(test.workspace / f"{name}.nam", "w") as f:
        f.write("BEGIN options\n")
        f.write("  SAVE_FLOWS\n")
        f.write(f"  NETCDF_STRUCTURED  FILEOUT  {name}.nc\n")
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
        f.write("END options\n\n")
        for i in range(nper):
            f.write(f"BEGIN period {i + 1}\n")
            f.write("  recharge NETCDF\n")
            f.write(f"END period {i + 1}\n\n")

    success, buff = flopy.run_model(
        test.targets["mf6"],
        test.workspace / "mfsim.nam",
        model_ws=test.workspace,
        report=True,
    )
    assert success, "Pass 2 run with NETCDF recharge failed"
    test.success = success

    # verify head output from Pass 2 matches a fresh binary headfile
    hds_fpth = test.workspace / f"{name}.hds"
    hds = flopy.utils.HeadFile(str(hds_fpth), precision="double")
    with xa.open_dataset(test.workspace / f"{name}.nc") as xds:
        for kper in range(nper):
            rec = hds.get_data(kstpkper=(0, kper))
            assert np.allclose(np.array(rec), xds["head"][kper].values), (
                f"Head mismatch at period {kper + 1}"
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
