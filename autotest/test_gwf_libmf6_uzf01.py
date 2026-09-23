"""
Test that the pre-refactor uzf variable names resolve through the API.

The uzf variables were renamed, so the old names are checked in to the memory
manager as aliases of the renamed arrays. A program that addresses uzf
variables through the API must resolve either name, and both names must refer
to the same memory.

Cases:
  - libgwf_uzf01 : every aliased name resolves, matches the renamed variable it
                   aliases, and shares memory with it.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_uzf01"]

# aliases that did not resolve, filled in by api_func and checked afterward
alias_failures = []

# pre-refactor name : renamed variable it aliases
aliases = {
    "UZDPST": "WAVE_DEPTH",
    "UZTHST": "WAVE_THETA",
    "UZFLST": "WAVE_FLUX",
    "UZSPST": "WAVE_SPEED",
    "NWAVST": "NWAVES",
    "NWAV_PVAR": "NWAVES_MAX",
    "THTR": "THETA_RES",
    "THTS": "THETA_SAT",
    "THTI": "THETA_INIT",
    "EPS": "BC_EPS",
    "HA": "AIR_ENTRY",
    "HROOT": "ROOT_POT",
    "ROOTACT": "ROOT_ACT",
    "EXTWC": "THETA_EXT",
    "ETACT": "ET_UZ",
    "TOTFLUX": "FLUX_TO_WT",
    "SINF": "FINF_SPEC",
    "GWET_PVAR": "GWET",
    "PETMAX": "PET_MAX",
    "GWPET": "GW_PET",
    "EXTDP": "EXT_DEPTH",
    "EXTDPUZ": "EXT_DEPTH_UZ",
    "WATAB": "WATER_TABLE",
    "WATABOLD": "WATER_TABLE_OLD",
    "SURFLUX": "SURF_INFIL",
    "SURFLUXBELOW": "SURF_INFIL_BELOW",
    "SURFSEEP": "SURF_SEEP",
    "IVERTCON": "CELL_BELOW",
    "SINF_PVAR": "FINF_INPUT",
    "PET_PVAR": "PET_INPUT",
    "EXDP_PVAR": "EXTDP_INPUT",
    "EXTWC_PVAR": "EXTWC_INPUT",
    "HA_PVAR": "HA_INPUT",
    "HROOT_PVAR": "HROOT_INPUT",
    "ROOTACT_PVAR": "ROOTACT_INPUT",
    "NTRAIL_PVAR": "NTRAIL_INPUT",
    "NSETS": "NWAVESETS",
}

# model spatial dimensions
nlay, nrow, ncol = 3, 3, 3

# cell spacing
delr = delc = 10.0

# top and bottom of the aquifer
top = 10.0
botm = [6.0, 3.0, 0.0]

# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def get_model(ws, name):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=2.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(nlay - 1, 0, 0), 1.0]])

    # -- one uzf cell per model cell, connected vertically
    packagedata = []
    perioddata = []
    iuzno = 0
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                landflag = 1 if k == 0 else 0
                ivertcon = iuzno + nrow * ncol if k < nlay - 1 else -1
                packagedata.append(
                    [
                        iuzno,
                        (k, i, j),
                        landflag,
                        ivertcon,
                        0.05,
                        1.0,
                        0.05,
                        0.35,
                        0.05,
                        4.0,
                    ]
                )
                if landflag == 1:
                    perioddata.append([iuzno, 0.01, 0.005, 1.0, 0.1, 0.0, 0.0, 0.0])
                iuzno += 1
    flopy.mf6.ModflowGwfuzf(
        gwf,
        pname="uzf-1",
        nuzfcells=len(packagedata),
        ntrailwaves=7,
        nwavesets=40,
        simulate_et=True,
        unsat_etwc=True,
        linear_gwet=True,
        packagedata=packagedata,
        perioddata={0: perioddata},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sim = get_model(test.workspace, name)

    # -- the api_func is only invoked for a comparison model in a workspace
    #    named libmf6
    mc = get_model(os.path.join(test.workspace, "libmf6"), name)

    return sim, mc


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx].upper()
    if model_ws is None:
        model_ws = "."

    output_file_path = os.path.join(model_ws, "mfsim.stdout")

    try:
        mf6 = ModflowApi(exe, working_directory=model_ws)
    except Exception as e:
        print("Failed to load " + str(exe))
        print("with message: " + str(e))
        return False, open(output_file_path).readlines()

    try:
        mf6.initialize()
        mf6.update()
    except:
        return False, open(output_file_path).readlines()

    # -- an alias must resolve and refer to the same values as the renamed
    #    variable it aliases
    failures = alias_failures
    for old, new in aliases.items():
        try:
            vold = mf6.get_value_ptr(mf6.get_var_address(old, name, "UZF-1"))
            vnew = mf6.get_value_ptr(mf6.get_var_address(new, name, "UZF-1"))
        except Exception as e:
            failures.append(f"{old} -> {new}: {e}")
            continue
        if vold.shape != vnew.shape:
            failures.append(f"{old} -> {new}: shape {vold.shape} != {vnew.shape}")
        elif not np.array_equal(vold, vnew):
            failures.append(f"{old} -> {new}: values differ")

    # -- an alias must share memory with the renamed variable, not copy it
    etact = mf6.get_value_ptr(mf6.get_var_address("ETACT", name, "UZF-1"))
    et_uz = mf6.get_value_ptr(mf6.get_var_address("ET_UZ", name, "UZF-1"))
    saved = etact[0]
    etact[0] = -12345.0
    if et_uz[0] != -12345.0:
        failures.append("ETACT does not share memory with ET_UZ")
    etact[0] = saved

    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    return True, open(output_file_path).readlines()


def check_output(idx, test):
    assert not alias_failures, "aliased uzf variable names failed:\n  " + "\n  ".join(
        alias_failures
    )


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()
