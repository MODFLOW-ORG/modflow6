"""
Test the API package in a GWE model. A linear exchange term,
k * (text - t), is applied through the API package (hcof = -k,
rhs = -k * text) to a subset of cells. The model has no advection or
conduction and no flow, so each cell is isolated and the temperature
follows the backward Euler recursion

    t1 = (s * t0 - rhs) / (s - hcof),  s = heat capacity / delt

Cells without an API boundary must keep their initial temperature.
The API budget term is checked against hcof * t - rhs.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwe_api01"]

# temporal discretization
nper = 1
perlen = 10.0
nstp = 10
delt = perlen / nstp

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 10
delr = delc = 1.0
top, botm = 1.0, 0.0
cell_volume = delr * delc * (top - botm)

# energy transport
porosity = 0.2
rhow, cpw = 1000.0, 4184.0
rhos, cps = 2500.0, 800.0
heat_capacity = (porosity * rhow * cpw + (1.0 - porosity) * rhos * cps) * cell_volume
t0 = 10.0

# api boundaries: non-contiguous cells, varying k and external temperature
api_cells = np.array([1, 4, 5, 8])  # zero-based columns
k_exchange = 0.2 * heat_capacity * np.arange(1, len(api_cells) + 1)
text = np.array([20.0, 0.0, 35.0, 15.0])
api_packname = "MYAPI"

# solver data
nouter, ninner = 100, 300
hclose, rclose = 1e-9, 1e-6


def build_models(idx, test):
    name = cases[idx]
    gwfname = "gwf"
    gwename = name
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=test.workspace, exe_name="libmf6"
    )
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, nstp, 1.0)])

    # flow model, uniform head so there is no flow
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        filename=f"{gwfname}.ims",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
    )
    sim.register_ims_package(imsgwf, [gwfname])
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
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary="TEMPERATURE",
        stress_period_data=[[(0, 0, 0), 1.0, t0]],
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )

    # energy transport model, no advection or conduction
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, save_flows=True)
    imsgwe = flopy.mf6.ModflowIms(
        sim,
        filename=f"{gwename}.ims",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(imsgwe, [gwename])
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGweic(gwe, strt=t0)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=porosity,
        density_water=rhow,
        heat_capacity_water=cpw,
        density_solid=rhos,
        heat_capacity_solid=cps,
    )
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "TEMPERATURE")])
    flopy.mf6.ModflowGweapi(
        gwe, maxbound=ncol, pname=api_packname, print_flows=True, save_flows=True
    )
    flopy.mf6.ModflowGweoc(
        gwe,
        temperature_filerecord=f"{gwename}.ucn",
        budget_filerecord=f"{gwename}.cbc",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )

    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    return sim


def check_output(idx, test):
    name = cases[idx]
    ws = test.workspace
    temp = flopy.utils.HeadFile(ws / f"{name}.ucn", text="TEMPERATURE").get_alldata()
    temp = temp.reshape(nstp, ncol)

    # expected temperatures
    hcof = -k_exchange
    rhs = -k_exchange * text
    s = heat_capacity / delt
    expected = np.full((nstp, ncol), t0)
    tcur = np.full(len(api_cells), t0)
    for kstp in range(nstp):
        tcur = (s * tcur - rhs) / (s - hcof)
        expected[kstp, api_cells] = tcur

    assert np.allclose(temp, expected, atol=1e-6), (
        "API package temperatures do not match the expected solution, "
        f"max difference = {np.abs(temp - expected).max()}"
    )

    # API budget term
    cbc = flopy.utils.CellBudgetFile(ws / f"{name}.cbc", precision="double")
    api_flows = cbc.get_data(text="API")
    assert len(api_flows) == nstp
    for kstp, rec in enumerate(api_flows):
        assert np.array_equal(rec["node"], api_cells + 1)
        q = hcof * temp[kstp, api_cells] - rhs
        assert np.allclose(rec["q"], q, rtol=1e-6), (
            f"API budget term does not match hcof * t - rhs in step {kstp + 1}"
        )


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

    # initialize the model
    try:
        mf6.initialize()
    except:
        return False, open(output_file_path).readlines()

    # get pointers to API package data
    nbound = mf6.get_value_ptr(mf6.get_var_address("NBOUND", name, api_packname))
    nodelist = mf6.get_value_ptr(mf6.get_var_address("NODELIST", name, api_packname))
    hcof = mf6.get_value_ptr(mf6.get_var_address("HCOF", name, api_packname))
    rhs = mf6.get_value_ptr(mf6.get_var_address("RHS", name, api_packname))

    # apply the exchange term to the api cells
    nb = len(api_cells)
    nbound[0] = nb
    nodelist[:nb] = api_cells + 1

    # model time loop
    nsol = mf6.get_subcomponent_count()
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        mf6.prepare_time_step(mf6.get_time_step())
        for sol_id in range(1, nsol + 1):
            max_iter = mf6.get_value(mf6.get_var_address("MXITER", f"SLN_{sol_id}"))
            mf6.prepare_solve(sol_id)
            for _ in range(max_iter[0]):
                hcof[:nb] = -k_exchange
                rhs[:nb] = -k_exchange * text
                has_converged = mf6.solve(sol_id)
                if has_converged:
                    break
            mf6.finalize_solve(sol_id)
            if not has_converged:
                print(f"solution {sol_id} did not converge")
                break
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()
        if not has_converged:
            break

    # cleanup
    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        api_func=lambda exe, ws: api_func(exe, idx, ws),
        compare=None,
    )
    test.run()
