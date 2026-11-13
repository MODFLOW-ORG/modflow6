"""
Test the API package in a GWT model. A linear exchange term,
k * (cext - c), is applied to every cell through the API package
(hcof = -k, rhs = -k * cext). The comparison simulation, run without
the API, represents the same term with first-order decay in MST
(k = decay * porosity * cell volume) and an SRC mass rate of k * cext.
Regression test for https://github.com/MODFLOW-ORG/modflow6/issues/3154
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwt_api01"]

# temporal discretization
nper = 1
perlen = 10.0
nstp = 10

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 10
delr = delc = 1.0
top, botm = 1.0, 0.0
cell_volume = delr * delc * (top - botm)

# flow
hk = 1.0
chd_spd = [[(0, 0, 0), 1.0, 1.0], [(0, 0, ncol - 1), 0.0, 0.0]]

# transport
porosity = 0.2
decay = 0.5
cext = 2.0
k_exchange = decay * porosity * cell_volume
api_packname = "MYAPI"

# solver data
nouter, ninner = 100, 300
hclose, rclose = 1e-9, 1e-6


def get_model(ws, name, api=False):
    gwfname = "gwf"
    gwtname = name
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, nstp, 1.0)])

    # flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
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
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, k=hk)
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary="CONCENTRATION",
        stress_period_data=chd_spd,
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # transport model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        filename=f"{gwtname}.ims",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    if api:
        flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
        flopy.mf6.ModflowGwtapi(
            gwt, maxbound=ncol, pname=api_packname, print_flows=True
        )
    else:
        flopy.mf6.ModflowGwtmst(
            gwt, porosity=porosity, first_order_decay=True, decay=decay
        )
        src_spd = [[(0, 0, j), k_exchange * cext] for j in range(ncol)]
        flopy.mf6.ModflowGwtsrc(gwt, stress_period_data=src_spd)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        budget_filerecord=f"{gwtname}.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sim = get_model(test.workspace, name)
    mc = get_model(os.path.join(test.workspace, "libmf6"), name, api=True)
    return sim, mc


def check_output(idx, test):
    name = cases[idx]
    fpth = test.workspace / f"{name}.ucn"
    conc = flopy.utils.HeadFile(fpth, text="CONCENTRATION").get_alldata()
    fpth = test.workspace / "libmf6" / f"{name}.ucn"
    conc_api = flopy.utils.HeadFile(fpth, text="CONCENTRATION").get_alldata()
    assert conc.max() > 0.1, "exchange term had no effect on concentrations"
    assert np.allclose(conc, conc_api, atol=1e-6), (
        "API package concentrations do not match the decay + SRC solution, "
        f"max difference = {np.abs(conc - conc_api).max()}"
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

    # apply the exchange term to every cell
    nbound[0] = ncol
    nodelist[:] = np.arange(1, ncol + 1)

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
                hcof[:] = -k_exchange
                rhs[:] = -k_exchange * cext
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
    )
    test.run()
