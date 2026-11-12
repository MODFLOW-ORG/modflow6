"""
Test the API package in a GWT model.

Two versions of a one-dimensional transport model are run. In the first, run
with mf6, a mass source is specified with the SRC package and first-order
decay with the MST package. In the second, run through libmf6, the same two
terms are supplied from outside the program through an API package: the mass
source goes into the package's RHS and the decay coefficient into its HCOF.
The concentrations from the two runs are compared.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwt_api01"]

# API package name
api_pname = "API-1"

# model spatial dimensions
nlay, nrow, ncol = 1, 1, 10
delr = delc = 1.0
top, botm = 1.0, 0.0

# temporal discretization
tdis_rc = [(10.0, 10, 1.0)]

# flow
hk = 1.0
chd_spd = {0: [[(0, 0, 0), 2.0], [(0, 0, ncol - 1), 1.0]]}

# transport
porosity = 0.2
decay = 0.1
source_cell = 2
mass_rate = 1.0

# coefficient the MST package adds to the diagonal for first-order decay in a
# saturated cell
decay_hcof = -decay * delr * delc * (top - botm) * porosity

# solver data
nouter, ninner = 100, 300
hclose, rclose = 1e-10, 1e-10


def get_model(ws, name, api=False):
    gwfname, gwtname = f"gwf_{name}", f"gwt_{name}"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=tdis_rc)

    # flow model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        filename=f"{gwfname}.ims",
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
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL")],
    )

    # transport model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
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
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtssm(gwt)
    if api:
        flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
        flopy.mf6.ModflowGwtapi(gwt, maxbound=ncol, pname=api_pname)
    else:
        flopy.mf6.ModflowGwtmst(
            gwt, porosity=porosity, first_order_decay=True, decay=decay
        )
        flopy.mf6.ModflowGwtsrc(
            gwt, stress_period_data={0: [[(0, 0, source_cell), mass_rate]]}
        )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    return sim


def build_models(idx, test):
    name = cases[idx]

    # model with the SRC package and MST decay
    sim = get_model(test.workspace, name)

    # model with the API package
    mc = get_model(os.path.join(test.workspace, "libmf6"), name, api=True)

    return sim, mc


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    gwtname = f"gwt_{cases[idx]}".upper()
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

    # the API package variables
    nbound = mf6.get_value_ptr(mf6.get_var_address("NBOUND", gwtname, api_pname))
    nodelist = mf6.get_value_ptr(mf6.get_var_address("NODELIST", gwtname, api_pname))
    hcof = mf6.get_value_ptr(mf6.get_var_address("HCOF", gwtname, api_pname))
    rhs = mf6.get_value_ptr(mf6.get_var_address("RHS", gwtname, api_pname))

    # put every cell in the package
    nbound[0] = ncol
    nodelist[:] = np.arange(1, ncol + 1)

    # maximum outer iterations of each solution
    nsolutions = mf6.get_subcomponent_count()
    max_iter = [
        mf6.get_value(mf6.get_var_address("MXITER", f"SLN_{isol}"))[0]
        for isol in range(1, nsolutions + 1)
    ]

    # model time loop
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)

        # decay in every cell and the mass source in one cell
        hcof[:] = decay_hcof
        rhs[:] = 0.0
        rhs[source_cell] = -mass_rate

        for isol in range(1, nsolutions + 1):
            mf6.prepare_solve(isol)
            kiter = 0
            has_converged = False
            while kiter < max_iter[isol - 1]:
                has_converged = mf6.solve(isol)
                kiter += 1
                if has_converged:
                    break
            if not has_converged:
                return False, open(output_file_path).readlines()
            mf6.finalize_solve(isol)

        # finalize time step and update time
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

    # cleanup
    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    # cleanup and return
    return True, open(output_file_path).readlines()


def check_output(idx, test):
    fname = f"gwt_{cases[idx]}.ucn"
    conc = []
    for ws in (test.workspace, os.path.join(test.workspace, "libmf6")):
        cobj = flopy.utils.HeadFile(
            os.path.join(ws, fname), precision="double", text="CONCENTRATION"
        )
        conc.append(cobj.get_alldata())
        cobj.close()
    conc_src, conc_api = conc

    # the source and the decay must both have had an effect
    assert conc_src.max() > 1.0, "the mass source did not raise the concentration"

    assert conc_api.shape == conc_src.shape
    assert np.allclose(conc_api, conc_src, rtol=0.0, atol=1e-8), (
        "concentrations with the API package differ from those with the SRC "
        f"package and MST decay by up to {np.abs(conc_api - conc_src).max()}"
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
