"""
Tests Eisenstat-Walker forcing terms for the linear solve and backtracking.

Uses the bump problem with backtracking.

Cases:
  - ew_nr  : Newton-Raphson with FORCING_TERMS; heads match the reference
             solution and the outer CSV output reports the forcing terms.
  - ew_sd  : Newton-Raphson with BACKTRACKING_SUFFICIENT_DECREASE; heads match the
             reference solution.
  - ew_std : standard formulation with FORCING_TERMS, which has no effect and
             is reported as a warning.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = ["ew_nr", "ew_sd", "ew_std"]
data_path = project_root_path / "autotest/data/ex-gwf-bump/"
nlay, nrow, ncol = 1, 51, 51
xlen, ylen = 100.0, 100.0
top = 25.0
k11 = 1.0
H1, H2 = 7.5, 2.5
delr = xlen / float(ncol)
delc = ylen / float(nrow)
shape3d = (nlay, nrow, ncol)
nouter = 200
ninner = 100
hclose = 1e-9
hclose_outer = hclose * 10.0
rclose = 1e-3
botm = np.loadtxt(data_path / "bottom.txt").reshape(shape3d)
chd_spd = [[0, i, 0, H1] for i in range(nrow)]
chd_spd += [[0, i, ncol - 1, H2] for i in range(nrow)]

FORCING_MSG = "Eisenstat-Walker forcing terms will be used"
WARNING_MSG = "BACKTRACKING_SUFFICIENT_DECREASE have no effect"


def build_models(idx, test):
    name = cases[idx]
    newton = name != "ew_std"

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=f"{name}.outer.csv",
        forcing_terms=True,
        backtracking_sufficient_decrease=name == "ew_sd",
        outer_maximum=nouter,
        outer_dvclose=hclose_outer,
        backtracking_number=20,
        backtracking_tolerance=1.05,
        backtracking_reduction_factor=0.1,
        backtracking_residual_limit=0.002,
        linear_acceleration="bicgstab",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="newton under_relaxation" if newton else None,
    )
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
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1 if newton else 0, k=k11)
    flopy.mf6.ModflowGwfic(gwf, strt=H1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    mfsim_lst = (test.workspace / "mfsim.lst").read_text()

    if name == "ew_std":
        assert WARNING_MSG in mfsim_lst, "missing FORCING_TERMS warning"
        assert FORCING_MSG not in mfsim_lst, "forcing terms used without Newton"
        df = pd.read_csv(test.workspace / f"{name}.outer.csv")
        assert "forcing_term" not in df.columns, "forcing columns without Newton"
        return

    assert WARNING_MSG not in mfsim_lst, "unexpected FORCING_TERMS warning"
    assert FORCING_MSG in mfsim_lst, "forcing terms not used"

    df = pd.read_csv(test.workspace / f"{name}.outer.csv")
    eta = df["forcing_term"].to_numpy()
    assert eta[0] == 0.5, "initial forcing term is not 0.5"
    assert np.all((eta >= 0.0) & (eta <= 0.9)), "forcing term out of range"
    assert np.all(df["residual_l2norm"] > 0.0), "missing residual norms"

    base_heads = flopy.utils.HeadFile(data_path / "results.hds.cmp").get_data()
    heads = flopy.utils.HeadFile(test.workspace / f"{name}.hds").get_data()
    assert np.allclose(base_heads, heads), "head comparison failed"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
