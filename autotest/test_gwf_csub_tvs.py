"""
Test time-varying storage (TVS) changes to specific storage with the CSUB
package. CSUB requires specific storage in the STO package to be zero, so
a TVS change of SS to a non-zero value is an input error.

Transient flow in a row of confined cells with CSUB coarse-grained
storage and STO specific storage zero.

Case 1 (csub-tvs01): TVS sets SS to a non-zero value at the start of
    period 2. Input error.
Case 2 (csub-tvs02): TVS sets SS from a time series that becomes non-zero
    at the second time step of period 2. Input error.
Case 3 (csub-tvs03): TVS sets SS to zero at the start of period 2.
Case 4 (csub-tvs04): convertible cells; TVS changes SY at the start of
    period 2 (CSUB does not use SY).
"""

import flopy
import pytest
from framework import TestFramework

cases = ["csub-tvs01", "csub-tvs02", "csub-tvs03", "csub-tvs04"]
xfail = {"csub-tvs01", "csub-tvs02"}
errors = ["cannot change SS for cell (1,1,3) to", "because the CSUB package is active."]

nlay, nrow, ncol = 1, 1, 5
cell = (0, 0, 2)


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 2, 1.0), (1.0, 2, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowIms(sim, inner_dvclose=1e-10, outer_dvclose=1e-10)
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    iconvert, sy = (1, 0.1) if name == "csub-tvs04" else (0, 0.0)
    sto = flopy.mf6.ModflowGwfsto(
        gwf, iconvert=iconvert, ss=0.0, sy=sy, transient={0: True}
    )
    if name == "csub-tvs02":
        tvs = flopy.mf6.ModflowUtltvs(
            sto, perioddata={1: [(cell, "SS", "sts")]}, filename=f"{name}.tvs"
        )
        tvs.ts.initialize(
            filename=f"{name}.tvs.ts",
            timeseries=[(0.0, 0.0), (1.5, 1e-3), (2.0, 1e-3)],
            time_series_namerecord="sts",
            interpolation_methodrecord="stepwise",
        )
    elif name == "csub-tvs04":
        flopy.mf6.ModflowUtltvs(
            sto, perioddata={1: [(cell, "SY", 0.2)]}, filename=f"{name}.tvs"
        )
    else:
        ss = 1e-3 if name == "csub-tvs01" else 0.0
        flopy.mf6.ModflowUtltvs(
            sto, perioddata={1: [(cell, "SS", ss)]}, filename=f"{name}.tvs"
        )
    flopy.mf6.ModflowGwfcsub(
        gwf, ninterbeds=0, cg_theta=0.3, cg_ske_cr=1e-4, packagedata=None
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={
            0: [[(0, 0, 0), 5.0], [(0, 0, ncol - 1), 5.0]],
            1: [[(0, 0, 0), 8.0], [(0, 0, ncol - 1), 5.0]],
        },
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    if name in xfail:
        text = " ".join((test.workspace / "mfsim.lst").read_text().split())
        for error in errors:
            assert error in text, f"expected error not reported: {error}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=name in xfail,
    )
    test.run()
