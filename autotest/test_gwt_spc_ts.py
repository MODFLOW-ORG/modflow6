"""SPC CONCENTRATION persistence, covering: a value linked to a time
series for BNDNO 1 continuing to track that series in a later period
whose own PERIOD block reappears (setting BNDNO 2) without repeating
BNDNO 1 (spccontinue); and a plain (non-TS) value for BNDNO 1 persisting
across a later, present-but-empty PERIOD block (spcempty).

SPC's BNDNO is period-relative -- tied to whatever order the flow
package's own list happens to be in for that period, not a stable
identity like IFNO or a grid node number. This test uses a WEL list
that stays identical and in the same order across all periods, so the
persistence check isn't confounded by BNDNO reordering.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["spccontinue", "spcempty"]

TS_VAL_P1 = 10.0
TS_VAL_P2 = 20.0
TS_VAL_P3 = 30.0
BNDNO2_LITERAL = 99.0
EMPTY_CONC_VAL = 10.0


def _build_base(name, ws):
    nlay, nrow, ncol = 1, 1, 5
    delr, delc = 1.0, 1.0
    top, botm = 0.0, [-1.0]
    idomain = np.ones((nlay, nrow, ncol), dtype=int)

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=3, perioddata=[(1.0, 5, 1.0)] * 3
    )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)
    flopy.mf6.ModflowIms(sim, print_option="NONE", linear_acceleration="BICGSTAB")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=10.0, k33=10.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 0.0], [(0, 0, ncol - 1), 0.0]],
        pname="CHD-1",
    )
    # 2 wells, same order, same rate, defined once, unchanged across all
    # periods -- avoids the BNDNO-reordering complication entirely
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=[[(0, 0, 1), 1.0], [(0, 0, 3), 1.0]],
        pname="WEL-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim, model_type="gwt6", modelname=gwtname, model_nam_file=f"{gwtname}.nam"
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(ims, [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30)

    return sim, gwf, gwt, gwfname, gwtname


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim, gwf, gwt, gwfname, gwtname = _build_base(name, ws)

    if name == "spccontinue":
        spc = flopy.mf6.ModflowUtlspc(
            gwt,
            print_input=True,
            maxbound=2,
            perioddata={
                0: [(0, "concentration", "conc_ts")],
                # period 1: SPC's own PERIOD block reappears (setting
                # BNDNO 2, literal) without repeating BNDNO 1's TS-linked
                # setting
                1: [(1, "concentration", BNDNO2_LITERAL)],
                2: [(1, "concentration", BNDNO2_LITERAL)],
            },
            filename=f"{gwtname}.wel1.spc",
        )
        ts_data = [
            (0.0, TS_VAL_P1),
            (1.0, TS_VAL_P2),
            (2.0, TS_VAL_P3),
            (3.0, TS_VAL_P3),
        ]
        spc.ts.initialize(
            filename="conc.ts",
            timeseries=ts_data,
            time_series_namerecord="conc_ts",
            interpolation_methodrecord="stepwise",
        )
    else:  # spcempty
        flopy.mf6.ModflowUtlspc(
            gwt,
            print_input=True,
            maxbound=2,
            # periods 1, 2: SPC's own PERIOD block reappears present but
            # empty, without repeating BNDNO 1's plain (non-TS) setting
            perioddata={
                0: [(0, "concentration", EMPTY_CONC_VAL)],
                1: [],
                2: [],
            },
            filename=f"{gwtname}.wel1.spc",
        )

    flopy.mf6.ModflowGwtssm(
        gwt,
        print_flows=True,
        sources=[()],
        fileinput=[("WEL-1", f"{gwtname}.wel1.spc")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name
    lst_fname = test.workspace / f"{gwtname}.lst"
    lines = lst_fname.read_text().splitlines()

    # parse "INPUT VALUES FOR CONCENTRATION PACKAGE" blocks; each occurrence
    # is one apply_input_values() call (once per timestep once ts_active).
    # capture BNDNO 1's value from each block, in order.
    bndno1_vals = []
    i = 0
    while i < len(lines):
        if "INPUT VALUES FOR CONCENTRATION PACKAGE" in lines[i]:
            # next line is the "NO.  CONCENTRATION" header; then rows
            j = i + 2
            while j < len(lines):
                m = re.match(r"\s*(\d+)\s+([\d.eE+-]+)\s*$", lines[j])
                if not m:
                    break
                if int(m.group(1)) == 1:
                    bndno1_vals.append(float(m.group(2)))
                j += 1
            i = j
        else:
            i += 1

    print(f"SPC BNDNO 1 concentration values, in order: {bndno1_vals}")

    if name == "spccontinue":
        # TS-linked: apply_input_values() re-evaluates every timestep
        assert len(bndno1_vals) == 15, (
            f"Expected 15 applied-value blocks (5/period, 3 periods), "
            f"got {len(bndno1_vals)}: {bndno1_vals}"
        )
        period1, period2, period3 = (
            bndno1_vals[:5],
            bndno1_vals[5:10],
            bndno1_vals[10:],
        )
        assert np.allclose(period1, TS_VAL_P1), (
            f"Period 1 BNDNO 1 concentration expected {TS_VAL_P1}, got {period1}"
        )
        assert np.allclose(period2, TS_VAL_P2), (
            f"Period 2 BNDNO 1 concentration expected {TS_VAL_P2} (TS should "
            f"still be tracked even though period 2's SPC PERIOD block "
            f"reappears for BNDNO 2 without repeating BNDNO 1), got {period2}"
        )
        assert np.allclose(period3, TS_VAL_P3), (
            f"Period 3 BNDNO 1 concentration expected {TS_VAL_P3}, got {period3}"
        )
    else:  # spcempty
        # plain value: apply_input_values() only echoes once per period,
        # not per timestep, since there's no series to re-evaluate
        assert len(bndno1_vals) == 3, (
            f"Expected 3 applied-value blocks (1/period, 3 periods), "
            f"got {len(bndno1_vals)}: {bndno1_vals}"
        )
        assert np.allclose(bndno1_vals, EMPTY_CONC_VAL), (
            f"BNDNO 1 concentration expected {EMPTY_CONC_VAL} in every "
            f"period -- plain (non-TS) value did not persist across "
            f"periods 2/3's present-but-empty PERIOD blocks, "
            f"got {bndno1_vals}"
        )


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
