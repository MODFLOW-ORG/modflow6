"""
Confirm the solver output reflects stress-period-varying IMS linear settings.

Two GWF models are coupled by a GWF-GWF exchange (the same layout used by the
parallel tests) so this case can be reused in parallel. Period 1 uses the base
LINEAR settings -- a loose INNER_DVCLOSE (1e-3) with INNER_MAXIMUM set to the
number of inner iterations that tolerance actually needs -- and a PERIOD 2 block
tightens INNER_DVCLOSE and raises INNER_MAXIMUM.

The constant heads change in period 2 so the tighter settings are exercised by a
fresh solve. Without that change period 2 would start from the period 1 solution,
which the outer loop has already converged to OUTER_DVCLOSE, leaving the inner
solve too little work to distinguish the two tolerances.

The inner/outer iteration data are written to both mfsim.lst (PRINT_OPTION ALL)
and the IMS csv files, and the check asserts the recorded iteration counts and
dependent-variable changes actually track the per-period settings.

flopy does not support the option yet, so LINEAR_PERIODDATA and the period-data
file are injected into the written input.
"""

import glob
import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ims_period_out"]

name_left = "leftmodel"
name_right = "rightmodel"

nper = 2
nlay, nrow, ncol = 1, 10, 10
delr = delc = 100.0

# constant heads by stress period. the right-hand boundary is raised in period 2
# so that period has a fresh system to solve
h_left = 1.0
h_right = [10.0, 20.0]

# solver data. rclose is deliberately non-binding so that INNER_DVCLOSE is the
# criterion controlling the inner iterations in both periods
nouter, rclose = 200, 1.0e3
outer_dvclose = 1.0e-9

# period 1 (base LINEAR settings): loose tolerance, and an inner maximum equal to
# the number of inner iterations that tolerance needs to converge
dvclose_loose, inner_max_loose = 1.0e-3, 11
# period 2 (PERIOD block): tight tolerance and a larger inner maximum. 1e-10 is
# tight enough to need many more inner iterations but stays above the level where
# BiCGSTAB stagnates on this system
dvclose_tight, inner_max_tight = 1.0e-10, 500


def inject_period_data(ws, period_text):
    """Add LINEAR_PERIODDATA to the IMS options and write the period-data file."""
    # flopy writes "BEGIN options" (mixed case), so match case-insensitively
    ims_files = glob.glob(os.path.join(ws, "*.ims"))
    assert len(ims_files) == 1, f"expected one .ims file, found {ims_files}"
    ims = ims_files[0]
    impd = os.path.splitext(os.path.basename(ims))[0] + ".impd"
    with open(ims) as f:
        txt = f.read()
    new = re.sub(
        r"(?i)(BEGIN\s+OPTIONS[^\n]*\n)",
        rf"\1  LINEAR_PERIODDATA FILEIN {impd}\n",
        txt,
        count=1,
    )
    assert new != txt, "could not inject LINEAR_PERIODDATA into the IMS OPTIONS block"
    with open(ims, "w") as f:
        f.write(new)
    with open(os.path.join(ws, impd), "w") as f:
        f.write(period_text)


def get_model(name, ws):
    """Build two GWF models coupled by a GWF-GWF exchange."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(1.0, 1, 1.0)] * nper
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=outer_dvclose,
        outer_maximum=nouter,
        inner_maximum=inner_max_loose,
        inner_dvclose=dvclose_loose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        csv_outer_output_filerecord=f"{name}.outer.csv",
        csv_inner_output_filerecord=f"{name}.inner.csv",
        pname="ims",
    )

    shift_x = ncol * delr
    for modelname, chd_col, xorigin in (
        (name_left, 0, 0.0),
        (name_right, ncol - 1, shift_x),
    ):
        gwf = flopy.mf6.ModflowGwf(sim, modelname=modelname, save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            xorigin=xorigin,
            top=0.0,
            botm=[-100.0],
        )
        flopy.mf6.ModflowGwfic(gwf, strt=0.0)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
        if modelname == name_left:
            chd_spd = {0: [[(0, irow, chd_col), h_left] for irow in range(nrow)]}
        else:
            chd_spd = {
                kper: [[(0, irow, chd_col), h] for irow in range(nrow)]
                for kper, h in enumerate(h_right)
            }
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=f"{modelname}.hds",
            saverecord=[("HEAD", "ALL")],
        )

    angldegx, cdist = 0.0, delr
    gwfgwf_data = [
        [
            (0, irow, ncol - 1),
            (0, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for irow in range(nrow)
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=name_left,
        exgmnameb=name_right,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
    )
    return sim


def build_models(idx, test):
    sim = get_model(cases[idx], test.workspace)
    sim.write_simulation()
    inject_period_data(
        test.workspace,
        f"BEGIN PERIOD 2\n"
        f"  INNER_MAXIMUM {inner_max_tight}\n"
        f"  INNER_DVCLOSE {dvclose_tight}\n"
        f"END PERIOD\n",
    )
    return sim, None


def _rank_files(ws, serial, parallel):
    """Return the serial file, or the per-rank files written by a parallel run."""
    paths = sorted(glob.glob(os.path.join(ws, parallel)))
    if not paths:
        paths = [os.path.join(ws, serial)]
    for path in paths:
        assert os.path.isfile(path), f"{path} was not written"
    return paths


def _period_block_settings(lst, kper, path):
    """Parse the settings reported in a PROCESSING LINEAR PERIOD DATA block."""
    match = re.search(
        rf"PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD {kper}\n(.*?)\n"
        rf"\s*END PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD {kper}",
        lst,
        re.DOTALL,
    )
    assert match, (
        f"the period {kper} block in {path} reports no settings; the listing "
        f"would keep advertising the base values after the period changed them"
    )
    reported = {}
    for line in match.group(1).splitlines():
        key, sep, value = line.partition("=")
        if sep:
            reported[key.strip()] = value.strip()
    return reported


def _read_inner_csv(path):
    """Group the inner-iteration csv rows by stress period and outer iteration."""
    import csv

    with open(path) as f:
        rows = list(csv.DictReader(f))
    assert rows, f"{path} is empty"
    per = {}
    for row in rows:
        kper, nouter = int(row["kper"]), int(row["nouter"])
        per.setdefault(kper, {}).setdefault(nouter, []).append(row)
    return per


def check_output(idx, test):
    ws = test.workspace
    name = cases[idx]

    # the period data must have been read and applied, on every rank
    for path in _rank_files(ws, "mfsim.lst", "mfsim.p*.lst"):
        with open(path) as f:
            lst = f.read()
        assert "LINEAR_PERIODDATA input will be read" in lst, (
            f"the period-data file was not opened ({path})"
        )
        assert "PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD 2" in lst, (
            f"the period 2 linear settings were not processed ({path})"
        )

        # the updated settings are echoed, so the listing reports the values that
        # are actually in effect rather than only the base LINEAR block values
        reported = _period_block_settings(lst, 2, path)
        assert int(reported["INNER_MAXIMUM"]) == inner_max_tight, (
            f"the listing reports INNER_MAXIMUM {reported['INNER_MAXIMUM']} for "
            f"period 2, expected {inner_max_tight} ({path})"
        )
        assert float(reported["INNER_DVCLOSE"]) == pytest.approx(dvclose_tight), (
            f"the listing reports INNER_DVCLOSE {reported['INNER_DVCLOSE']} for "
            f"period 2, expected {dvclose_tight} ({path})"
        )

        # PRINT_OPTION ALL writes an inner and outer summary for each period
        for tag in ("OUTER ITERATION SUMMARY", "INNER ITERATION SUMMARY"):
            assert lst.count(tag) == nper, (
                f"expected {nper} '{tag}' tables in {path}, found {lst.count(tag)}"
            )

    # the outer csv is written alongside the inner csv
    _rank_files(ws, f"{name}.outer.csv", f"{name}.outer.p*.csv")

    for path in _rank_files(ws, f"{name}.inner.csv", f"{name}.inner.p*.csv"):
        per = _read_inner_csv(path)
        assert sorted(per) == list(range(1, nper + 1)), (
            f"expected inner csv data for periods 1-{nper}, found {sorted(per)} "
            f"({path})"
        )

        ninner1 = max(len(rows) for rows in per[1].values())
        ninner2 = max(len(rows) for rows in per[2].values())

        # period 1 honors the base INNER_MAXIMUM, and period 2 exceeds it, which
        # is only possible if the raised PERIOD 2 INNER_MAXIMUM was applied
        assert ninner1 <= inner_max_loose, (
            f"period 1 used {ninner1} inner iterations, exceeding INNER_MAXIMUM "
            f"{inner_max_loose} ({path})"
        )
        assert ninner2 > inner_max_loose, (
            f"period 2 used at most {ninner2} inner iterations, so the PERIOD 2 "
            f"INNER_MAXIMUM {inner_max_tight} was not applied ({path})"
        )

        # the loose period stops short of the tight tolerance and the tight period
        # drives the dependent-variable change down to the period 2 INNER_DVCLOSE,
        # so the tolerance in effect tracks the stress period
        dvmax1 = abs(float(per[1][1][-1]["solution_inner_dvmax"]))
        dvmax2 = abs(float(per[2][1][-1]["solution_inner_dvmax"]))
        assert dvmax1 > 1.0e-6, (
            f"period 1 reached dvmax {dvmax1:.3e}, which is tighter than the loose "
            f"INNER_DVCLOSE {dvclose_loose} should allow ({path})"
        )
        assert dvmax2 <= dvclose_tight, (
            f"period 2 reached dvmax {dvmax2:.3e}, which does not meet the tight "
            f"INNER_DVCLOSE {dvclose_tight} ({path})"
        )

        # the loose linear solve needs more Picard iterations than the tight one
        assert len(per[1]) > len(per[2]), (
            f"period 1 used {len(per[1])} outer iterations and period 2 used "
            f"{len(per[2])}; the loose inner tolerance should require more ({path})"
        )

    # both periods still produce the expected uniform flow field
    heads = [
        flopy.utils.HeadFile(os.path.join(ws, f"{m}.hds")).get_alldata()
        for m in (name_left, name_right)
    ]
    ncol_tot = 2 * ncol
    for kper in range(nper):
        simulated = np.concatenate(
            [h[kper].flatten()[0:ncol] for h in heads],
        )
        expected = h_left + (h_right[kper] - h_left) * np.arange(ncol_tot) / (
            ncol_tot - 1
        )
        np.testing.assert_allclose(simulated, expected, atol=1.0e-6)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        overwrite=False,  # keep the injected LINEAR_PERIODDATA + period file
    )
    test.run()
