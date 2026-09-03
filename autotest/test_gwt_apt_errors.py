"""
Error/edge-case coverage for TspAptType%apt_source_cvs / apt_source_options,
shared by the 8 GWT/GWE APT (advanced package transport) packages.  None of
these are exercised by any other regression test, since they are all
error paths:

  - apt_source_cvs must reject a PACKAGEDATA IFNO outside 1..ncv, and must
    reject a PACKAGEDATA block that doesn't cover every feature exactly
    once, instead of silently reading/writing out of bounds.
  - apt_source_cvs must source STRT by feature number (via PACKAGEDATA_IFNO),
    not by PACKAGEDATA row order, since rows may be supplied in any order.
  - find_*_package's "could not find flow package" error must be attributed
    to the correct input file.
  - a PERIOD block's keystring dispatch must reject a compound record's own
    sub-member name (e.g. AUXVAL, a sub-member of the AUXILIARY record) used
    directly as a top-level dispatch keyword, the same as any other
    unrecognized keyword, instead of silently accepting it.

Uses a minimal single-cell-row GWF+LAK / GWT+LKT model; LKT was chosen
arbitrarily among the 8 APT packages since the code under test is shared
base-class behavior (TspAptType), not package-specific.  The IFNO/missing-
feature/row-order cases are also duplicated onto GWE+LKE, since that code
path is shared but had only ever been exercised via a GWT package.
"""

import re
import subprocess

import flopy
import pytest


def run_mf6(argv, ws):
    buff = []
    proc = subprocess.Popen(
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, _ = proc.communicate()
    if result is not None:
        c = result.decode("utf-8").rstrip("\r\n")
        print(f"{c}")
        buff.append(c)
    return proc.returncode, buff


def run_mf6_error(ws, exe, err_str_list):
    returncode, buff = run_mf6([exe], ws)
    msg = "mf terminated with error"
    if returncode != 0:
        if not isinstance(err_str_list, list):
            err_str_list = [err_str_list]
        for err_str in err_str_list:
            if any(err_str in s for s in buff):
                raise RuntimeError(msg)
            else:
                msg += " but did not print correct error message."
                msg += f'  Correct message should have been "{err_str}"'
                raise ValueError(msg)
    else:
        raise ValueError("mf6 terminated successfully but was expected to fail")


def _build_gwf_lak(sim, name, nlakes=1):
    ncol = 5 if nlakes == 1 else 6
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam")
    flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", complexity="SIMPLE", filename=f"{gwfname}.ims"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=20.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), -0.5, 0.0],
            [(0, 0, ncol - 1), -0.5, 0.0],
        ],
        pname="CHD-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.chd",
    )
    connlen = connwidth = 0.5
    if nlakes == 1:
        packagedata = [(0, -0.4, 3, 0.0)]
        connectiondata = [
            (0, 0, (0, 0, 1), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", 0.0, 10, 10, connlen, connwidth),
        ]
        perioddata = [(0, "STATUS", "CONSTANT"), (0, "STAGE", -0.4)]
    else:
        packagedata = [(0, -0.4, 1, 0.0), (1, -0.4, 1, 0.0)]
        connectiondata = [
            (0, 0, (0, 0, 2), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (1, 0, (0, 0, 3), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
        ]
        perioddata = [
            (0, "STATUS", "CONSTANT"),
            (0, "STAGE", -0.4),
            (1, "STATUS", "CONSTANT"),
            (1, "STAGE", -0.4),
        ]
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=nlakes,
        noutlets=0,
        ntables=0,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="LAK-1",
        auxiliary=["CONCENTRATION"],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{gwfname}.cbc")
    return gwf, gwfname, ncol


def _build_gwt_lkt(
    sim,
    name,
    gwfname,
    ncol,
    flow_package_name="LAK-1",
    packagedata=None,
    lakeperioddata=None,
    auxiliary=None,
):
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, model_nam_file=f"{gwtname}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwtname}.ims"), [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=ncol, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(gwt, budget_filerecord=f"{gwtname}.cbc")
    if packagedata is None:
        packagedata = [(0, 35.0, 0.0, "mylake")] if auxiliary else [(0, 35.0, "mylake")]
    if lakeperioddata is None:
        lakeperioddata = [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 100.0)]
    lkt = flopy.mf6.ModflowGwtlkt(
        gwt,
        boundnames=True,
        auxiliary=auxiliary,
        packagedata=packagedata,
        lakeperioddata=lakeperioddata,
        flow_package_name=flow_package_name,
        pname="LKT-1",
        print_concentration=True,
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return gwt, gwtname


def test_lkt_flow_package_not_found(function_tmpdir, targets):
    """find_lkt_package must attribute its error to the LKT input file."""
    mf6 = targets["mf6"]
    name = "lkt_nopkg"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    _build_gwt_lkt(sim, name, gwfname, ncol, flow_package_name="nonexistent-pkg")
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "Could not find flow package with name NONEXISTENT-PKG",
        )


def test_lkt_packagedata_ifno_out_of_range(function_tmpdir, targets):
    """apt_source_cvs must reject an out-of-range PACKAGEDATA IFNO."""
    mf6 = targets["mf6"]
    name = "lkt_ifno"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    # only 1 lake exists (ncv=1), so IFNO=1 (0-based: 1) is out of range
    _build_gwt_lkt(sim, name, gwfname, ncol, packagedata=[(1, 35.0, "mylake")])
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "LKT PACKAGEDATA IFNO (2) must be greater than 0 and less than or equal",
        )


def test_lkt_packagedata_missing_feature(function_tmpdir, targets):
    """apt_source_cvs must reject a PACKAGEDATA block missing a feature."""
    mf6 = targets["mf6"]
    name = "lkt_missing"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    gwt, gwtname = _build_gwt_lkt(sim, name, gwfname, ncol)
    sim.write_simulation()

    # empty the PACKAGEDATA block (ncv=1, so feature 1 is now missing)
    lkt_fname = function_tmpdir / f"gwt_{name}.lkt"
    text = lkt_fname.read_text()
    lines = text.splitlines(keepends=True)
    out = []
    skipping = False
    for line in lines:
        if line.strip().upper().startswith("BEGIN PACKAGEDATA"):
            skipping = True
            out.append(line)
            continue
        if line.strip().upper().startswith("END PACKAGEDATA"):
            skipping = False
            out.append(line)
            continue
        if skipping:
            continue
        out.append(line)
    lkt_fname.write_text("".join(out))

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "LKT PACKAGEDATA no data specified for feature 1",
        )


def test_lkt_packagedata_out_of_order(function_tmpdir, targets):
    """apt_source_cvs must source STRT by feature number, not row order."""
    mf6 = targets["mf6"]
    name = "lkt_order"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0e-6, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name, nlakes=2)
    # rows deliberately out of IFNO order: row 1 -> feature 2, row 2 -> feature 1
    _build_gwt_lkt(
        sim,
        name,
        gwfname,
        ncol,
        packagedata=[(1, 222.0, "lake-two"), (0, 111.0, "lake-one")],
        lakeperioddata=[],
    )
    sim.write_simulation()

    returncode, _ = run_mf6([mf6], str(function_tmpdir))
    assert returncode == 0, "mf6 did not terminate successfully"

    lst_text = (function_tmpdir / f"gwt_{name}.lst").read_text()
    idx_one = lst_text.index("LAKE-ONE")
    idx_two = lst_text.index("LAKE-TWO")
    conc_one = float(lst_text[idx_one : idx_one + 60].split()[2])
    conc_two = float(lst_text[idx_two : idx_two + 60].split()[2])
    assert conc_one == pytest.approx(111.0, abs=0.5), (
        f"LAKE-ONE (feature 1) concentration {conc_one} does not match its "
        "own PACKAGEDATA STRT (111.0) -- STRT was likely sourced by row "
        "order instead of by feature number"
    )
    assert conc_two == pytest.approx(222.0, abs=0.5), (
        f"LAKE-TWO (feature 2) concentration {conc_two} does not match its "
        "own PACKAGEDATA STRT (222.0) -- STRT was likely sourced by row "
        "order instead of by feature number"
    )


def test_lkt_packagedata_aux_out_of_order(function_tmpdir, targets):
    """allocate_featureauxvar must source PACKAGEDATA AUX by feature
    number, not row order, the same as apt_source_cvs does for STRT."""
    mf6 = targets["mf6"]
    name = "lkt_auxorder"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0e-6, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name, nlakes=2)

    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, model_nam_file=f"{gwtname}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwtname}.ims"), [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=ncol, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt, budget_filerecord=f"{gwtname}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    # rows deliberately out of IFNO order: row 1 -> feature 2, row 2 -> feature 1
    flopy.mf6.ModflowGwtlkt(
        gwt,
        boundnames=True,
        save_flows=True,
        auxiliary=["myaux"],
        packagedata=[(1, 0.0, 999.0, "lake-two"), (0, 0.0, 111.0, "lake-one")],
        lakeperioddata=[],
        flow_package_name="LAK-1",
        pname="LKT-1",
        budget_filerecord=f"{gwtname}.lkt.bud",
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
    )
    sim.write_simulation()

    returncode, _ = run_mf6([mf6], str(function_tmpdir))
    assert returncode == 0, "mf6 did not terminate successfully"

    bobj = flopy.utils.CellBudgetFile(
        str(function_tmpdir / f"{gwtname}.lkt.bud"), precision="double"
    )
    rec = bobj.get_data(text="AUXILIARY")[-1]
    aux_one = rec["MYAUX"][0]
    aux_two = rec["MYAUX"][1]
    assert aux_one == pytest.approx(111.0, abs=0.5), (
        f"Feature 1 (lake-one) AUX {aux_one} does not match its own "
        "PACKAGEDATA AUX (111.0) -- AUX was likely sourced by row order "
        "instead of by feature number"
    )
    assert aux_two == pytest.approx(999.0, abs=0.5), (
        f"Feature 2 (lake-two) AUX {aux_two} does not match its own "
        "PACKAGEDATA AUX (999.0) -- AUX was likely sourced by row order "
        "instead of by feature number"
    )


def test_lkt_period_ifno_out_of_range(function_tmpdir, targets):
    """apt_rp must reject an out-of-range PERIOD IFNO instead of silently
    skipping the row.
    """
    mf6 = targets["mf6"]
    name = "lkt_pifno"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    # only 1 lake exists (ncv=1); PACKAGEDATA is valid, but the PERIOD row
    # below targets feature 5 (0-based), i.e. IFNO=6, which is out of range
    _build_gwt_lkt(sim, name, gwfname, ncol, lakeperioddata=[(5, "STATUS", "CONSTANT")])
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "Featureno",
        )


def test_lkt_period_auxiliary_ifno_out_of_range(function_tmpdir, targets):
    """apply_period_auxiliary must reject an out-of-range PERIOD IFNO the
    same as apply_period_settings does for non-AUX settings.
    """
    mf6 = targets["mf6"]
    name = "lkt_pauxifno"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    # only 1 lake exists (ncv=1); the PERIOD row below targets feature 5
    # (0-based), i.e. IFNO=6, which is out of range
    _build_gwt_lkt(
        sim,
        name,
        gwfname,
        ncol,
        auxiliary=["myaux"],
        lakeperioddata=[(5, "AUXILIARY", "myaux", 5.0)],
    )
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "IFNO must be greater than 0 and less than or equal to",
        )


CPW = 4183.0
RHOW = 999.728
LHV = 2500.0
CPS = 800.0
RHOS = 2650.0


def _build_gwe_lke(
    sim,
    name,
    gwfname,
    ncol,
    packagedata=None,
    lakeperioddata=None,
):
    gwename = "gwe_" + name
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwename}.ims"), [gwe.name])
    flopy.mf6.ModflowGwedis(
        gwe, nlay=1, nrow=1, ncol=ncol, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=CPW,
        density_water=RHOW,
        latent_heat_vaporization=LHV,
        heat_capacity_solid=CPS,
        density_solid=RHOS,
    )
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.2700)
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGweoc(gwe, budget_filerecord=f"{gwename}.cbc")
    if packagedata is None:
        packagedata = [(0, 35.0, 0.5, 0.1, "mylake")]
    if lakeperioddata is None:
        lakeperioddata = [(0, "STATUS", "CONSTANT"), (0, "TEMPERATURE", 100.0)]
    flopy.mf6.ModflowGwelke(
        gwe,
        boundnames=True,
        packagedata=packagedata,
        lakeperioddata=lakeperioddata,
        flow_package_name="LAK-1",
        pname="LKE-1",
        print_temperature=True,
    )
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{name}.gwfgwe",
    )
    return gwe, gwename


def test_lke_packagedata_ifno_out_of_range(function_tmpdir, targets):
    """Same as test_lkt_packagedata_ifno_out_of_range, but for GWE/LKE."""
    mf6 = targets["mf6"]
    name = "lke_ifno"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    _build_gwe_lke(
        sim, name, gwfname, ncol, packagedata=[(1, 35.0, 0.5, 0.1, "mylake")]
    )
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "LKE PACKAGEDATA IFNO (2) must be greater than 0 and less than or equal",
        )


def test_lke_packagedata_missing_feature(function_tmpdir, targets):
    """Same as test_lkt_packagedata_missing_feature, but for GWE/LKE."""
    mf6 = targets["mf6"]
    name = "lke_missing"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    _build_gwe_lke(sim, name, gwfname, ncol)
    sim.write_simulation()

    lke_fname = function_tmpdir / f"gwe_{name}.lke"
    text = lke_fname.read_text()
    lines = text.splitlines(keepends=True)
    out = []
    skipping = False
    for line in lines:
        if line.strip().upper().startswith("BEGIN PACKAGEDATA"):
            skipping = True
            out.append(line)
            continue
        if line.strip().upper().startswith("END PACKAGEDATA"):
            skipping = False
            out.append(line)
            continue
        if skipping:
            continue
        out.append(line)
    lke_fname.write_text("".join(out))

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            "LKE PACKAGEDATA no data specified for feature 1",
        )


def test_lke_packagedata_out_of_order(function_tmpdir, targets):
    """Same as test_lkt_packagedata_out_of_order, but for GWE/LKE."""
    mf6 = targets["mf6"]
    name = "lke_order"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0e-6, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name, nlakes=2)
    _build_gwe_lke(
        sim,
        name,
        gwfname,
        ncol,
        packagedata=[
            (1, 222.0, 0.5, 0.1, "lake-two"),
            (0, 111.0, 0.5, 0.1, "lake-one"),
        ],
        lakeperioddata=[],
    )
    sim.write_simulation()

    returncode, _ = run_mf6([mf6], str(function_tmpdir))
    assert returncode == 0, "mf6 did not terminate successfully"

    # GWE's much larger eqnsclfac (rhow*cpw) collapses absolute temperature
    # magnitudes for this near-zero timestep, unlike the GWT/LKT case -- but
    # the ratio between features is scale-invariant and preserved exactly,
    # so check that instead of absolute values.
    lst_text = (function_tmpdir / f"gwe_{name}.lst").read_text()
    idx_one = lst_text.index("LAKE-ONE")
    idx_two = lst_text.index("LAKE-TWO")
    temp_one = float(lst_text[idx_one : idx_one + 60].split()[2])
    temp_two = float(lst_text[idx_two : idx_two + 60].split()[2])
    assert temp_two / temp_one == pytest.approx(222.0 / 111.0, rel=1e-3), (
        f"LAKE-TWO/LAKE-ONE temperature ratio {temp_two / temp_one} does not "
        "match the ratio of their own PACKAGEDATA STRT values (2.0) -- STRT "
        "was likely sourced by row order instead of by feature number"
    )


def test_lkt_period_submember_keyword_rejected(function_tmpdir, targets):
    """A compound PERIOD dispatch's own sub-member name must not be usable
    as a top-level dispatch keyword on its own -- it must be reported the
    same as any other unrecognized keyword."""
    mf6 = targets["mf6"]
    name = "lkt_auxval"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf, gwfname, ncol = _build_gwf_lak(sim, name)
    _build_gwt_lkt(sim, name, gwfname, ncol, auxiliary=["myaux"])
    sim.write_simulation()

    # replace the well-formed PERIOD block with one that addresses the
    # AUXILIARY record's own AUXVAL sub-member directly, skipping the
    # AUXILIARY dispatch keyword itself
    lkt_fname = function_tmpdir / f"gwt_{name}.lkt"
    text = lkt_fname.read_text()
    text = re.sub(
        r"BEGIN period  1\n.*?\nEND period  1",
        "BEGIN period  1\n  1  AUXVAL  5.0\nEND period  1",
        text,
        flags=re.DOTALL | re.IGNORECASE,
    )
    lkt_fname.write_text(text)

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(function_tmpdir),
            mf6,
            'Unrecognized keystring keyword "AUXVAL"',
        )
