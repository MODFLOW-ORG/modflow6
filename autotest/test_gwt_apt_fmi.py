"""
FMI (flow-model-interface) coverage for the GWT/GWE APT packages.

Every other test exercising any of the 8 APT packages (SFT/MWT/LKT/UZT,
SFE/MWE/LKE/UZE) couples to a live, in-sim GWF model, where find_*_package
resolves flowbudptr from that live package's own BudgetObjectType and ncv
traces back to a declared DIMENSIONS value the flow package itself read.

FMI coupling is architecturally different: the GWF model runs separately
(here, as its own simulation) and the transport model reads its saved
budget file after the fact.  find_*_package's flows_from_file branch
resolves flowbudptr via fmi%set_aptbudobj_pointer, which derives ncv from
the binary budget file's own metadata -- there is no live flow package and
no declared DIMENSIONS value backing it.

These tests check that PACKAGEDATA validation against ncv holds up just as
well when ncv's provenance is a budget file instead of a live package's
DIMENSIONS block.  Deferred ncv also means the PERIOD block's own
row-count sizing can't rely on a DIMENSIONS value either;
test_lkt_fmi_maxbound checks that case specifically.
"""

import flopy
import pytest
from test_gwt_apt_errors import run_mf6, run_mf6_error


def _build_and_run_gwf_lak(ws, exe, name, nlakes=1, nnodes=5):
    """Standalone GWF+LAK model; returns the flow-file paths FMI needs.

    nlakes=1 keeps the original 3-connection single-lake layout; nlakes>1
    gives each lake one connection, legitimately sharing the nnodes cells.
    """
    gwfname = "gwf_" + name
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 5, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY", complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam")
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=nnodes, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=0,
        k=20.0,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), -0.5], [(0, 0, nnodes - 1), -0.6]],
        pname="CHD-1",
        filename=f"{gwfname}.chd",
    )
    connlen = connwidth = 0.5
    if nlakes == 1:
        packagedata = [(0, -0.4, 3)]
        connectiondata = [
            (0, 0, (0, 0, 1), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", 0.0, 10, 10, connlen, connwidth),
        ]
        perioddata = [(0, "STATUS", "CONSTANT"), (0, "STAGE", -0.4)]
    else:
        packagedata = [(n, -0.4, 1) for n in range(nlakes)]
        connectiondata = [
            (n, 0, (0, 0, n % nnodes), "VERTICAL", 0.0, 10, 10, connlen, connwidth)
            for n in range(nlakes)
        ]
        perioddata = [(n, "STATUS", "CONSTANT") for n in range(nlakes)]
    flopy.mf6.ModflowGwflak(
        gwf,
        save_flows=True,
        nlakes=nlakes,
        noutlets=0,
        ntables=0,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="LAK-1",
        budget_filerecord=f"{gwfname}.lak.bud",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation()
    returncode, _ = run_mf6([exe], ws)
    assert returncode == 0, "standalone GWF+LAK run did not terminate successfully"
    return {
        "GWFHEAD": f"{gwfname}.hds",
        "GWFBUDGET": f"{gwfname}.bud",
        "LAK-1": f"{gwfname}.lak.bud",
    }


def _build_gwt_lkt_fmi(
    ws, exe, name, flow_files, packagedata=None, lakeperioddata=None, nnodes=5
):
    """nnodes must match the flow model's own grid -- FMI has no way to
    detect a mismatch unless a GWFGRID entry is supplied, and a mismatch
    is invalid input (see mf6io's FMI section), not something to test here.
    """
    gwtname = "gwt_" + name
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 5, 1.0)])
    flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", complexity="SIMPLE", linear_acceleration="BICGSTAB"
    )
    gwt = flopy.mf6.MFModel(
        sim, model_type="gwt6", modelname=gwtname, model_nam_file=f"{gwtname}.nam"
    )
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=nnodes, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.30)
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtoc(gwt, budget_filerecord=f"{gwtname}.cbc")
    pd = [
        ("GWFHEAD", flow_files["GWFHEAD"], None),
        ("GWFBUDGET", flow_files["GWFBUDGET"], None),
        ("LAK-1", flow_files["LAK-1"], None),
    ]
    flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)
    if packagedata is None:
        packagedata = [(0, 35.0, "mylake")]
    if lakeperioddata is None:
        lakeperioddata = [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 100.0)]
    flopy.mf6.ModflowGwtlkt(
        gwt,
        boundnames=True,
        packagedata=packagedata,
        lakeperioddata=lakeperioddata,
        flow_package_name="LAK-1",
        pname="LKT-1",
        print_concentration=True,
    )
    sim.write_simulation()
    return sim, gwtname


def test_lkt_fmi_success(function_tmpdir, targets):
    """LKT must resolve ncv and PACKAGEDATA correctly when flow is FMI-only."""
    mf6 = targets["mf6"]
    flow_ws = function_tmpdir / "flow"
    flow_ws.mkdir()
    flow_files = _build_and_run_gwf_lak(str(flow_ws), mf6, "fmiflow")
    flow_files = {k: f"../flow/{v}" for k, v in flow_files.items()}

    transport_ws = function_tmpdir / "transport"
    transport_ws.mkdir()
    _build_gwt_lkt_fmi(str(transport_ws), mf6, "fmilkt", flow_files)

    returncode, _ = run_mf6([mf6], str(transport_ws))
    assert returncode == 0, "FMI-coupled LKT simulation did not terminate successfully"

    lst_text = (transport_ws / "gwt_fmilkt.lst").read_text()
    assert "NUMBER OF CONTROL VOLUMES = 1" in lst_text, (
        "LKT did not resolve ncv=1 from the FMI-supplied LAK-1 budget file"
    )


def test_lkt_fmi_packagedata_ifno_out_of_range(function_tmpdir, targets):
    """PACKAGEDATA IFNO validation must still reject an out-of-range value
    when ncv comes from an FMI budget file instead of a live flow
    package's declared dimension."""
    mf6 = targets["mf6"]
    flow_ws = function_tmpdir / "flow"
    flow_ws.mkdir()
    flow_files = _build_and_run_gwf_lak(str(flow_ws), mf6, "fmiflow")
    flow_files = {k: f"../flow/{v}" for k, v in flow_files.items()}

    transport_ws = function_tmpdir / "transport"
    transport_ws.mkdir()
    # the FMI budget file has 1 lake (ncv=1); IFNO=2 (0-based: 1) is out of range
    _build_gwt_lkt_fmi(
        str(transport_ws),
        mf6,
        "fmibad",
        flow_files,
        packagedata=[(1, 35.0, "mylake")],
    )

    with pytest.raises(RuntimeError):
        run_mf6_error(
            str(transport_ws),
            mf6,
            "LKT PACKAGEDATA IFNO (2) must be greater than 0 and less than or equal",
        )


def test_lkt_fmi_maxbound(function_tmpdir, targets):
    """PERIOD-block sizing must scale with PACKAGEDATA's own feature count,
    not grid node count, when ncv is FMI-deferred: many lakes legitimately
    sharing a few cells must not overflow a grid-sized allocation."""
    mf6 = targets["mf6"]
    nlakes, nnodes = 20, 3
    flow_ws = function_tmpdir / "flow"
    flow_ws.mkdir()
    flow_files = _build_and_run_gwf_lak(
        str(flow_ws), mf6, "fmiflow", nlakes=nlakes, nnodes=nnodes
    )
    flow_files = {k: f"../flow/{v}" for k, v in flow_files.items()}

    transport_ws = function_tmpdir / "transport"
    transport_ws.mkdir()
    packagedata = [(n, 35.0, f"lake{n}") for n in range(nlakes)]
    lakeperioddata = []
    for n in range(nlakes):
        lakeperioddata.append((n, "STATUS", "CONSTANT"))
        lakeperioddata.append((n, "CONCENTRATION", 100.0 + n))
    _build_gwt_lkt_fmi(
        str(transport_ws),
        mf6,
        "fmibig",
        flow_files,
        packagedata=packagedata,
        lakeperioddata=lakeperioddata,
        nnodes=nnodes,
    )

    returncode, _ = run_mf6([mf6], str(transport_ws))
    assert returncode == 0, "FMI-coupled LKT simulation did not terminate successfully"

    lst_text = (transport_ws / "gwt_fmibig.lst").read_text()
    assert "exceeds pre-allocated maxbound" not in lst_text
