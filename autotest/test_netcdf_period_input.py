"""
NetCDF period input exported in validate mode (EXPORT_ARRAY_NETCDF) and
read back. Every exported period slice must hold the input given for that
period (state):
    - an array given as values holds them until given again;
    - grid array (compressed) input holds values only at the period's
      boundary cells; an empty block removes all boundaries;
    - a time-array series and the periods it is linked holds the fill
      value; a layer array (READASARRAYS) series link ends at the next
      PERIOD block, after which the array holds the fill value until
      given values;
    - an array never given holds the fill value, except the layer index
      arrays (IRCH, IEVT), which hold 1 until given.

The exported file is then read back with two package files, each run
compared with the ASCII run:
    netcdf: the original period blocks, array lines converted to NETCDF
        (TIMEARRAYSERIES lines kept);
    state_full: a PERIOD block for every period giving every array with a
        value as NETCDF, and every linked series.

Each case runs on a DIS grid with structured and layered mesh NetCDF
files, and on a DISV grid with layered mesh files.

Cases (boundary cells of the grid array cases do not change between
periods while an array carries forward):
    rcha_multi: RCHA, IRCH and RECHARGE in period 1, no block in 2,
        RECHARGE in 3, RECHARGE series in 4, no block in 5, IRCH only in 6.
    rcha_aux: RCHA with auxiliary variables TAG1 (values), TAG2 (series,
        ends at the period 3 block) and TAG3 (never given).
    evta_multi: EVTA, SURFACE, RATE and DEPTH in period 1, RATE only in 2,
        no block in 3.
    welg_carry: WELG, Q and AUX in period 1, no block in period 2, Q at the
        same cells in period 3 (AUX carries forward).
    welg_aux2: WELG with two auxiliary variables, TAG2 never given.
    rivg_empty: RIVG, all arrays in period 1, empty block in period 2, all
        arrays in period 3.
    chdg_change: CHDG, HEAD in period 1, HEAD at the same cells in period 2.
"""

import re

import flopy
import netCDF4
import numpy as np
import pytest
from framework import TestFramework

DNODATA = 3.0e30
nlay, nrow, ncol = 2, 3, 3
ncpl = nrow * ncol
cell_a, cell_b, cell_chd = (0, 0, 0), (1, 2, 2), (0, 2, 2)
NOBLOCK = object()

# case: (package file extension, budget text, mode)
cases = {
    "rcha_multi": ("rcha", "RCH", "layer"),
    "rcha_aux": ("rcha", "RCH", "layer"),
    "evta_multi": ("evta", "EVT", "layer"),
    "welg_carry": ("welg", "WEL", "compressed"),
    "welg_aux2": ("welg", "WEL", "compressed"),
    "rivg_empty": ("rivg", "RIV", "compressed"),
    "chdg_change": ("chdg", "CHD", "compressed"),
}
# (discretization, NetCDF file format)
grids = [("dis", "structured"), ("dis", "mesh"), ("disv", "mesh")]


def grid_values(values):
    """Full grid array (layer, row, column), DNODATA except at the given
    cells."""
    a = np.full((nlay, nrow, ncol), DNODATA)
    for cell, v in values.items():
        a[cell] = v
    return a


def layer_values(v):
    return np.full((nrow, ncol), v, dtype=type(v))


def canonical(a):
    """Grid array as (layer, cell in layer), layer array as (cell in
    layer); cells in layer are ordered by row, then column (also the DISV
    cell order)."""
    a = np.asarray(a)
    return a.reshape(nlay, ncpl) if a.ndim == 3 else a.reshape(ncpl)


def cellid(cell, dis):
    k, i, j = cell
    return cell if dis == "dis" else (k, i * ncol + j)


def model_array(a, dis):
    """Array in the shape the flopy package expects for the grid."""
    return a if dis == "dis" else canonical(a)


def array_text(a, integer=False):
    """INTERNAL array control record and values, one row per line."""
    fmt = "{:d}" if integer else "{:.8e}"
    rows = np.asarray(a).reshape(-1, ncol)
    return "    INTERNAL\n" + "".join(
        "      " + " ".join(fmt.format(v) for v in row) + "\n" for row in rows
    )


def package_text(options, spec, nper, integer=()):
    """Package file from a period spec: {kper: {tag: array | ("series", name)}}."""
    text = "BEGIN options\n" + "".join(f"  {o}\n" for o in options)
    text += "END options\n\n"
    for kper in range(1, nper + 1):
        blk = spec.get(kper, NOBLOCK)
        if blk is NOBLOCK:
            continue
        text += f"BEGIN period {kper}\n"
        for tag, v in blk.items():
            if isinstance(v, tuple):
                text += f"  {tag} TIMEARRAYSERIES {v[1]}\n"
            else:
                text += f"  {tag}\n" + array_text(v, tag in integer)
        text += f"END period {kper}\n\n"
    return text


def series_text(name, records):
    """Time-array series file with constant records {time: value}."""
    text = f"BEGIN attributes\n  NAME {name}\n  METHOD STEPWISE\nEND attributes\n\n"
    for t, v in records.items():
        text += f"BEGIN time {t}\n  CONSTANT {v}\nEND time\n\n"
    return text


def remove_array(text, tag):
    """Remove an array (tag line and its data) from the period blocks."""
    out, skip = [], False
    for line in text.splitlines(keepends=True):
        w = line.split()
        if w and len(line) - len(line.lstrip()) <= 2:
            skip = w[0].lower() == tag
        if not skip:
            out.append(line)
    return "".join(out)


def add_discretization(gwf, dis):
    if dis == "dis":
        flopy.mf6.ModflowGwfdis(
            gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=10.0, botm=[5.0, 0.0]
        )
        return
    # unit cells ordered by row (from the top), then column
    vertices = [
        [i * (ncol + 1) + j, float(j), float(nrow - i)]
        for i in range(nrow + 1)
        for j in range(ncol + 1)
    ]
    cell2d = []
    for i in range(nrow):
        for j in range(ncol):
            iv = i * (ncol + 1) + j
            verts = [iv, iv + 1, iv + ncol + 2, iv + ncol + 1]
            cell2d.append([i * ncol + j, j + 0.5, nrow - i - 0.5, 4, *verts])
    flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=nlay,
        ncpl=ncpl,
        nvert=len(vertices),
        top=10.0,
        botm=[5.0, 0.0],
        vertices=vertices,
        cell2d=cell2d,
    )


def base_model(ws, name, nper, dis):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-9, inner_dvclose=1e-9)
    add_discretization(gwf, dis)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    return sim, gwf


def add_oc(gwf, name):
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )


def build_layer(ws, name, dis):
    """RCHA or EVTA with hand-written package (and series) files."""
    nper = {"rcha_multi": 6, "rcha_aux": 4}.get(name, 3)
    sim, gwf = base_model(ws, name, nper, dis)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[cellid(cell_chd, dis), 10.0]]}
    )
    series = {}
    if name == "rcha_multi":
        flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, filename=f"{name}.rcha")
        spec = {
            1: {"irch": layer_values(2), "recharge": layer_values(0.001)},
            3: {"recharge": layer_values(0.002)},
            4: {"recharge": ("series", "rts")},
            6: {"irch": layer_values(1)},
        }
        options = ["READASARRAYS", f"TAS6 FILEIN {name}.rcha.tas"]
        series = {"rts": {3.0: 0.003, 4.0: 0.004, 6.0: 0.004}}
        integer = ("irch",)
    elif name == "rcha_aux":
        flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, filename=f"{name}.rcha")
        spec = {
            1: {
                "recharge": layer_values(0.001),
                "tag1": layer_values(5.0),
                "tag2": ("series", "ats"),
            },
            3: {"recharge": layer_values(0.002)},
        }
        options = [
            "READASARRAYS",
            "AUXILIARY TAG1 TAG2 TAG3",
            f"TAS6 FILEIN {name}.rcha.tas",
        ]
        series = {"ats": {0.0: 1.0, 1.0: 2.0, 4.0: 2.0}}
        integer = ("irch",)
    else:
        flopy.mf6.ModflowGwfevta(
            gwf, surface=10.0, rate=0.0, depth=1.0, filename=f"{name}.evta"
        )
        spec = {
            1: {
                "surface": layer_values(9.0),
                "rate": layer_values(0.001),
                "depth": layer_values(2.0),
            },
            2: {"rate": layer_values(0.002)},
        }
        options = ["READASARRAYS"]
        integer = ("ievt",)
    add_oc(gwf, name)
    sim.write_simulation()
    ext = cases[name][0]
    (ws / f"{name}.{ext}").write_text(package_text(options, spec, nper, integer))
    for sname, records in series.items():
        (ws / f"{name}.{ext}.tas").write_text(series_text(sname, records))
    return sim, spec, nper


def build_grid(ws, name, dis):
    nper = 3 if name in ("welg_carry", "rivg_empty") else 2
    sim, gwf = base_model(ws, name, nper, dis)
    chd = {0: [[cellid(cell_chd, dis), 10.0]]}
    ab = {cell_a: 1.0, cell_b: 1.0}
    if name == "welg_carry":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
        q1 = grid_values({cell_a: -1.0, cell_b: -2.0})
        q3 = grid_values({cell_a: -3.0, cell_b: -4.0})
        mult = grid_values({cell_a: 0.5, cell_b: 2.0})
        flopy.mf6.ModflowGwfwelg(
            gwf,
            maxbound=2,
            auxiliary="WEL_MULT",
            auxmultname="WEL_MULT",
            q={0: model_array(q1, dis), 2: model_array(q3, dis)},
            aux={0: [model_array(mult, dis)]},
        )
        spec = {1: {"q": q1, "wel_mult": mult}, 3: {"q": q3}}
    elif name == "welg_aux2":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
        q = grid_values({cell_a: -1.0, cell_b: -2.0})
        tag1 = grid_values({cell_a: 5.0, cell_b: 6.0})
        # flopy writes every auxiliary array; TAG2 is removed below
        flopy.mf6.ModflowGwfwelg(
            gwf,
            maxbound=2,
            auxiliary=["TAG1", "TAG2"],
            q={0: model_array(q, dis)},
            aux={0: [model_array(tag1, dis), model_array(tag1, dis)]},
        )
        spec = {1: {"q": q, "tag1": tag1}}
    elif name == "rivg_empty":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
        stage = grid_values({cell_a: 9.0, cell_b: 8.0})
        cond = grid_values(ab)
        rbot = grid_values({cell_a: 7.0, cell_b: 6.0})
        flopy.mf6.ModflowGwfrivg(
            gwf,
            maxbound=2,
            stage={0: model_array(stage, dis), 1: model_array(stage + 0.5, dis)},
            cond={0: model_array(cond, dis), 1: model_array(cond, dis)},
            rbot={0: model_array(rbot, dis), 1: model_array(rbot, dis)},
        )
        spec = {
            1: {"stage": stage, "cond": cond, "rbot": rbot},
            2: {},
            3: {"stage": stage + 0.5, "cond": cond, "rbot": rbot},
        }
    else:
        h1 = grid_values({cell_a: 10.0, cell_b: 9.0})
        h2 = grid_values({cell_a: 10.0, cell_b: 8.0})
        flopy.mf6.ModflowGwfchdg(
            gwf, maxbound=2, head={0: model_array(h1, dis), 1: model_array(h2, dis)}
        )
        spec = {1: {"head": h1}, 2: {"head": h2}}
    add_oc(gwf, name)
    sim.write_simulation()
    if name == "rivg_empty":
        # flopy wrote periods 1 and 2; period 2 becomes 3 after an empty 2
        fpth = ws / f"{name}.rivg"
        text = fpth.read_text().replace("BEGIN period  2", "BEGIN period  3")
        text = text.replace("END period  2", "END period  3")
        tag = "BEGIN period  3"
        fpth.write_text(
            text.replace(tag, f"BEGIN period  2\nEND period  2\n\n{tag}", 1)
        )
    elif name == "welg_aux2":
        fpth = ws / f"{name}.welg"
        fpth.write_text(remove_array(fpth.read_text(), "tag2"))
    return sim, spec, nper


def build_models(name, dis, test):
    if cases[name][2] == "layer":
        sim, spec, nper = build_layer(test.workspace, name, dis)
    else:
        sim, spec, nper = build_grid(test.workspace, name, dis)
    # expected input in the canonical shape
    test.spec = {
        k: {t: v if isinstance(v, tuple) else canonical(v) for t, v in blk.items()}
        for k, blk in spec.items()
    }
    test.nper = nper
    return sim, None


def expected_state(spec, nper, tags, mode, aux=(), index_tag=None):
    """Expected input given per period: {kper: {tag: array or None}},
    and the series linked per period: {kper: {tag: name}}."""
    state = dict.fromkeys(tags)
    if index_tag is not None:
        state[index_tag] = np.ones(ncpl, dtype=int)
    linked = {}
    out, links = {}, {}
    for kper in range(1, nper + 1):
        blk = spec.get(kper, NOBLOCK)
        if blk is not NOBLOCK:
            if mode == "layer":
                # a series link ends at the next block
                for t in list(linked):
                    if t not in blk:
                        state[t] = None
                        linked.pop(t)
            if mode == "compressed" and not blk:
                state = dict.fromkeys(tags)
                linked = {}
            for t, v in blk.items():
                if isinstance(v, tuple):
                    state[t] = None
                    linked[t] = v[1]
                else:
                    state[t] = v
                    linked.pop(t, None)
            mains = [
                t for t, v in blk.items() if t not in aux and not isinstance(v, tuple)
            ]
            if mode == "compressed" and mains:
                bnd = np.logical_and.reduce([blk[t] != DNODATA for t in mains])
                for t in tags:
                    if state[t] is not None:
                        state[t] = np.where(bnd, state[t], DNODATA)
        out[kper] = dict(state)
        links[kper] = dict(linked)
    return out, links


def period_data(ds):
    """{tag: (data, fill, aux)} for the period arrays, data as (time, layer,
    cell in layer) or (time, cell in layer); layered mesh variables are
    joined by layer and aux arrays named by auxiliary name."""
    data, fill, aux = {}, {}, {}
    for v in ds.variables.values():
        if "modflow_input" not in v.ncattrs() or v.dimensions[0] != "time":
            continue
        pkg = v.modflow_input.split("/")[-2]
        if "modflow_iaux" in v.ncattrs():
            tag = re.sub(r"_l\d+$", "", v.name[len(pkg) + 1 :]).lower()
        else:
            tag = v.modflow_input.split("/")[-1].lower()
        a = np.asarray(v[:])
        nt = a.shape[0]
        if "layer" in v.ncattrs():
            data.setdefault(tag, {})[v.layer] = a.reshape(nt, ncpl)
        elif a.ndim == 4:
            data[tag] = a.reshape(nt, nlay, ncpl)
        else:
            data[tag] = a.reshape(nt, ncpl)
        fill[tag] = v._FillValue
        aux[tag] = "modflow_iaux" in v.ncattrs()
    for tag, a in data.items():
        if isinstance(a, dict):
            data[tag] = np.stack([a[k] for k in sorted(a)], axis=1)
    return {t: (data[t], fill[t], aux[t]) for t in data}


def check_state(ncfile, expected):
    """Every exported period slice holds the expected input."""
    with netCDF4.Dataset(ncfile) as ds:
        ds.set_auto_mask(False)
        pdata = period_data(ds)
    for kper, state in expected.items():
        for tag, (a, fill, _) in pdata.items():
            slc, exp = a[kper - 1], state.get(tag)
            if exp is None:
                assert np.all(slc == fill), (
                    f"{tag} period {kper}: expected fill, got {np.unique(slc)}"
                )
            else:
                assert np.allclose(slc, exp), (
                    f"{tag} period {kper}: expected {np.unique(exp)}, "
                    f"got {np.unique(slc)}"
                )


def to_netcdf_lines(text, block):
    """Convert array lines of the named blocks to NETCDF, dropping their
    data; TIMEARRAYSERIES lines are kept."""
    out, inblk, skip = [], False, False
    for line in text.splitlines(keepends=True):
        w = line.split()
        if len(w) > 1 and w[0].upper() == "BEGIN" and w[1].lower() == block:
            inblk, skip = True, False
        elif len(w) > 1 and w[0].upper() == "END" and w[1].lower() == block:
            inblk, skip = False, False
        elif inblk and w:
            indent = len(line) - len(line.lstrip())
            if indent <= 2:
                skip = False
                if len(w) > 1 and w[1].upper() == "TIMEARRAYSERIES":
                    out.append(line)
                    continue
                out.append(f"  {w[0]}  NETCDF\n")
                skip = True
                continue
            if skip:
                continue
        out.append(line)
    return "".join(out)


def full_blocks(text, expected, links, nper):
    """Package file with a PERIOD block for every period giving every array
    with a value as NETCDF and every linked series."""
    out = [text[: text.upper().index("BEGIN PERIOD")]]
    for kper in range(1, nper + 1):
        lines = [f"  {t} TIMEARRAYSERIES {n}" for t, n in links[kper].items()]
        lines += [
            f"  {t}  NETCDF"
            for t, a in expected[kper].items()
            if a is not None and t not in links[kper]
        ]
        out.append(f"BEGIN period {kper}\n" + "".join(f"{ln}\n" for ln in lines))
        out.append(f"END period {kper}\n\n")
    return "".join(out)


def compare(ws, ref, budtext, label):
    hds = next(ref.glob("*.hds")).name
    h = flopy.utils.HeadFile(ws / hds).get_alldata()
    href = flopy.utils.HeadFile(ref / hds).get_alldata()
    assert np.allclose(h, href, atol=1e-10), (
        f"{label}: heads differ from the ASCII run: max diff {np.abs(h - href).max()}"
    )
    cbc = next(ref.glob("*.cbc")).name
    q = flopy.utils.CellBudgetFile(ws / cbc, precision="double").get_data(text=budtext)
    qref = flopy.utils.CellBudgetFile(ref / cbc, precision="double").get_data(
        text=budtext
    )
    assert len(q) == len(qref) > 0, f"{label}: {budtext} budget records differ"
    for a, b in zip(q, qref):
        assert np.allclose(a["q"], b["q"], atol=1e-10), (
            f"{label}: {budtext} rates differ"
        )


def copy_files(src, dst, skip):
    dst.mkdir()
    for f in src.iterdir():
        if f.is_file() and f.suffix not in skip:
            (dst / f.name).write_bytes(f.read_bytes())


def add_option(text, option):
    return text.replace("END options", f"  {option}\nEND options", 1)


def check_output(name, dis, fmt, test):
    ext, budtext, mode = cases[name]
    exe = test.targets["mf6"]
    ref = test.workspace
    pkgfile = f"{name}.{ext}"
    ncfile = f"{name}.nc"
    nam = (ref / f"{name}.nam").read_text()
    grid_packages = (dis, "ic", "npf")
    fileout = "NETCDF_STRUCTURED" if fmt == "structured" else "NETCDF_MESH2D"

    # export in validate mode
    exp = ref / "export"
    copy_files(ref, exp, (".hds", ".cbc", ".lst", ".grb"))
    (exp / f"{name}.nam").write_text(add_option(nam, f"{fileout}  FILEOUT  {ncfile}"))
    for f in [f"{name}.{p}" for p in grid_packages] + [pkgfile]:
        (exp / f).write_text(add_option((ref / f).read_text(), "EXPORT_ARRAY_NETCDF"))
    success, _ = flopy.run_model(
        exe, exp / "mfsim.nam", model_ws=exp, cargs=["--mode=validate"], silent=True
    )
    assert success, "validate mode NetCDF export failed"

    # exported slices hold the input given
    with netCDF4.Dataset(exp / ncfile) as ds:
        pdata = period_data(ds)
    tags = list(pdata)
    aux = [t for t, (_, _, isaux) in pdata.items() if isaux]
    index_tag = {"rcha": "irch", "evta": "ievt"}.get(ext)
    expected, links = expected_state(
        test.spec, test.nper, tags, mode, aux=aux, index_tag=index_tag
    )
    check_state(exp / ncfile, expected)

    # NetCDF input; TIMEARRAYSERIES lines stay in the package file
    pkgtext = to_netcdf_lines((ref / pkgfile).read_text(), "period")
    variants = {
        "netcdf": pkgtext,
        "state_full": full_blocks(pkgtext, expected, links, test.nper),
    }
    for label, text in variants.items():
        ws = ref / label
        copy_files(exp, ws, (".lst",))
        (ws / f"{name}.nam").write_text(add_option(nam, f"NETCDF  FILEIN  {ncfile}"))
        for p in grid_packages:
            f = f"{name}.{p}"
            (ws / f).write_text(to_netcdf_lines((ref / f).read_text(), "griddata"))
        (ws / pkgfile).write_text(text)
        success, _ = flopy.run_model(exe, ws / "mfsim.nam", model_ws=ws, silent=True)
        assert success, f"{label}: NetCDF input run failed"
        compare(ws, ref, budtext, label)


@pytest.mark.netcdf
@pytest.mark.parametrize("dis, fmt", grids)
@pytest.mark.parametrize("name", list(cases))
def test_mf6model(name, dis, fmt, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(name, dis, t),
        check=lambda t: check_output(name, dis, fmt, t),
        targets=targets,
        overwrite=False,
    )
    test.run()
