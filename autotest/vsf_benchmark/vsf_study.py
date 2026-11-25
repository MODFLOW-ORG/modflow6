"""
Parametrised VSF performance study at regional scale.

Builds a family of regional models (fine unsaturated-zone layers on top of a
heterogeneous aquifer, units cm and s) over a matrix of grid size, soil,
scenario, formulation and solver settings, runs them with a given mf6
executable and stores timings, iteration counts and results under a tag. Two
tags, for example two builds, are compared with ``compare``.

Not collected by pytest. Run the timings serially on an otherwise idle
machine.

Examples::

    python vsf_study.py list --n 18 58 --soil loam sand
    python vsf_study.py run --tag base --exe /path/to/mf6 --n 18 58 \\
        --soil loam --scenario shallow deep --form newton --solver milu0
    python vsf_study.py compare base opt1
"""

import argparse
import itertools
import json
import os
import re
import subprocess
import time
from pathlib import Path

import flopy
import numpy as np
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).resolve().parent
MF6 = HERE.parents[1] / "bin" / "mf6"
OUT = HERE / "study"
DAY = 86400.0
NAME = "vsf"

# Soils, parameters in cm and s. Van Genuchten values after Carsel and
# Parrish (1988), Haverkamp sand after Celia et al. (1990).
SOILS = {
    "sand": dict(model="vg", ks=8.25e-3, ts=0.43, tr=0.045, alpha=0.145, n=2.68),
    "loam": dict(model="vg", ks=2.89e-4, ts=0.43, tr=0.078, alpha=0.036, n=1.56),
    "siltloam": dict(model="vg", ks=1.25e-4, ts=0.45, tr=0.067, alpha=0.020, n=1.41),
    "haverkamp": dict(
        model="hvk",
        ks=9.44e-3,
        ts=0.287,
        tr=0.075,
        alpha=0.027074,
        n=3.96,
        beta=0.052408,
        k=4.74,
    ),
}

# Scenarios: unsaturated-zone discretisation, initial water table depth and
# a two-period recharge history (rate as a multiple of the soil Ks, length).
SCENARIOS = {
    # as shallow, with fixed rates of 20 mm/d then 1 mm/d instead of
    # fractions of Ks (absolute=True: rates in cm/s)
    "rain20": dict(
        nfine=10,
        dzfine=30.0,
        wt_depth=150.0,
        rch=[(2.0 / DAY, 2 * DAY), (0.1 / DAY, 3 * DAY)],
        absolute=True,
        drn_layer=4,
        drn_depth=120.0,
        layered=False,
    ),
    # shallow water table, rain then drainage to drains at 1.2 m
    "shallow": dict(
        nfine=10,
        dzfine=30.0,
        wt_depth=150.0,
        rch=[(0.05, 2 * DAY), (0.0025, 3 * DAY)],
        drn_layer=4,
        drn_depth=120.0,
        layered=False,
    ),
    # deep unsaturated zone, an infiltration front travels down
    "deep": dict(
        nfine=20,
        dzfine=50.0,
        wt_depth=900.0,
        rch=[(0.2, 1 * DAY), (0.001, 4 * DAY)],
        drn_layer=None,
        drn_depth=None,
        layered=False,
    ),
    # rain above Ks on a shallow profile, surface drains take the excess
    "wet": dict(
        nfine=10,
        dzfine=30.0,
        wt_depth=100.0,
        rch=[(1.5, 1 * DAY), (0.01, 2 * DAY)],
        drn_layer=0,
        drn_depth=0.0,
        layered=False,
    ),
    # as shallow, with alternating layers of the soil and ten times less
    # conductive material (same retention, Ks / 10)
    "layered": dict(
        nfine=10,
        dzfine=30.0,
        wt_depth=150.0,
        rch=[(0.05, 2 * DAY), (0.0025, 3 * DAY)],
        drn_layer=4,
        drn_depth=120.0,
        layered=True,
    ),
}

# IMS settings: COMPLEX with BICGSTAB, overridden per solver variant
IMS_BASE = dict(
    print_option="SUMMARY",
    complexity="COMPLEX",
    outer_dvclose=0.01,
    outer_maximum=100,
    inner_maximum=300,
    inner_dvclose=0.001,
    rcloserecord=[10.0, "strict"],
    linear_acceleration="BICGSTAB",
    csv_outer_output_filerecord="ims_outer.csv",
)
MILU0 = dict(
    preconditioner_levels=0,
    preconditioner_drop_tolerance=0.0,
    relaxation_factor=0.97,
    backtracking_number=0,
)
SOLVERS = {
    "complex": dict(),  # ILUT(5, 1e-4), backtracking 20
    "milu0": MILU0,
    # inner (linear) convergence variants on top of MILU0
    "milu0_dv01": dict(MILU0, inner_dvclose=0.01),
    "milu0_rel1": dict(MILU0, rcloserecord=[0.1, "relative_rclose"]),
    "milu0_rel2": dict(MILU0, rcloserecord=[0.01, "relative_rclose"]),
    # nonlinear (under-relaxation, backtracking) variants on top of MILU0;
    # COMPLEX defaults are DBD with theta 0.8, kappa 1e-4, gamma 0
    "milu0_urnone": dict(MILU0, under_relaxation="NONE"),
    "milu0_ursimple": dict(
        MILU0, under_relaxation="SIMPLE", under_relaxation_gamma=0.7
    ),
    "milu0_dbd95": dict(MILU0, under_relaxation_theta=0.95),
    "milu0_dbdg": dict(MILU0, under_relaxation_gamma=0.2),
    "milu0_bt": dict(MILU0, backtracking_number=20),
    # DBD with fast recovery after a flip-flop
    "milu0_dbdk": dict(MILU0, under_relaxation_kappa=0.2),
    "milu0_dbdk_dv01": dict(MILU0, under_relaxation_kappa=0.2, inner_dvclose=0.01),
    "milu0_dbdk5": dict(MILU0, under_relaxation_theta=0.5, under_relaxation_kappa=0.2),
}

NCOARSE, DZCOARSE, TOP = 20, 500.0, 10000.0


def kfield(n, seed, mean_k, sigma_ln, corr_cells):
    rng = np.random.default_rng(seed)
    z = gaussian_filter(rng.standard_normal((n, n)), corr_cells, mode="wrap")
    z = (z - z.mean()) / z.std()
    return mean_k * np.exp(sigma_ln * z - 0.5 * sigma_ln**2)


def disu_grid(nlay, n, delr, botm, order):
    """DISU arrays for a regular grid with layer-first or column-first node
    numbers; returns the DISU keyword arguments and the node number of each
    (layer, row, column) cell."""
    shape = (nlay, n, n)
    if order == "layer":
        nodeid = np.arange(nlay * n * n).reshape(shape)
    else:
        nodeid = np.arange(nlay * n * n).reshape((n, n, nlay)).transpose(2, 0, 1)
    top = np.concatenate([[TOP], botm[:-1]])
    dz = top - botm
    kk, ii, jj = np.meshgrid(np.arange(nlay), np.arange(n), np.arange(n), indexing="ij")
    conn = []  # node, neighbour, ihc, cl12, hwva
    for dk, di, dj in [
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    ]:
        k2, i2, j2 = kk + dk, ii + di, jj + dj
        ok = (k2 >= 0) & (k2 < nlay) & (i2 >= 0) & (i2 < n) & (j2 >= 0) & (j2 < n)
        vert = dk != 0
        conn.append(
            np.stack(
                [
                    nodeid[ok],
                    nodeid[k2[ok], i2[ok], j2[ok]],
                    np.full(ok.sum(), 0 if vert else 1),
                    0.5 * dz[kk[ok]] if vert else np.full(ok.sum(), 0.5 * delr),
                    np.full(ok.sum(), delr * delr if vert else delr),
                ],
                axis=1,
            )
        )
    nodes = nlay * n * n
    diag = np.stack(
        [
            np.arange(nodes),
            np.arange(nodes),
            np.ones(nodes),
            np.zeros(nodes),
            np.zeros(nodes),
        ],
        axis=1,
    )
    c = np.concatenate([diag] + conn)
    # diagonal first in each row, then the neighbours in ascending order
    isdiag = c[:, 0] == c[:, 1]
    c = c[np.lexsort((c[:, 1], ~isdiag, c[:, 0]))]
    iac = np.bincount(c[:, 0].astype(int), minlength=nodes)

    def bynode(a):
        out = np.empty(nodes)
        out[nodeid.ravel()] = np.broadcast_to(a, shape).ravel()
        return out

    kw = dict(
        nodes=nodes,
        nja=len(c),
        iac=iac,
        ja=c[:, 1].astype(int),
        ihc=c[:, 2].astype(int),
        cl12=c[:, 3],
        hwva=c[:, 4],
        top=bynode(top[:, None, None]),
        bot=bynode(botm[:, None, None]),
        area=np.full(nodes, delr * delr),
    )
    return kw, nodeid, bynode


# VSF options set from the command line, with their case id suffixes; None
# leaves the option out (program default)
VSF_OPTS = {"dsat": "newton_dsat_max", "krsmooth": "kr_smoothing", "hae": "haevgn"}


def case_id(c):
    cid = (
        f"n{c['n']}_{c['soil']}_{c['scenario']}_{c['form']}_{c['solver']}"
        f"_dt{c['dt0']:g}" + ("" if c.get("grid", "dis") == "dis" else f"_{c['grid']}")
    )
    for k in VSF_OPTS:
        if c.get(k) is not None:
            cid += f"_{k}{c[k]:g}"
    return cid


def build(
    ws,
    n,
    soil,
    scenario,
    form,
    solver,
    delr=1000.0,
    dt0=300.0,
    grid="dis",
    dsat=None,
    krsmooth=None,
    hae=None,
):
    sc, so = SCENARIOS[scenario], SOILS[soil]
    newton = form == "newton"
    nfine = sc["nfine"]
    thick = np.array([sc["dzfine"]] * nfine + [DZCOARSE] * NCOARSE)
    botm = TOP - np.cumsum(thick)
    zc = TOP - (np.cumsum(thick) - 0.5 * thick)
    nlay = nfine + NCOARSE
    wt0 = TOP - sc["wt_depth"]
    shape = (nlay, n, n)

    # conductivity: soil in the fine layers, lognormal aquifer, aquitard
    k = np.empty(shape)
    k[:nfine] = so["ks"]
    if sc["layered"]:
        k[1:nfine:2] = 0.1 * so["ks"]
    for ilay in range(NCOARSE):
        k[nfine + ilay] = kfield(n, 100 + ilay, 0.01, 1.0, 4.0)
    k[nfine + 6 : nfine + 8] = 1.0e-5
    k33 = 0.1 * k
    k33[:nfine] = k[:nfine]

    def full(v):
        return arr(np.full(shape, v))

    pattern = kfield(n, 7, 1.0, 0.7, 3.0)
    tend = sum(p[1] for p in sc["rch"])

    sim = flopy.mf6.MFSimulation(sim_name=NAME, sim_ws=ws, exe_name=str(MF6))
    sim.name_file.profile_option = "detail"
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="SECONDS",
        nper=len(sc["rch"]),
        perioddata=[(p[1], 1, 1.0) for p in sc["rch"]],
    )
    ats = [(i, dt0, 1.0, 6 * 3600.0, 2.0, 2.0) for i in range(len(sc["rch"]))]
    tdis.ats.initialize(maxats=len(ats), perioddata=ats, filename=f"{NAME}.ats")

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=NAME, newtonoptions="newton" if newton else None
    )
    ims = flopy.mf6.ModflowIms(sim, **dict(IMS_BASE, **SOLVERS[solver]))
    sim.register_ims_package(ims, [gwf.name])
    if grid == "dis":
        flopy.mf6.ModflowGwfdis(
            gwf, nlay=nlay, nrow=n, ncol=n, delr=delr, delc=delr, top=TOP, botm=botm
        )

        def arr(a):
            return a

        def cell(kk, i, j):
            return (kk, i, j)
    else:
        kw, nodeid, arr = disu_grid(nlay, n, delr, botm, grid.split("_")[1])
        flopy.mf6.ModflowGwfdisu(gwf, **kw)
        np.save(Path(ws) / "nodeid.npy", nodeid)

        def cell(kk, i, j):
            return (int(nodeid[kk, i, j]),)

    flopy.mf6.ModflowGwfic(gwf, strt=wt0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=arr(k), k33=arr(k33))
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-7, sy=0.0, iconvert=0, transient={0: True})
    common = dict(
        ivarsat=1,
        storage_scheme=None if newton else "chord-slope",
        kr_averaging="arithmetic",
        porosity=full(so["ts"]),
        satres=full(so["tr"] / so["ts"]),
        sat_filerecord=f"{NAME}.sat",
    )
    if dsat is not None:
        common["newton_dsat_max"] = dsat
    if so["model"] == "vg":
        if krsmooth is not None:
            common["kr_smoothing"] = krsmooth
        if hae is not None:
            common["haevgn"] = full(hae)
    if so["model"] == "vg":
        flopy.mf6.ModflowGwfvsf(
            gwf,
            soil_model="vangenuchten",
            alphavgn=full(so["alpha"]),
            nvgn=full(so["n"]),
            **common,
        )
    else:
        flopy.mf6.ModflowGwfvsf(
            gwf,
            soil_model="haverkamp",
            alphahvk=full(so["alpha"]),
            nhvk=full(so["n"]),
            betahvk=full(so["beta"]),
            khvk=full(so["k"]),
            **common,
        )
    scale = 1.0 if sc.get("absolute") else so["ks"]
    rch = {i: f * scale * pattern for i, (f, _) in enumerate(sc["rch"])}
    if grid == "dis":
        flopy.mf6.ModflowGwfrcha(gwf, recharge=rch)
    else:
        spd = {
            ip: [[cell(0, i, j), r[i, j]] for i in range(n) for j in range(n)]
            for ip, r in rch.items()
        }
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data=spd)
    if sc["drn_layer"] is not None:
        cdrn = delr * delr / (100 * DAY)
        drn = [
            [cell(sc["drn_layer"], i, j), TOP - sc["drn_depth"], cdrn]
            for i in range(n)
            for j in range(n)
        ]
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: drn})
    edge = [
        (i, j) for i in range(n) for j in range(n) if i in (0, n - 1) or j in (0, n - 1)
    ]
    chd = [
        [cell(kk, i, j), wt0] for kk in range(nlay) if zc[kk] < wt0 for (i, j) in edge
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd})
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{NAME}.hds", saverecord=[("HEAD", "LAST")]
    )
    sim.write_simulation(silent=True)
    return nlay * n * n, tend


def parse_listing(ws):
    lst = (Path(ws) / "mfsim.lst").read_text()
    out = dict(converged="Normal termination" in lst)
    # accepted steps and their iterations from the IMS step summaries
    out["steps"] = len(re.findall(r"(\d+) CALLS TO NUMERICAL", lst))
    out["outer_ok"] = sum(int(c) for c in re.findall(r"(\d+) CALLS TO NUMERICAL", lst))
    # all outer and inner iterations, including failed step attempts, from
    # the per-iteration table (print_option SUMMARY)
    rows = re.findall(r"^ Model\s+(\d+)\s+(\d+)\s", lst, flags=re.M)
    out["outer"] = len(rows)
    out["inner"] = sum(int(r[1]) for r in rows)
    out["inner_max_hits"] = sum(int(r[1]) >= 300 for r in rows)
    out["failed_attempts"] = lst.count("Solution 1 did not converge")
    out["pivot_warnings"] = lst.count("ADDED SMALL VALUE TO PIVOT")
    for key, label in [
        ("run", "Run"),
        ("solve", r"Solve \(1\)"),
        ("formulate", "Formulate"),
        ("linear", "Linear solve"),
        ("flows", "Calculate flows"),
    ]:
        m = re.search(rf"\]\s+{label}:\s+([0-9.]+)s", lst)
        out[f"t_{key}"] = float(m.group(1)) if m else None
    return out


def cases(a):
    keys = ["n", "soil", "scenario", "form", "solver", "dt0", "grid", *VSF_OPTS]
    vals = [a.n, a.soil, a.scenario, a.form, a.solver, a.dt0, a.grid]
    vals += [getattr(a, k) for k in VSF_OPTS]
    return [dict(zip(keys, v)) for v in itertools.product(*vals)]


def run(a):
    exe = Path(a.exe).resolve()
    tagdir = OUT / a.tag
    tagdir.mkdir(parents=True, exist_ok=True)
    version = subprocess.run(
        [str(exe), "-v"], capture_output=True, text=True
    ).stdout.strip()
    # extra environment for the executable, e.g. experiment switches
    env = dict(kv.split("=", 1) for kv in a.env)
    for c in cases(a):
        cid = case_id(c)
        ws = tagdir / cid
        if (ws / "result.json").exists() and not a.force:
            print(f"{cid}: exists, skipped")
            continue
        ws.mkdir(parents=True, exist_ok=True)
        ncells, tend = build(
            ws,
            c["n"],
            c["soil"],
            c["scenario"],
            c["form"],
            c["solver"],
            dt0=c["dt0"],
            grid=c["grid"],
            **{k: c[k] for k in VSF_OPTS},
        )
        t0 = time.perf_counter()
        try:
            p = subprocess.run(
                [str(exe)],
                cwd=ws,
                capture_output=True,
                text=True,
                timeout=a.timeout,
                env=dict(os.environ, **env),
            )
            timeout = False
        except subprocess.TimeoutExpired:
            timeout = True
        wall = time.perf_counter() - t0
        r = dict(
            c,
            case=cid,
            tag=a.tag,
            exe=str(exe),
            version=version,
            env=env,
            ncells=ncells,
            wall=wall,
            timeout=timeout,
        )
        if not timeout:
            r.update(parse_listing(ws))
            r["returncode"] = p.returncode
        (ws / "result.json").write_text(json.dumps(r, indent=1))
        with open(tagdir / "results.jsonl", "a") as f:
            f.write(json.dumps(r) + "\n")
        print(
            json.dumps(
                {
                    k: r.get(k)
                    for k in (
                        "case",
                        "wall",
                        "converged",
                        "steps",
                        "outer",
                        "inner",
                        "failed_attempts",
                        "t_formulate",
                        "t_linear",
                    )
                }
            )
        )


def last_record(path, text):
    reader = flopy.utils.HeadFile(path, text=text, precision="double")
    data = reader.get_data()
    nodeid = path.parent / "nodeid.npy"
    if nodeid.exists():
        data = data.ravel()[np.load(nodeid)]
    return data.reshape(-1)


def compare(a):
    rows = []
    base, other = OUT / a.base, OUT / a.other
    for rb in sorted(base.glob("*/result.json")):
        ro = other / (rb.parent.name + a.suffix) / "result.json"
        if not ro.exists():
            continue
        b, o = json.loads(rb.read_text()), json.loads(ro.read_text())
        row = dict(
            case=b["case"],
            ncells=b["ncells"],
            wall_base=b["wall"],
            wall_other=o["wall"],
            speedup=b["wall"] / o["wall"],
            outer=f"{b.get('outer')}/{o.get('outer')}",
            inner=f"{b.get('inner')}/{o.get('inner')}",
        )
        if b.get("converged") and o.get("converged"):
            hb = last_record(rb.parent / f"{NAME}.hds", "HEAD")
            ho = last_record(ro.parent / f"{NAME}.hds", "HEAD")
            sb = last_record(rb.parent / f"{NAME}.sat", "SATURATION")
            so = last_record(ro.parent / f"{NAME}.sat", "SATURATION")
            row["max_dh_cm"] = float(np.abs(hb - ho).max())
            row["max_dsat"] = float(np.abs(sb - so).max())
        rows.append(row)
    for r in rows:
        print(json.dumps(r))
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("run", "list"):
        s = sub.add_parser(name)
        s.add_argument("--n", type=int, nargs="+", default=[18])
        s.add_argument("--soil", nargs="+", default=["loam"], choices=SOILS)
        s.add_argument("--scenario", nargs="+", default=["shallow"], choices=SCENARIOS)
        s.add_argument(
            "--form", nargs="+", default=["newton"], choices=["newton", "picard"]
        )
        s.add_argument("--solver", nargs="+", default=["milu0"], choices=SOLVERS)
        s.add_argument(
            "--dt0",
            type=float,
            nargs="+",
            default=[300.0],
            help="initial ATS time step [s]",
        )
        s.add_argument(
            "--grid",
            nargs="+",
            default=["dis"],
            choices=["dis", "disu_layer", "disu_column"],
            help="DIS, or DISU with layer- or column-first node numbers",
        )
        for k, opt in VSF_OPTS.items():
            s.add_argument(
                f"--{k}",
                type=float,
                nargs="+",
                default=[None],
                help=f"VSF {opt.upper()} (default: left out)",
            )
        if name == "run":
            s.add_argument("--tag", required=True)
            s.add_argument("--exe", default=str(MF6))
            s.add_argument("--timeout", type=float, default=3600.0)
            s.add_argument("--force", action="store_true")
            s.add_argument(
                "--env",
                nargs="*",
                default=[],
                metavar="KEY=VALUE",
                help="environment variables for the executable",
            )
    s = sub.add_parser("compare")
    s.add_argument("base")
    s.add_argument("other")
    s.add_argument(
        "--suffix", default="", help="case suffix in other, e.g. _disu_column"
    )
    a = p.parse_args()
    if a.cmd == "list":
        for c in cases(a):
            print(case_id(c))
    elif a.cmd == "run":
        run(a)
    else:
        compare(a)
