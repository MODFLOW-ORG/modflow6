"""
Regional-scale VSF benchmark: full VSF vs VSF_1D_ONLY and IMS settings.

Not collected by pytest. Results are discussed in README.md.

Units: cm and s. Grid: n x n cells of 100 m (--delr, in cm), 10 fine
layers (30 cm) for the unsaturated zone on top of 20 coarse aquifer layers (5 m) with a
heterogeneous K field and an aquitard. Patchy recharge (rain event followed
by a dry period), drains at 1.2 m depth and fixed heads on the lateral
boundary.

Usage: python vsf_big.py WS --n 183 [--delr 1000] [--solver milu0_nobt]
       [--one-d] [--picard] [--build-only]
"""

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

import flopy
import numpy as np
from scipy.ndimage import gaussian_filter

MF6 = Path(__file__).resolve().parents[2] / "bin" / "mf6"

# Celia et al. (1990) Haverkamp sand, parameters in 1/cm
SOIL = dict(
    porosity=0.287, satres=0.26132, alpha=0.027074, n=3.96, beta=0.052408, k=4.74
)
KSOIL = 0.00944  # cm/s
DAY = 86400.0


def kfield(n, seed, mean_k, sigma_ln, corr_cells):
    rng = np.random.default_rng(seed)
    z = gaussian_filter(rng.standard_normal((n, n)), corr_cells, mode="wrap")
    z = (z - z.mean()) / z.std()
    return mean_k * np.exp(sigma_ln * z - 0.5 * sigma_ln**2)


SOLVERS = {
    "ilut5": dict(),  # COMPLEX default: ILUT level 5, droptol 1e-4
    "ilu0": dict(preconditioner_levels=0, preconditioner_drop_tolerance=0.0),
    "milu0": dict(
        preconditioner_levels=0,
        preconditioner_drop_tolerance=0.0,
        relaxation_factor=0.97,
    ),
    "milu0_nobt": dict(
        preconditioner_levels=0,
        preconditioner_drop_tolerance=0.0,
        relaxation_factor=0.97,
        backtracking_number=0,
    ),
    "ilut1": dict(preconditioner_levels=1, preconditioner_drop_tolerance=1e-3),
    "ilu0rcm": dict(
        preconditioner_levels=0,
        preconditioner_drop_tolerance=0.0,
        reordering_method="rcm",
    ),
}


def build(ws, n, one_d, newton, delr=10000.0, solver="ilut5"):
    name = "big"
    nfine, ncoarse = 10, 20
    dz_fine, dz_coarse = 30.0, 500.0
    top = 10000.0
    thick = [dz_fine] * nfine + [dz_coarse] * ncoarse
    botm = top - np.cumsum(thick)
    nlay = nfine + ncoarse
    wt0 = top - 150.0

    # hydraulic conductivity
    k = np.empty((nlay, n, n))
    k[:nfine] = KSOIL
    for ilay in range(ncoarse):
        k[nfine + ilay] = kfield(n, 100 + ilay, 0.01, 1.0, 4.0)
    k[nfine + 6 : nfine + 8] = 1.0e-5  # aquitard
    k33 = 0.1 * k
    k33[:nfine] = KSOIL

    # recharge pattern
    pattern = kfield(n, 7, 1.0, 0.7, 3.0)
    rain = 20.0 / 10.0 / DAY * pattern  # 20 mm/d in cm/s
    dry = 1.0 / 10.0 / DAY * pattern  # 1 mm/d

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=str(MF6))
    sim.name_file.profile_option = "detail"
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="SECONDS",
        nper=2,
        perioddata=[(2 * DAY, 1, 1.0), (3 * DAY, 1, 1.0)],
    )
    ats = [(0, 300.0, 1.0, 6 * 3600.0, 2.0, 2.0), (1, 300.0, 1.0, 6 * 3600.0, 2.0, 2.0)]
    tdis.ats.initialize(maxats=len(ats), perioddata=ats, filename=f"{name}.ats")

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="newton" if newton else None
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=0.01,
        outer_maximum=100,
        inner_maximum=300,
        inner_dvclose=0.001,
        rcloserecord=[10.0, "strict"],
        linear_acceleration="BICGSTAB",
        **SOLVERS[solver],
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=n, ncol=n, delr=delr, delc=delr, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=wt0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k, k33=k33)
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-7, sy=0.0, iconvert=0, transient={0: True})
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        soil_model="haverkamp",
        storage_scheme=None if newton else "chord-slope",
        kr_averaging="arithmetic",
        vsf_1d_only=one_d,
        porosity=SOIL["porosity"],
        satres=SOIL["satres"],
        alphahvk=SOIL["alpha"],
        nhvk=SOIL["n"],
        betahvk=SOIL["beta"],
        khvk=SOIL["k"],
    )
    flopy.mf6.ModflowGwfrcha(gwf, recharge={0: rain, 1: dry})

    # drains at 1.2 m depth (layer 4, center at 1.35 m), ~100 d resistance
    cdrn = delr * delr / (100 * DAY)
    drn = [[(4, i, j), top - 120.0, cdrn] for i in range(n) for j in range(n)]
    flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: drn})

    # fixed heads on the lateral boundary, saturated layers
    zc = top - (np.cumsum(thick) - 0.5 * np.array(thick))
    edge = [
        (i, j) for i in range(n) for j in range(n) if i in (0, n - 1) or j in (0, n - 1)
    ]
    chd = [[(kk, i, j), wt0] for kk in range(nlay) if zc[kk] < wt0 for (i, j) in edge]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd})

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )
    sim.write_simulation(silent=True)
    return sim


def stats(ws):
    lst = (Path(ws) / "mfsim.lst").read_text()
    out = {}
    calls = [int(c) for c in re.findall(r"(\d+) CALLS TO NUMERICAL", lst)]
    inner = [int(c) for c in re.findall(r"(\d+) TOTAL ITERATIONS", lst)]
    out["n_steps"] = len(calls)
    out["n_outer"] = int(sum(calls))
    out["n_inner"] = int(sum(inner))
    m = re.search(r"Elapsed run time:\s*(.*)", lst)
    out["elapsed"] = m.group(1).strip() if m else None
    out["converged"] = "Normal termination" in lst
    for key, label in [
        ("solve", r"Solve \(1\)"),
        ("formulate", "Formulate"),
        ("linear_solve", "Linear solve"),
        ("calc_flows", "Calculate flows"),
    ]:
        m = re.search(rf"\]\s+{label}:\s+([0-9.]+)s", lst)
        out[f"t_{key}"] = float(m.group(1)) if m else None
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("ws")
    p.add_argument("--n", type=int, default=183)
    p.add_argument("--one-d", action="store_true")
    p.add_argument("--picard", action="store_true")
    p.add_argument("--build-only", action="store_true")
    p.add_argument("--delr", type=float, default=10000.0)
    p.add_argument("--solver", default="ilut5", choices=list(SOLVERS))
    a = p.parse_args()
    ws = Path(a.ws)
    ws.mkdir(parents=True, exist_ok=True)
    build(ws, a.n, a.one_d, not a.picard, a.delr, a.solver)
    if a.build_only:
        raise SystemExit
    t0 = time.perf_counter()
    subprocess.run(
        [str(MF6)],
        cwd=ws,
        stdout=open(ws / "stdout.txt", "w"),
        stderr=subprocess.STDOUT,
    )
    s = stats(ws)
    s["wall"] = time.perf_counter() - t0
    s["ncells"] = 30 * a.n * a.n
    (ws / "stats.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))
