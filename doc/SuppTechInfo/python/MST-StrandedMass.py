"""
Figures showing the effect of the MST stranded mass option on a cell with an
oscillating water table.

A single active cell is drained and rewetted by a sinusoidal head in the cell
next to it. Five transport simulations are run on the same flow field:

  base            no sorption, no decay
  sorption        sorption, stranded mass off
  sorption+SM     sorption, stranded mass on
  sorp+decay      sorption and first-order decay, stranded mass off
  sorp+decay+SM   sorption and first-order decay, stranded mass on

Three figures are written: the oscillating water table, the comparison without
decay, and the comparison with decay.
"""

import subprocess
import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np

MF6_EXE = Path(__file__).resolve().parent.parent.parent.parent / "bin" / "mf6"

# ---------------------------------------------------------------------------
# Problem parameters
# ---------------------------------------------------------------------------
TOP, BOTM = 40.0, 0.0
DELR = DELC = 100.0
AREA = DELR * DELC
VCELL = AREA * (TOP - BOTM)
POROSITY = 0.3
SY = 0.15
RHOB = 1600.0
KD = 1.0e-4
CINIT = 100.0
DECAY = 5.0e-3

PERLEN = 1.0
NPER = 400
TCYCLE = 100.0
HMEAN, HAMP = 24.0, 14.0


def head_schedule():
    """Sinusoidal head in the constant head cell, one value per stress period.

    The first period holds the water table still, which lets the saturation of
    the previous time step initialize before any mass is stranded.
    """
    t = np.arange(1, NPER + 1) * PERLEN
    h = HMEAN + HAMP * np.sin(2.0 * np.pi * t / TCYCLE)
    return np.concatenate(([h[0]], h))


def build_model(ws, sorption, decay, stranded):
    name = "sm"
    gwfname, gwtname = "gwf-" + name, "gwt-" + name
    hds = head_schedule()

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=str(MF6_EXE), sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=len(hds), perioddata=[(PERLEN, 1, 1.0)] * len(hds)
    )

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwfname}.ims",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, newtonoptions="NEWTON"
    )
    sim.register_ims_package(imsgwf, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwfic(gwf, strt=hds[0])
    flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=True, save_saturation=True, icelltype=1, k=100.0
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=SY, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={k: [[(0, 0, 1), h, 0.0]] for k, h in enumerate(hds)},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        budget_filerecord=f"{gwfname}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=2, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwtic(gwt, strt=CINIT)

    kwargs = {"porosity": POROSITY, "save_flows": True}
    if sorption:
        kwargs.update(sorption="linear", bulk_density=RHOB, distcoef=KD)
    if decay:
        kwargs.update(first_order_decay=True, decay=DECAY, decay_sorbed=DECAY)
    if stranded:
        kwargs.update(
            stranded_mass=True,
            stranded_filerecord=[(f"{gwtname}.strand.bin",)],
        )
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)

    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        budget_filerecord=f"{gwtname}.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    sim.write_simulation(silent=True)
    return sim


def run_case(tmp, tag, sorption, decay, stranded):
    ws = Path(tmp) / tag
    build_model(ws, sorption, decay, stranded)
    r = subprocess.run([str(MF6_EXE)], cwd=ws, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise RuntimeError(f"mf6 failed for case {tag}")

    conc = flopy.utils.HeadFile(ws / "gwt-sm.ucn", text="CONCENTRATION")
    times = np.array(conc.get_times())
    c = np.array([conc.get_data(totim=t).flatten()[0] for t in times])
    head = flopy.utils.HeadFile(ws / "gwf-sm.hds")
    h = np.array([head.get_data(totim=t).flatten()[0] for t in times])
    s = None
    if stranded:
        sf = flopy.utils.HeadFile(ws / "gwt-sm.strand.bin", text="STRANDED")
        s = np.array([sf.get_data(totim=t).flatten()[0] for t in times])
    return times, h, c, s, case_mass(ws, times, sorption, stranded)


def case_mass(ws, times, sorption, stranded):
    """Solute mass in the model, dissolved, sorbed, and stranded."""
    conc = flopy.utils.HeadFile(ws / "gwt-sm.ucn", text="CONCENTRATION")
    head = flopy.utils.HeadFile(ws / "gwf-sm.hds")
    sf = None
    if stranded:
        sf = flopy.utils.HeadFile(ws / "gwt-sm.strand.bin", text="STRANDED")
    fac = POROSITY + (RHOB * KD if sorption else 0.0)
    m = np.zeros(times.shape, dtype=float)
    for i, t in enumerate(times):
        c = conc.get_data(totim=t).flatten()
        h = head.get_data(totim=t).flatten()
        sat = np.clip((h - BOTM) / (TOP - BOTM), 0.0, 1.0)
        m[i] = float((VCELL * sat * fac * c).sum())
        if sf is not None:
            m[i] = m[i] + float(sf.get_data(totim=t).sum())
    return m


CASES = [
    ("base", False, False, False, "current approach", "#1f77b4", "-"),
    ("basesm", False, False, True, "stranded mass", "#d62728", "-"),
    ("sorb", True, False, False, "sorption", "#1f77b4", "-"),
    ("sorbsm", True, False, True, "sorption, stranded mass", "#d62728", "-"),
    ("dcy", True, True, False, "sorption and decay", "#1f77b4", "--"),
    ("dcysm", True, True, True, "sorption and decay, stranded mass", "#d62728", "--"),
]

results = {}
with tempfile.TemporaryDirectory() as tmp:
    for tag, sorb, dcy, sm, *_ in CASES:
        results[tag] = run_case(tmp, tag, sorb, dcy, sm)
        print(f"ran {tag}")

times = results["base"][0]
heads = results["base"][1]

figpth = Path(__file__).resolve().parent.parent / "Figures"


#: the case each stranded mass case is compared against
PAIRED = {"basesm": "base", "sorbsm": "sorb", "dcysm": "dcy"}


def comparison_figure(fname, tags, caption_tags):
    """Concentration, solute mass, and stranded mass for a set of cases."""
    with styles.USGSPlot():
        fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(6.8, 5.6),
            sharex=True,
            constrained_layout=True,
        )

        ax = axes[0]
        for tag in tags:
            label, color, ls = LOOKUP[tag]
            ax.plot(times, results[tag][2], color=color, ls=ls, lw=1.2, label=label)
        ax.set_ylabel("Concentration, in grams\nper cubic meter")
        ax.set_ylim(0.0, 265.0)
        ax.tick_params(direction="in", top=True, right=True)
        styles.heading(ax=ax, letter="A")
        styles.graph_legend(
            ax=ax,
            ncols=1,
            loc="upper right",
            fontsize=7,
            handlelength=2.0,
            framealpha=0.9,
            edgecolor="0.7",
        )

        ax = axes[1]
        for tag in tags:
            label, color, ls = LOOKUP[tag]
            ax.plot(
                times, results[tag][4] / 1000.0, color=color, ls=ls, lw=1.2, label=label
            )
        ax.set_ylabel("Solute mass in the\nmodel, in kilograms")
        # -- headroom above the curves for the explanation
        mmax = max(results[tag][4].max() for tag in tags) / 1000.0
        ax.set_ylim(0.0, 1.42 * mmax)
        ax.tick_params(direction="in", top=True, right=True)
        styles.heading(ax=ax, letter="B")
        styles.graph_legend(
            ax=ax,
            ncols=1,
            loc="upper right",
            fontsize=7,
            handlelength=2.0,
            framealpha=0.9,
            edgecolor="0.7",
        )

        ax = axes[2]
        for tag in tags:
            if results[tag][3] is None:
                continue
            label, color, ls = LOOKUP[tag]
            ax.plot(
                times,
                results[tag][3] / 1000.0,
                color=color,
                ls=ls,
                lw=1.2,
                label=label,
            )
        ax.set_ylabel("Stranded mass, in kilograms")
        ax.set_ylim(0.0, 9800.0)
        ax.set_xlabel("Time, in days")
        ax.set_xlim(0.0, times[-1])
        ax.tick_params(direction="in", top=True, right=True)
        styles.heading(ax=ax, letter="C")
        styles.graph_legend(
            ax=ax,
            ncols=1,
            loc="upper right",
            fontsize=7,
            handlelength=2.0,
            framealpha=0.9,
            edgecolor="0.7",
        )

    fig.savefig(figpth / fname, dpi=300)
    print(f"Saved {figpth / fname}")


LOOKUP = {tag: (label, color, ls) for tag, _, _, _, label, color, ls in CASES}

# -- figure 1, the oscillating water table
with styles.USGSPlot():
    fig, ax = plt.subplots(figsize=(6.8, 2.6), constrained_layout=True)
    ax.plot(times, heads, color="black", lw=1.2)
    ax.set_ylabel("Head, in meters")
    ax.set_xlabel("Time, in days")
    ax.set_ylim(0.0, TOP)
    ax.set_xlim(0.0, times[-1])
    ax.tick_params(direction="in", top=True, right=True)
    ax2 = ax.twinx()
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel("Saturation")
fig.savefig(figpth / "MSTStrandedMassHead.pdf", dpi=300)
print(f"Saved {figpth / 'MSTStrandedMassHead.pdf'}")

# -- figure 2, without decay; figure 3, with decay
comparison_figure("MSTStrandedMassNoSorption.pdf", ["base", "basesm"], None)
comparison_figure(
    "MSTStrandedMassSorption.pdf", ["sorb", "sorbsm", "dcy", "dcysm"], None
)

# -- summary numbers quoted in the chapter text
for tag, _, _, _, label, _, _ in CASES:
    c = results[tag][2]
    s = results[tag][3]
    smax = "n/a" if s is None else f"{s.max() / 1000.0:9.1f} kg"
    print(f"{label:38s} c_max={c.max():8.3f}  c_end={c[-1]:8.3f}  strand_max={smax}")

for tag, _, _, _, label, _, _ in CASES:
    m = results[tag][4] / 1000.0
    print(
        f"{label:38s} mass min={m.min():9.1f} kg  max={m.max():9.1f} kg  "
        f"start={m[0]:9.1f} kg  end={m[-1]:9.1f} kg"
    )

for tag, ref in PAIRED.items():
    d = 100.0 * (results[tag][4] - results[ref][4]) / results[ref][4][0]
    print(
        f"{LOOKUP[tag][0]:38s} mass difference min={d.min():6.2f} %  "
        f"max={d.max():6.2f} %  end={d[-1]:6.2f} %"
    )
