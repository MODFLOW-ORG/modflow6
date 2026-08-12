"""
Figure showing the effect of the MST stranded mass option on a cell with an
oscillating water table.

A single active cell is drained and rewetted by a sinusoidal head in the cell
next to it. Five transport simulations are run on the same flow field:

  base            no sorption, no decay
  sorption        sorption, stranded mass off
  sorption+SM     sorption, stranded mass on
  sorp+decay      sorption and first-order decay, stranded mass off
  sorp+decay+SM   sorption and first-order decay, stranded mass on

Panel A shows the head, panel B the simulated concentration, and panel C the
mass held out of the mobile domain.
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
SY = 0.3
RHOB = 1600.0
KD = 1.0e-4
CINIT = 100.0
DECAY = 5.0e-3

PERLEN = 5.0
NPER = 80
TCYCLE = 100.0
HMEAN, HAMP = 24.0, 14.0


def head_schedule():
    """Sinusoidal head in the constant head cell, one value per stress period."""
    t = np.arange(1, NPER + 1) * PERLEN
    return HMEAN + HAMP * np.sin(2.0 * np.pi * t / TCYCLE)


def build_model(ws, sorption, decay, stranded):
    name = "sm"
    gwfname, gwtname = "gwf-" + name, "gwt-" + name
    hds = head_schedule()

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=str(MF6_EXE), sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=NPER, perioddata=[(PERLEN, 1, 1.0)] * NPER
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
    flopy.mf6.ModflowGwfic(gwf, strt=TOP)
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
    return times, h, c, s


CASES = [
    ("base", False, False, False, "no sorption", "0.45", "-"),
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

with styles.USGSPlot():
    fig, axes = plt.subplots(
        nrows=3, ncols=1, figsize=(6.8, 7.0), sharex=True, constrained_layout=True
    )

    # -- A, head and saturation
    ax = axes[0]
    ax.plot(times, heads, color="black", lw=1.2)
    ax.set_ylabel("Head, in meters")
    ax.set_ylim(0.0, TOP)
    ax.tick_params(direction="in", top=True, right=True)
    ax2 = ax.twinx()
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel("Saturation")
    styles.heading(ax=ax, letter="A")

    # -- B, concentration
    ax = axes[1]
    for tag, _, _, _, label, color, ls in CASES:
        ax.plot(times, results[tag][2], color=color, ls=ls, lw=1.2, label=label)
    ax.set_ylabel("Concentration, in grams per cubic meter")
    ax.set_ylim(0.0, 480.0)
    ax.tick_params(direction="in", top=True, right=True)
    styles.heading(ax=ax, letter="B")
    styles.graph_legend(
        ax=ax,
        ncols=2,
        loc="upper center",
        fontsize=7,
        handlelength=2.0,
        framealpha=0.9,
        edgecolor="0.7",
    )

    # -- C, mass held out of the mobile domain
    ax = axes[2]
    for tag, _, _, _, label, color, ls in CASES:
        if results[tag][3] is None:
            continue
        ax.plot(
            times,
            results[tag][3] / 1000.0,
            color=color,
            ls=ls,
            lw=1.2,
            label=label,
        )
    ax.set_ylabel("Stranded mass, in kilograms")
    ax.set_ylim(0.0, 3200.0)
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

figpth = Path(__file__).resolve().parent.parent / "Figures"
fig.savefig(figpth / "MSTStrandedMass.pdf", dpi=300)
print(f"Saved {figpth / 'MSTStrandedMass.pdf'}")

# -- summary numbers quoted in the chapter text
for tag, _, _, _, label, _, _ in CASES:
    c = results[tag][2]
    s = results[tag][3]
    smax = "n/a" if s is None else f"{s.max() / 1000.0:9.1f} kg"
    print(f"{label:38s} c_max={c.max():8.3f}  c_end={c[-1]:8.3f}  strand_max={smax}")
