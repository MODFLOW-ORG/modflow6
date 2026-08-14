"""
Figure showing what the MST stranded mass option does for a cell that drains
completely and is deactivated for transport.

The water table falls below the bottom of the tested cell, stays there, and
returns. The cell takes no part in the solution while it is inactive, and the
mass its reservoirs hold is retained through that period, decaying only if the
solute decays, and is returned when the cell is rewet.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np

MF6_EXE = Path(__file__).resolve().parent.parent.parent.parent / "bin" / "mf6"

NLAY, NROW, NCOL = 3, 1, 2
DELR = DELC = 100.0
TOP = 30.0
BOTM = [20.0, 10.0, 0.0]
POROSITY, SY, RHOB, KD, CINIT = 0.3, 0.15, 1600.0, 1.0e-4, 100.0
LAMBDA = 0.01
PERLEN = 10.0

#: the water table falls below the bottom of the top layer and then returns
HEADS = [30.0, 22.0, 14.0, 5.0, 5.0, 14.0, 22.0, 30.0, 30.0]
DRY = -1.0e29


def build_model(ws, decay):
    sim = flopy.mf6.MFSimulation(
        sim_name="d", version="mf6", exe_name=str(MF6_EXE), sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(HEADS),
        perioddata=[(PERLEN, 1, 1.0)] * len(HEADS),
    )
    solver = {
        "complexity": "complex",
        "outer_dvclose": 1e-8,
        "inner_dvclose": 1e-9,
        "linear_acceleration": "bicgstab",
        "outer_maximum": 200,
        "inner_maximum": 300,
    }

    gwf = flopy.mf6.ModflowGwf(sim, modelname="f", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, filename="f.ims", **solver), ["f"]
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwfic(gwf, strt=TOP)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_saturation=True,
        icelltype=1,
        k=50.0,
        rewet_record=[("WETFCT", 1.0, "IWETIT", 1, "IHDWET", 0)],
        wetdry=1.0,
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=0.0, sy=SY, transient={0: True}, save_flows=True
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={
            k: [[(NLAY - 1, 0, 1), h, 0.0]] for k, h in enumerate(HEADS)
        },
        auxiliary="CONCENTRATION",
        pname="CHD-1",
        save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="f.hds",
        budget_filerecord="f.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="t", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, filename="t.ims", **solver), ["t"]
    )
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwtic(gwt, strt=CINIT)
    flopy.mf6.ModflowGwtadv(gwt)
    kwargs = {
        "porosity": POROSITY,
        "save_flows": True,
        "sorption": "linear",
        "bulk_density": RHOB,
        "distcoef": KD,
        "stranded_mass": True,
        "stranded_filerecord": [("t.strand.bin",)],
    }
    if decay:
        kwargs.update(first_order_decay=True, decay=LAMBDA, decay_sorbed=LAMBDA)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="t.ucn",
        budget_filerecord="t.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="f", exgmnameb="t")
    sim.write_simulation(silent=True)
    return sim


def run_case(tmp, tag, decay):
    ws = Path(tmp) / tag
    shutil.rmtree(ws, ignore_errors=True)
    build_model(ws, decay)
    r = subprocess.run([str(MF6_EXE)], cwd=ws, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:])
        raise RuntimeError(f"mf6 failed for {tag}")
    cf = flopy.utils.HeadFile(ws / "t.ucn", text="CONCENTRATION")
    times = np.array(cf.get_times())
    conc = np.array([cf.get_data(totim=t)[0, 0, 0] for t in times])
    head = np.array(
        [flopy.utils.HeadFile(ws / "f.hds").get_data(totim=t)[0, 0, 0] for t in times]
    )
    sf = flopy.utils.HeadFile(ws / "t.strand.bin", text="STRANDED")
    st = np.array([sf.get_data(totim=t)[0, 0, 0] for t in times])
    return times, head, conc, st


with tempfile.TemporaryDirectory() as tmp:
    times, head, conc, st = run_case(tmp, "nodecay", False)
    _, _, conc_d, st_d = run_case(tmp, "decay", True)

dry = head < DRY
print(f"the cell is dry for {dry.sum()} of {len(dry)} output times")
print(
    f"stranded held   {st[dry][0] / 1e3:9.3f} -> {st[dry][-1] / 1e3:9.3f} kg (no decay)"
)
print(
    f"stranded decays {st_d[dry][0] / 1e3:9.3f} -> {st_d[dry][-1] / 1e3:9.3f} kg (decay)"
)

# the cell has no head or concentration while it is inactive
hplot = np.where(dry, np.nan, head)
cplot = np.where(dry, np.nan, conc)
cplot_d = np.where(dry, np.nan, conc_d)
t0, t1 = times[dry][0] - PERLEN, times[dry][-1]

figpth = Path(__file__).resolve().parent.parent / "Figures"

with styles.USGSPlot():
    fig, axes = plt.subplots(
        nrows=3, figsize=(6.8, 5.0), sharex=True, constrained_layout=True
    )
    for ax in axes:
        ax.axvspan(t0, t1, color="0.88", zorder=0)
        ax.tick_params(direction="in", top=True, right=True)

    ax = axes[0]
    ax.plot(times, hplot, color="black", lw=1.2)
    ax.axhline(BOTM[0], color="0.5", lw=0.8, ls=":")
    ax.set_ylabel("Head in the tested\ncell, in meters")
    ax.set_ylim(0.0, 34.0)
    styles.heading(ax=ax, letter="A")

    ax = axes[1]
    ax.plot(times, cplot, color="#d62728", lw=1.2, label="no decay")
    ax.plot(times, cplot_d, color="#d62728", lw=1.2, ls="--", label="first-order decay")
    ax.set_ylabel("Concentration, in grams\nper cubic meter")
    styles.heading(ax=ax, letter="B")
    styles.graph_legend(
        ax=ax,
        ncols=2,
        loc="upper right",
        fontsize=7,
        handlelength=2.0,
        framealpha=0.9,
        edgecolor="0.7",
    )

    ax = axes[2]
    ax.plot(times, st / 1000.0, color="#d62728", lw=1.2)
    ax.plot(times, st_d / 1000.0, color="#d62728", lw=1.2, ls="--")
    ax.set_ylabel("Stranded mass,\nin kilograms")
    ax.set_xlabel("Time, in days")
    ax.set_xlim(times[0], times[-1])
    styles.heading(ax=ax, letter="C")

    fig.savefig(figpth / "MSTStrandedMassDry.pdf", dpi=300)
print(f"Saved {figpth / 'MSTStrandedMassDry.pdf'}")
