"""
Figure showing what controls the size of the effect of the MST stranded mass
option, for the oscillating water table of MST-StrandedMass.py.

The retained water volume is the change in the mobile water volume less the
water released to storage, so the option is controlled by the specific yield
relative to the porosity, and the sorbed part by the distribution coefficient.
Each case is run with the option off and on, and the two quantities plotted are
the mass held out of the mobile domain and the departure in concentration from
the run without the option, both as percentages so that they share an axis.
"""

import subprocess
import tempfile
from pathlib import Path

import flopy
import flopy.plot.styles as styles
import matplotlib.pyplot as plt
import numpy as np

MF6_EXE = Path(__file__).resolve().parent.parent.parent.parent / "bin" / "mf6"

TOP, BOTM = 40.0, 0.0
DELR = DELC = 100.0
VCELL = DELR * DELC * (TOP - BOTM)
POROSITY = 0.3
RHOB = 1600.0
CINIT = 100.0

PERLEN = 1.0
NPER = 400
TCYCLE = 100.0
HMEAN, HAMP = 24.0, 14.0

#: specific yields, from the value at which no water is retained downward
SY_VALUES = [0.30, 0.25, 0.20, 0.15, 0.10, 0.05]
#: distribution coefficients, from no sorption upward
KD_VALUES = [0.0, 2.5e-5, 5.0e-5, 1.0e-4, 2.0e-4, 4.0e-4]
SY_BASE = 0.15
KD_BASE = 1.0e-4
#: the specific yield sweep omits sorption, so that it isolates the water
#: retained against drainage and shows it vanish where the two are equal
KD_SY_SWEEP = 0.0


def head_schedule():
    """Sinusoidal head in the constant head cell, one value per stress period."""
    t = np.arange(1, NPER + 1) * PERLEN
    h = HMEAN + HAMP * np.sin(2.0 * np.pi * t / TCYCLE)
    return np.concatenate(([h[0]], h))


def build_model(ws, sy, kd, stranded):
    hds = head_schedule()
    sim = flopy.mf6.MFSimulation(
        sim_name="s", version="mf6", exe_name=str(MF6_EXE), sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=len(hds), perioddata=[(PERLEN, 1, 1.0)] * len(hds)
    )
    solver = {
        "complexity": "simple",
        "outer_dvclose": 1e-9,
        "inner_dvclose": 1e-10,
        "linear_acceleration": "bicgstab",
    }

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="f", save_flows=True, newtonoptions="NEWTON"
    )
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, filename="f.ims", **solver), ["f"]
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwfic(gwf, strt=hds[0])
    flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=True, save_saturation=True, icelltype=1, k=100.0
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=sy, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={k: [[(0, 0, 1), h, 0.0]] for k, h in enumerate(hds)},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
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
        gwt, nlay=1, nrow=1, ncol=2, delr=DELR, delc=DELC, top=TOP, botm=BOTM
    )
    flopy.mf6.ModflowGwtic(gwt, strt=CINIT)
    kwargs = {"porosity": POROSITY, "save_flows": True}
    if kd > 0.0:
        kwargs.update(sorption="linear", bulk_density=RHOB, distcoef=kd)
    if stranded:
        kwargs.update(stranded_mass=True, stranded_filerecord=[("t.strand.bin",)])
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


def run_case(tmp, tag, sy, kd, stranded):
    ws = Path(tmp) / tag
    build_model(ws, sy, kd, stranded)
    r = subprocess.run([str(MF6_EXE)], cwd=ws, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:])
        raise RuntimeError(f"mf6 failed for {tag}")
    conc = flopy.utils.HeadFile(ws / "t.ucn", text="CONCENTRATION")
    times = conc.get_times()
    c = np.array([conc.get_data(totim=t).flatten()[0] for t in times])
    s = None
    if stranded:
        sf = flopy.utils.HeadFile(ws / "t.strand.bin", text="STRANDED")
        s = np.array([sf.get_data(totim=t).flatten()[0] for t in times])
    return c, s


def sweep(tmp, values, vary):
    """Peak stranded mass and peak concentration departure, both in percent."""
    mass, conc = [], []
    for i, v in enumerate(values):
        sy = v if vary == "sy" else SY_BASE
        kd = v if vary == "kd" else KD_SY_SWEEP
        c_off, _ = run_case(tmp, f"{vary}{i}off", sy, kd, False)
        c_on, s_on = run_case(tmp, f"{vary}{i}on", sy, kd, True)
        # mass of the tested cell when the simulation starts
        sat0 = (head_schedule()[0] - BOTM) / (TOP - BOTM)
        m0 = VCELL * sat0 * (POROSITY + RHOB * kd) * CINIT
        mass.append(100.0 * s_on.max() / m0)
        conc.append(100.0 * np.abs(c_on - c_off).max() / CINIT)
        print(
            f"  {vary}={v:<9.3g} stranded {mass[-1]:6.2f} %"
            f"  departure {conc[-1]:6.2f} %"
        )
    return np.array(mass), np.array(conc)


with tempfile.TemporaryDirectory() as tmp:
    print("specific yield:")
    sy_mass, sy_conc = sweep(tmp, SY_VALUES, "sy")
    print("distribution coefficient:")
    kd_mass, kd_conc = sweep(tmp, KD_VALUES, "kd")

figpth = Path(__file__).resolve().parent.parent / "Figures"

with styles.USGSPlot():
    fig, axes = plt.subplots(ncols=2, figsize=(6.8, 3.0), constrained_layout=True)

    ax = axes[0]
    ratio = np.array(SY_VALUES) / POROSITY
    ax.plot(
        ratio,
        sy_mass,
        color="#d62728",
        marker="o",
        ms=3.5,
        lw=1.2,
        label="mass held out of the mobile domain",
    )
    ax.plot(
        ratio,
        sy_conc,
        color="#1f77b4",
        marker="s",
        ms=3.5,
        lw=1.2,
        label="departure in concentration",
    )
    ax.set_xlabel("Specific yield divided by porosity\n(no sorption)")
    ax.set_ylabel("Percentage of the initial mass\nor of the initial concentration")
    ax.set_xlim(1.05, 0.10)
    ax.set_ylim(0.0, max(sy_mass.max(), sy_conc.max()) * 1.45)
    ax.tick_params(direction="in", top=True, right=True)
    styles.heading(ax=ax, letter="A")
    styles.graph_legend(
        ax=ax,
        ncols=1,
        loc="upper left",
        fontsize=6.5,
        handlelength=2.0,
        framealpha=0.9,
        edgecolor="0.7",
    )

    ax = axes[1]
    kd = np.array(KD_VALUES) * 1.0e4
    ax.plot(kd, kd_mass, color="#d62728", marker="o", ms=3.5, lw=1.2)
    ax.plot(kd, kd_conc, color="#1f77b4", marker="s", ms=3.5, lw=1.2)
    ax.set_xlabel(
        "Distribution coefficient, in units of\n"
        r"$1 \times 10^{-4}$ cubic meters per kilogram"
        "\n(specific yield half the porosity)"
    )
    ax.set_ylim(0.0, max(kd_mass.max(), kd_conc.max()) * 1.45)
    ax.tick_params(direction="in", top=True, right=True)
    styles.heading(ax=ax, letter="B")

    fig.savefig(figpth / "MSTStrandedMassSensitivity.pdf", dpi=300)
print(f"Saved {figpth / 'MSTStrandedMassSensitivity.pdf'}")
