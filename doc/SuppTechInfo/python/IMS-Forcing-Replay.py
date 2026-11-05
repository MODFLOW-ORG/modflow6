"""Generate the linear-solve residual history figure for the Inexact Newton
Forcing Terms chapter of the MODFLOW 6 Supplemental Technical Information
document.

The residual histories are from single linear solves of a fixed Newton system
for a continental-scale steady-state model, solved with five different
groupings of the model subdomains onto processors. They are archived in
IMS-Forcing-Replay-Residuals.csv, so the figure does not need the model runs.

The figure is written to ../Figures/IMSForcingReplay.pdf, and the inner
iterations needed to reach each forcing term are printed for the table.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.plot import styles

pth = os.path.dirname(os.path.abspath(__file__))
datafile = os.path.join(pth, "IMS-Forcing-Replay-Residuals.csv")
figpth = os.path.join(pth, "..", "Figures")

# processor groupings, the labels used in the text, and line styles
GROUPINGS = {
    "016p": ("16 processors", "#1b9e77", "-"),
    "032p": ("32 processors", "#d95f02", "-"),
    "032p_regroup": ("32 processors, regrouped", "#7570b3", "--"),
    "032p_placebo": ("32 processors, placebo", "#e7298a", ":"),
    "064p": ("64 processors", "#66a61e", "-"),
}
# outer iteration the starting heads were taken from, and the panel letter
STATES = {0: "A", 100: "B"}
# forcing terms marked on the figure and tabulated
ETA_LINES = (1e-1, 1e-3, 1e-6)
ETA_TABLE = (3e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-6)


def iterations_to_reach(rel, eta):
    """First inner iteration with a relative residual at or below eta."""
    idx = np.flatnonzero(rel <= eta)
    return int(idx[0]) if idx.size else None


df = pd.read_csv(datafile)

print("state grouping " + " ".join(f"{e:g}" for e in ETA_TABLE))
with styles.USGSPlot():
    plt.rcParams["ps.fonttype"] = 42
    fig, axd = plt.subplot_mosaic(
        [list(STATES.values())], figsize=(7.5, 3.9), layout="constrained"
    )
    for state, letter in STATES.items():
        ax = axd[letter]
        sel = df[df["state"] == state]
        imax = 0
        for key, (label, color, ls) in GROUPINGS.items():
            d = sel[sel["grouping"] == key].sort_values("iteration")
            k = d["iteration"].to_numpy()
            r = d["true_residual"].to_numpy()
            rel = r / r[0]
            imax = max(imax, k[-1])
            ax.semilogy(k, rel, color=color, ls=ls, lw=1.0, label=label)
            its = [iterations_to_reach(rel, e) for e in ETA_TABLE]
            print(f"{state:5d} {key:12s} " + " ".join(f"{i}" for i in its))
        for eta in ETA_LINES:
            ax.axhline(eta, color="0.5", lw=0.5, ls="-.", zorder=0)
            ax.text(
                0.98,
                eta * 1.6,
                f"$\\eta = 10^{{{int(np.log10(eta))}}}$",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=7,
            )
        ax.set_xlim(0, imax)
        ax.set_ylim(1e-13, 1e1)
        styles.heading(
            ax=ax, letter=letter, heading=f"Heads after {state} outer iterations"
        )
        styles.xlabel(ax=ax, label="Inner iteration")
        if letter == "A":
            styles.ylabel(ax=ax, label="Relative linear residual, dimensionless")
    # one explanation for both panels, below them
    handles, labels = axd["A"].get_legend_handles_labels()
    leg = fig.legend(
        handles, labels, loc="outside lower center", ncol=3, title="EXPLANATION"
    )
    leg.get_title().set_fontweight("bold")
    for ax in axd.values():
        styles.remove_edge_ticks(ax)
    fig.savefig(os.path.join(figpth, "IMSForcingReplay.pdf"))
