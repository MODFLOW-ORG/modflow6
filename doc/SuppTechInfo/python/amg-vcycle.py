"""Generate the AMG multigrid cycle-schedule figure.

Produces a single-panel figure showing the 4-level V-cycle (gamma=1) grid
schedule, analogous to Figure 3 in the MODFLOW GMG solver documentation
(Wilson and Naff, 2004).

Run from the repository root or from doc/SuppTechInfo/python/.
Output: doc/SuppTechInfo/Figures/amg-vcycle.pdf
"""

from pathlib import Path

import flopy.plot.styles as fstyles
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

figpth = Path(__file__).resolve().parent.parent / "Figures"
figpth.mkdir(parents=True, exist_ok=True)
width = 6.8
dpi = 300

# 4-level hierarchy; finest = level 0, coarsest = level 3 (= L-1)
nlevels = 4

# V-cycle (gamma=1): visit every level once on the way down and up
vcycle = list(range(nlevels)) + list(range(nlevels - 2, -1, -1))
# [0, 1, 2, 3, 2, 1, 0]

ytick_labels = [
    r"$\ell = 0$  (fine)",
    r"$\ell = 1$",
    r"$\ell = 2$",
    r"$\ell = L-1$  (coarse)",
]

with fstyles.USGSPlot():
    fig, ax = plt.subplots(1, 1, figsize=(width * 0.55, width * 0.42))

    x = np.arange(len(vcycle))
    ax.plot(x, vcycle, "ko-", markersize=5, linewidth=1.2, zorder=5)
    ax.set_title("V-cycle (one recursive call)")
    ax.set_xticks([])
    ax.set_xlabel("Cycle schedule")
    ax.set_xlim(-0.5, len(vcycle) - 0.5)
    for lev in range(nlevels):
        ax.axhline(lev, color="0.80", linewidth=0.5, zorder=0)

    # finest level (ell=0) at top, coarsest (ell=L-1) at bottom
    ax.set_ylim(nlevels - 0.5, -0.5)
    ax.set_yticks(range(nlevels))
    ax.set_yticklabels(ytick_labels)
    ax.set_ylabel("Multigrid level")

    outfile = figpth / "amg-vcycle.pdf"
    fig.savefig(outfile, bbox_inches="tight", dpi=dpi)
    print(f"Saved {outfile}")
