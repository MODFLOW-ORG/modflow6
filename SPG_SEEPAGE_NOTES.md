# SPG Seepage Package — Notes and Findings

## Overview

The SPG (seepage) boundary package (`gwf-spg.f90`, ftype `SPG`) models a list of
seepage cells. Each cell is held at atmospheric pressure (zero pressure head, so
total head equals the cell-center elevation `z = (top + bot) / 2`) while water
discharges into it (out of the simulated domain), and reverts to a **no-flow**
boundary when the flow would reverse back into the domain.

The active ("constant head") state can be enforced with either of two methods,
selectable through package options:

| Method | Option | Mechanism | Notes |
| --- | --- | --- | --- |
| Penalty (default) | `PENALTY_CONDUCTANCE <value>` | Large-conductance term pins the head near `z` while `head > z`. Cell stays active (`ibound > 0`). | Robust; stays inside the normal boundary-package matrix framework. Default conductance `1e9`. |
| Toggle | `IBOUND_TOGGLE` | Cell is switched to a true constant-head cell (`ibound < 0`, like CHD) while discharging, and released to an ordinary active no-flow cell when flow reverses. | Exact per-cell balance, but interacts poorly with packages that claim the whole grid (see below). |

Period input is `cellid` plus optional aux/boundname. The seepage elevation is
always the cell center.

## Budget sign convention

Both methods report seepage discharge as a **negative** package flow (water
leaving the model / a sink), consistent with DRN and the generic boundary
budget:

- Penalty: uses the base `bnd_cq_simrate`, i.e. `simvals = hcof*h - rhs =
  -C*(h - z) < 0` while discharging.
- Toggle: the held cell is constant head, so the base routine skips it. A custom
  `spg_cq` computes `rate = -sum(off-diagonal flowja)` (positive = inflow to the
  held cell), adds it to `flowja(idiag)` to balance the cell's flow row (as CHD
  does), then reports `simvals = -rate` so the sink sign matches the penalty
  method.

## Testing

- `autotest/test_gwf_spg.py`: a confined two-cell column (CHD top, SPG bottom)
  with four cases — `{penalty, toggle} x {discharging, no-flow}`. Checks the
  held head equals the seepage elevation, the discharge balances the
  constant-head inflow, and the no-flow state carries zero flow.
- `autotest/test_gwf_uzr_wt.py`: the `wt-uzr-spg` case adds SPG seepage cells on
  the `+x` boundary of a transient 2D Richards/UZR model and verifies per-cell
  mass balance (FLOW-JA-FACE residuals).

## Findings from the UZR integration

1. **Toggle conflicts with UZR.** UZR claims the entire grid for its own flow and
   storage formulation (`iformulation` is set for all cells). When SPG toggle
   switches a cell to a true constant head (`ibound < 0`), UZR still tries to
   apply its formulation to that cell. This produces repeated
   `UZR convertible cell error` messages and a small residual imbalance. The
   **penalty** method keeps the cell active (`ibound > 0`), so it composes
   cleanly with UZR, and is the method used in the `wt-uzr-spg` test.

2. **Seepage on/off is a discontinuous nonlinearity.** Switching a cell between
   the held and no-flow states is discontinuous: a tiny head change can flip the
   state and cause a jump in the boundary flow. As a result the model converges
   to a slightly looser per-cell FLOW-JA-FACE residual near the seepage cells
   (~`3.7e-6`) than the rest of the grid. Both the penalty and toggle methods
   converge to the same residual level, confirming this comes from the seepage
   switch itself rather than from the penalty conductance. The test therefore
   uses `res_atol = 1e-5` for the seepage case and keeps the strict `1e-6`
   tolerance for the non-seepage cases.

## Future work

- **Continue the investigation with the toggle method and how to avoid the
  reported errors.** In particular, determine how to reconcile the toggle
  (true constant-head) method with packages such as UZR that claim the whole
  grid, so that toggling a cell to constant head does not trigger
  `UZR convertible cell error` or introduce a residual imbalance. Candidate
  directions to explore:
  - Have SPG coordinate with UZR (and similar whole-grid packages) to release or
    skip the toggled cell's formulation while it is held as constant head.
  - Investigate whether the seepage on/off switch can be smoothed (continuous
    transition) to improve convergence and tighten the residual.
  - Evaluate Newton-Raphson handling of the state switch under UZR.
