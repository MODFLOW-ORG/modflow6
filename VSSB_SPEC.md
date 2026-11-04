# VSSB: Variably Saturated Surface Boundary package — v1 specification

Status: draft for review. Supersedes the earlier concept note (`VSSB_concept.md`).
Package type: `VSSB6` (GWF model). Requires the `VSF6` package in the same model.

## 1. Purpose

VSSB is the land-surface boundary for VSF (Richards) models. It applies precipitation to
the top face of a set of surface cells, limits infiltration by the soil, stores excess
water in a local pond, lets groundwater exfiltrate into that pond, and exports water that
cannot be stored as runoff through the MODFLOW 6 water mover (MVR).

Recharge, ponding, runoff, and exfiltration are results of the solution; none of them is
prescribed.

## 2. Scope

In scope for v1:

- precipitation on a fixed set of top-face cells,
- infiltration limited by the soil (`kr`-dependent conductance), driven by total head,
- exfiltration of groundwater through the same face (replaces the surface role of VSSP),
- per-cell pond storage with a static maximum depth,
- runoff as an MVR provider,
- Newton and Picard formulations,
- budget, observations, and optional pond-depth output.

Out of scope for v1 (documented gaps or later work):

- evaporation from ponded water (neither VSSB nor VSET handles it; a documented gap),
- negative precipitation,
- lateral movement of pond water; this belongs to the future surface-water (overland-flow)
  model that receives runoff,
- MVR as a receiver,
- GWT/GWE coupling (SSM registration),
- array-based input (RCH-style); a later addition that must define the full surface-cell set
  once, not per period,
- restart of pond state (pond starts at 0 on every run),
- changes to VSSP. VSSP is kept unchanged for lateral and bottom faces. Adopting the shared
  seepage kernel, Newton terms, and MVR provider support in VSSP is a follow-up.

## 3. Relationship to other packages

| Package | Relationship |
|---|---|
| VSF | Required. VSSB reads `KREL`/`DKRDH`, K and geometry from NPF, which VSF populates. VSF currently rejects mixed domains (`IVARSAT = 0` is fatal), so VSSB assumes every cell is a Richards cell. |
| VSET | Independent. Root uptake stays in VSET. Evaporation from ponded water is not covered by either package. |
| VSSP | Remains for lateral/bottom seepage faces. A node may not be listed in both VSSB and VSSP when the VSSP entry has `IHC = 0` (see §7). |
| RCH, EVT, UZF | Not replaced. VSSB is the Richards-domain surface boundary. |
| MVR | VSSB is a provider only. Runoff is the provided water. |
| Future overland-flow model | Receives runoff through MVR (or an equivalent exchange). VSSB does not route water laterally. |

## 4. Formulation

### 4.1 Notation (per surface cell *i*, node *n*)

| Symbol | Meaning |
|---|---|
| `A` | horizontal cell area (`dis%get_area(n)`) |
| `z_s` | surface elevation = cell top, `top(n)` |
| `h` | cell head, current iterate |
| `Δz` | cell thickness, `top(n) − bot(n)` |
| `K` | vertical conductivity: K33 when specified, otherwise K11 (as in VSSP) |
| `kr`, `dkr/dh` | relative permeability and its head derivative (NPF `KREL`, `DKRDH`) |
| `P` | precipitation rate, `[L/T]`, `P ≥ 0` (after `AUXMULTNAME` scaling) |
| `Qp` | `P·A`, `[L³/T]` |
| `d`, `d0` | pond depth at the end of the step / start of the step, `0 ≤ d ≤ dmax` |
| `dmax` | `MAX_POND_DEPTH`, static per cell |
| `Δt` | time step length |
| `Q_inf` | signed exchange with the aquifer through the top face: positive = infiltration (into the aquifer), negative = exfiltration |
| `Q_ro` | runoff, `≥ 0`, leaves the package |

### 4.2 Surface conductance

```
Cs = kr · K · A / (Δz / 2)
```

`Cs` is derived entirely from NPF; there is no user-specified surface or crust conductance.
This is a first-node approximation and depends on the top-cell thickness. Document it:
thick top cells under-resolve the capillary gradient at the surface.

### 4.3 Driving head

The surface-side head is the pond head:

```
Hs = z_s + d
Q_inf = Cs · (Hs − h)
```

A deeper pond pushes harder, and exfiltration stops when `h = Hs` (pond and soil are in
hydrostatic balance). Pond and exfiltration are therefore consistent.

### 4.4 Pond balance (implicit)

```
(d − d0) · A / Δt = Qp − Q_inf − Q_ro
```

with `0 ≤ d ≤ dmax`, `Q_ro ≥ 0`, and `Q_ro > 0` only when `d = dmax`.

### 4.5 Closed-form solution

With `Cs` frozen at the current iterate, the pond balance is linear in `d`, so no nested
iteration is required. Define:

```
β      = Δt · Cs / A
Hs0    = z_s + d0 + Δt · Qp / A
S      = Qp + d0 · A / Δt                  (water available this step)
Cs_eff = Cs / (1 + β)
```

Then, for transient steps:

```
Q_inf = min( S,
             Cs_eff · (Hs0 − h),           (pond as implicit series storage)
             Cs · (z_s + dmax − h) )       (pond full)
```

The three terms are the supply-limited, pond-storing, and pond-full regimes. Each is linear
in `h`, so `Q_inf(h)` is concave and piecewise linear.

Afterwards the state and runoff follow from mass balance:

```
d_unclamped = d0 + (Δt / A) · (Qp − Q_inf)
d     = clamp(d_unclamped, 0, dmax)
Q_ro  = Qp − Q_inf − (d − d0) · A / Δt       (≥ 0)
```

For steady-state stress periods the pond has no storage term:

```
Q_inf = min( Qp, Cs · (z_s + dmax − h) )
d     = clamp( (Qp − Cs·(z_s − h)) / Cs , 0, dmax )   (when Cs > 0)
Q_ro  = Qp − Q_inf
```

`dmax = 0` reduces to the no-pond case: all water that does not infiltrate immediately
becomes runoff.

### 4.6 Smoothing

The `min()` is smoothed with a fixed internal width (no user input in v1) so that the Newton
iteration converges. Invariants:

- The smoothed `Q_inf` used in the matrix is the same `Q_inf` used to compute `d`, `Q_ro`, and
  budget terms. Mass balance is then exact regardless of the smoothing.
- The smoothing never produces a `d` outside `[0, dmax]` (the clamp in §4.5 is applied after
  smoothing).
- The smoothed function has a continuous derivative with respect to `h`.

The width is an implementation constant to be tuned against the tests in §10.

### 4.7 Matrix terms

MODFLOW 6 sign convention: `Q_inf = hcof · h − rhs` (positive into the aquifer).

- `bnd_cf` selects the active (smoothed) branch from the current `h` and sets `hcof`, `rhs`.
  The branch conductance is `Cs`, `Cs_eff`, or 0; `kr` is frozen at the iterate. This is
  valid for Picard.
- `bnd_fn` adds the Newton term from `dkr/dh` (via `Cs`) using NPF's `DKRDH`. Pond depth is
  frozen when linearizing with respect to head (quasi-Newton); the outer iteration converges
  pond and head together.

### 4.8 Cells that are not active

If `ibound(n) ≤ 0` or `ibound(n) == IWETLAKE` (inactive, constant head, or lake-covered):

- `Q_inf = 0`, `d = 0`,
- all precipitation, plus any `d0 · A / Δt`, becomes runoff (budgeted and routable via MVR).

### 4.9 Relationship to empirical infiltration models

Horton, Green–Ampt, Philip, and SCS-type models prescribe an infiltration capacity `f_c(t)`
(or `f_c(F)`, with `F` the cumulative infiltration) from fitted parameters such as `f0`, `fc`,
and a decay constant, and compare it with the rainfall rate to split rainfall into infiltration
and infiltration-excess ("Hortonian") runoff. VSSB has no such parameters. The capacity is
`Cs · (Hs − h)`, which follows from the Richards solution in the top cell: it is high while the
near-surface soil is dry (steep gradient, small `kr`-limited conductance) and decays toward the
saturated vertical conductivity as the profile wets, so Horton-type behavior emerges rather than
being imposed. Runoff arises in the two ways these models distinguish: infiltration-excess
runoff when `P` exceeds the capacity and the pond is full, and saturation-excess ("Dunne")
runoff and exfiltration when the water table reaches the surface and `Cs · (Hs − h) → 0` or
changes sign. Empirical models cannot represent the second case, nor the interaction with a
rising water table, lateral flow in the soil, or root uptake. Two caveats apply. First, the
capacity is a first-node estimate (§4.2): the gradient is taken over `Δz / 2`, so a thick top
cell underestimates the very high early-time capacity (sorptivity-controlled, `f_c → ∞` as
`t → 0` in Green–Ampt) and the result depends on the vertical discretization of the top cell.
Second, `MAX_POND_DEPTH` is a storage/runoff threshold, not a depression-storage or
surface-sealing parameter; crusting must be represented through the soil properties of the top
cell. A fine vertical discretization of a uniform column with ponded surface is the intended way
to compare VSSB with Green–Ampt or Philip solutions as a verification exercise.

## 5. Mover

VSSB is an MVR provider only (`MOVER` option). The available water for cell `i` is `Q_ro`,
accumulated in `bnd_cf`/`bnd_fc` through `pakmvrobj`, as in DRN and UZF. The MVR provider id
is the `IFNO` of the cell. After MVR runs:

- `RUNOFF-TO-MVR` = amount taken by MVR,
- `RUNOFF` = remainder, which leaves the model unrouted.

Because the surface-cell set is fixed (§6), MVR ids are stable for the whole simulation.

## 6. Input

The cell set is fixed for the whole simulation, following UZF: `PACKAGEDATA` defines the cells
and `PERIOD` varies values by `IFNO`. This is why the set is fixed: a package with per-cell
state cannot drop cells between periods without losing water or breaking MVR ids.

### 6.1 Options

`AUXILIARY`, `AUXMULTNAME` (multiplier for `PRECIPITATION`), `BOUNDNAMES`, `PRINT_INPUT`,
`PRINT_FLOWS`, `SAVE_FLOWS`, `TS6`, `OBS6`, `MOVER`, and an optional pond-depth binary output
option that mirrors VSF's `PHEAD` output option naming.

### 6.2 DIMENSIONS

`NSURFCELLS` — number of surface cells (also sets `maxbound`).

### 6.3 PACKAGEDATA

| Field | Type | Notes |
|---|---|---|
| `IFNO` | integer | consecutive 1…`NSURFCELLS` |
| `CELLID` | integer(ncelldim) | surface cell |
| `MAX_POND_DEPTH` | double | static, `≥ 0`; `0` = no pond (default `0`) |
| `AUX` | double(naux) | optional |
| `BOUNDNAME` | string | optional |

### 6.4 PERIOD

`IFNO` followed by `PRECIPITATION` (time-series capable, `[L/T]`, `≥ 0`) and/or `AUX` values.
Values persist until changed. `PRECIPITATION` defaults to 0.

### 6.5 Validation (read time, once)

Errors:

- VSF package not active in the model,
- duplicate `CELLID` within the package,
- the same node listed in another VSSB package of the model,
- the same node listed in a VSSP package with `IHC = 0` (§7),
- an active cell directly above the listed cell (an active neighbor through a vertical
  connection, `ihc = 0`, whose bottom is at or above this cell's top); skipped where the grid
  has no vertical connections,
- `MAX_POND_DEPTH < 0`,
- `PRECIPITATION < 0` (checked when values are read, including time-series values each step),
- `IFNO` out of range or repeated in `PACKAGEDATA`.

## 7. VSSB/VSSP overlap

VSSB owns the top face of its cells, including exfiltration there. A VSSP entry with
`IHC = 0` on the same node would double-count that face, so it is an error. VSSP entries with
`IHC = 1` (lateral faces) on a VSSB cell are allowed.

Implementation note: the check needs both packages; implement it at the model level after all
packages are defined, not inside either package's input reading. `IHC = 0` does not say whether
the face is up or down; the conservative rule above treats every vertical VSSP face on a VSSB
cell as a conflict.

## 8. Budget, observations, output

### 8.1 Package budget (closed for each step)

```
PRECIPITATION = GWF + RUNOFF + RUNOFF-TO-MVR + STORAGE
```

| Term | Meaning |
|---|---|
| `PRECIPITATION` | `Qp` in |
| `GWF` | `Q_inf`; positive into the aquifer. Reported in/out by sign, so infiltration and exfiltration are both visible |
| `RUNOFF` | runoff not taken by MVR (leaves the model) |
| `RUNOFF-TO-MVR` | runoff taken by MVR |
| `STORAGE` | pond storage change rate, `(d − d0)·A/Δt` |

The model budget line for the package is `GWF` (the aquifer exchange). Precipitation, runoff,
and storage are package-internal terms.

### 8.2 Observations

`precipitation`, `infiltration` (`GWF` term, positive into the aquifer), `runoff`,
`runoff-to-mvr`, `pond-storage` (storage rate), `pond-depth`. `pond-depth` is a state, not a
flow. Boundnames and cell-based ids follow the existing package conventions.

### 8.3 Output

Optional binary pond-depth file, like VSF's `PHEAD`.

## 9. Time stepping and state

- `bnd_ad`: commit `d0 ← d` from the previous step. Pond depth starts at 0.
- `bnd_cf` (each outer iteration): compute `Cs`, `Q_inf`, `d`, `Q_ro` from the current `h` and
  `d0`; set `hcof`/`rhs`.
- Package convergence: the outer convergence check also includes the pond-depth change
  between iterations.
- `bnd_cq`/`bnd_bd`: compute final flows and the budget from the converged `h`.
- No pond state is saved or restored across runs in v1.

## 10. Tests (integration only, `autotest/`)

1. **Under capacity.** 1D column, `P` below `Cs`-limited capacity: no pond, no runoff;
   head profile matches the same flux applied with RCH.
2. **Over capacity.** 1D column, `P` above capacity: pond fills to `dmax`, then runoff equals
   `Qp − Q_inf` with exact mass balance; with `dmax = 0` all excess is runoff immediately.
3. **MVR.** Runoff routed to a receiver (e.g. a LAK or SFR reach); `RUNOFF-TO-MVR` equals the
   receiver inflow; an unrouted remainder appears in `RUNOFF`.
4. **Exfiltration.** Water table above the surface (e.g. CHD at depth above `z_s`): pond
   fills first, then runs off; steady state satisfies `Q_inf = Cs · (z_s + dmax − h)`.
5. **Rain, then dry.** A rainfall pulse followed by `P = 0`: the pond drains by infiltration,
   no chatter, bounded iteration counts.
6. **Newton and Picard.** The over-capacity case runs under both and agrees.
7. **Steady-state period.** The steady-state branch of §4.5 matches the transient solution at
   long times.
8. **Inactive and constant-head cells.** Precipitation becomes direct runoff.
9. **Validation errors,** where the framework supports expected failures: missing VSF,
   duplicate cell, active cell above, VSSP overlap, negative precipitation.

## 11. Implementation checklist

- New package source `gwf-vssb.f90` (`BndExtType`, UZF-style `PACKAGEDATA`/`PERIOD` reading),
  `vssb_create` factory.
- Register `VSSB6`: `gwf.f90` (package list and constructor case), `IdmGwfDfnSelector`, and
  `src/meson.build`; add `doc/mf6io/mf6ivar/dfn/gwf-vssb.dfn`. Regenerate IDM Fortran and
  mf6io input docs (`pixi run update-fortran-definitions`, `pixi run run-mf6ivar`); do not
  hand-edit generated files.
- Seepage/exchange kernel as a standalone module (smoothed, Newton-linearized, `kr`-scaled
  conductance), so VSSP can adopt it in the follow-up.
- Model-level VSSB/VSSP overlap check (§7).
- Docs: mf6io input description, release notes entry, a note that VSSB requires VSF, and a note
  that evaporation from ponded water is not modelled.
- Integration tests per §10; format Fortran with `fprettify` and Python with `ruff`.

## 12. Open details to settle during implementation

- Smoothing width and function for the `min()`.
- Treatment of `Cs` when `kr` is extremely small (floor on `Cs` for the closed form, and
  `β` overflow handling when `Cs` is large).
- Exact observation and budget label spellings.
- Whether `RUNOFF`/`RUNOFF-TO-MVR` accounting after MVR needs a one-iteration lag, as in UZF.
- Whether to expose pond depth through BMI in v1.
