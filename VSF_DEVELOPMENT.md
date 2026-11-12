# MODFLOW 6 VSF Development Specifications

Consolidated specifications and design notes for the Variably Saturated Flow (VSF) family of MODFLOW 6 packages. Each chapter was previously a standalone design document in the repository root; those documents have been consolidated into this file, which is now the single editable source.

## Table of contents

1. [VSSB: Variably Saturated Surface Boundary package — v1 specification](#vssb-variably-saturated-surface-boundary-package--v1-specification) — formerly `VSSB_SPEC.md`
2. [Design note: soil and pond evaporation in VSSB](#design-note-soil-and-pond-evaporation-in-vssb) — formerly `VSSB-EVAPORATION-DESIGN.md`
   - [Potential-rate partitioning: basis and choices](#7-potential-rate-partitioning-basis-and-choices) — literature and partitioning decisions
3. [VSSP (Seepage) Package — Status Handoff](#vssp-seepage-package--status-handoff) — formerly `VSSP_SEEPAGE_NOTES.md`
4. [Uncapped river leakage for Richards flow](#uncapped-river-leakage-for-richards-flow) — formerly `RIV_RICHARDS_LEAKAGE.md`
5. [Feddes transpiration for MODFLOW 6 VSF](#feddes-transpiration-for-modflow-6-vsf) — formerly `feddes-transpiration-proposal.md`
6. [Alternatives (not yet implemented)](#alternatives-not-yet-implemented) — superseded or deferred approaches
   - [SPG Seepage Package — Notes and Findings](#spg-seepage-package--notes-and-findings) — formerly `OLD_SPG_SEEPAGE_NOTES.md`
   - [VS2D-compatible transpiration for MODFLOW 6 VSF](#vs2d-compatible-transpiration-for-modflow-6-vsf) — formerly `vs2d-transpiration-proposal.md`

---

## VSSB: Variably Saturated Surface Boundary package — v1 specification

Status: draft for review. Supersedes the earlier concept note (`VSSB_concept.md`).
Package type: `VSSB6` (GWF model). Requires the `VSF6` package in the same model.

### 1. Purpose

VSSB is the land-surface boundary for VSF (Richards) models. It applies precipitation to
the top face of a set of surface cells, limits infiltration by the soil, stores excess
water in a local pond, lets groundwater exfiltrate into that pond, and exports water that
cannot be stored as runoff through the MODFLOW 6 water mover (MVR).

Recharge, ponding, runoff, and exfiltration are results of the solution; none of them is
prescribed.

### 2. Scope

In scope for v1:

- precipitation on a fixed set of top-face cells,
- infiltration limited by the soil (`kr`-dependent conductance), driven by total head,
- exfiltration of groundwater through the same face (replaces the surface role of VSSP),
- per-cell pond storage with a static maximum depth,
- runoff as an MVR provider,
- evaporation from ponded water and bare soil (potential `EVAPORATION` rate per surface cell;
  implemented, see the design note below),
- Newton and Picard formulations,
- budget, observations, and optional pond-depth output.

Out of scope for v1 (documented gaps or later work):

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

### 3. Relationship to other packages

| Package | Relationship |
|---|---|
| VSF | Required. VSSB reads `KREL`/`DKRDH`, K and geometry from NPF, which VSF populates. VSF currently rejects mixed domains (`IVARSAT = 0` is fatal), so VSSB assumes every cell is a Richards cell. |
| VSET | Independent. Root uptake (transpiration) stays in VSET. Evaporation from ponded water and bare soil is handled by VSSB; the potential `EVAPORATION` (VSSB) and potential transpiration (VSET) are independent and additive. |
| VSSP | Remains for lateral/bottom seepage faces. A node may not be listed in both VSSB and VSSP when the VSSP entry has `IHC = 0` (see §7). |
| RCH, EVT, UZF | Not replaced. VSSB is the Richards-domain surface boundary. |
| MVR | VSSB is a provider only. Runoff is the provided water. |
| Future overland-flow model | Receives runoff through MVR (or an equivalent exchange). VSSB does not route water laterally. |

### 4. Formulation

#### 4.1 Notation (per surface cell *i*, node *n*)

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

#### 4.2 Surface conductance

```
Cs = kr · K · A / (Δz / 2)
```

`Cs` is derived entirely from NPF; there is no user-specified surface or crust conductance.
This is a first-node approximation and depends on the top-cell thickness. Document it:
thick top cells under-resolve the capillary gradient at the surface.

#### 4.3 Driving head

The surface-side head is the pond head:

```
Hs = z_s + d
Q_inf = Cs · (Hs − h)
```

A deeper pond pushes harder, and exfiltration stops when `h = Hs` (pond and soil are in
hydrostatic balance). Pond and exfiltration are therefore consistent.

#### 4.4 Pond balance (implicit)

```
(d − d0) · A / Δt = Qp − Q_inf − Q_ro
```

with `0 ≤ d ≤ dmax`, `Q_ro ≥ 0`, and `Q_ro > 0` only when `d = dmax`.

#### 4.5 Closed-form solution

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

#### 4.6 Smoothing

The `min()` is smoothed with a fixed internal width (no user input in v1) so that the Newton
iteration converges. Invariants:

- The smoothed `Q_inf` used in the matrix is the same `Q_inf` used to compute `d`, `Q_ro`, and
  budget terms. Mass balance is then exact regardless of the smoothing.
- The smoothing never produces a `d` outside `[0, dmax]` (the clamp in §4.5 is applied after
  smoothing).
- The smoothed function has a continuous derivative with respect to `h`.

The width is an implementation constant to be tuned against the tests in §10.

#### 4.7 Matrix terms

MODFLOW 6 sign convention: `Q_inf = hcof · h − rhs` (positive into the aquifer).

- `bnd_cf` selects the active (smoothed) branch from the current `h` and sets `hcof`, `rhs`.
  The branch conductance is `Cs`, `Cs_eff`, or 0; `kr` is frozen at the iterate. This is
  valid for Picard.
- `bnd_fn` adds the Newton term from `dkr/dh` (via `Cs`) using NPF's `DKRDH`. Pond depth is
  frozen when linearizing with respect to head (quasi-Newton); the outer iteration converges
  pond and head together.

#### 4.8 Cells that are not active

If `ibound(n) ≤ 0` or `ibound(n) == IWETLAKE` (inactive, constant head, or lake-covered):

- `Q_inf = 0`, `d = 0`,
- all precipitation, plus any `d0 · A / Δt`, becomes runoff (budgeted and routable via MVR).

#### 4.9 Relationship to empirical infiltration models

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

### 5. Mover

VSSB is an MVR provider only (`MOVER` option). The available water for cell `i` is `Q_ro`,
accumulated in `bnd_cf`/`bnd_fc` through `pakmvrobj`, as in DRN and UZF. The MVR provider id
is the `IFNO` of the cell. After MVR runs:

- `RUNOFF-TO-MVR` = amount taken by MVR,
- `RUNOFF` = remainder, which leaves the model unrouted.

Because the surface-cell set is fixed (§6), MVR ids are stable for the whole simulation.

### 6. Input

The cell set is fixed for the whole simulation, following UZF: `PACKAGEDATA` defines the cells
and `PERIOD` varies values by `IFNO`. This is why the set is fixed: a package with per-cell
state cannot drop cells between periods without losing water or breaking MVR ids.

#### 6.1 Options

`AUXILIARY`, `AUXMULTNAME` (multiplier for `PRECIPITATION`), `BOUNDNAMES`, `PRINT_INPUT`,
`PRINT_FLOWS`, `SAVE_FLOWS`, `TS6`, `OBS6`, `MOVER`, and an optional pond-depth binary output
option that mirrors VSF's `PHEAD` output option naming.

#### 6.2 DIMENSIONS

`NSURFCELLS` — number of surface cells (also sets `maxbound`).

#### 6.3 PACKAGEDATA

| Field | Type | Notes |
|---|---|---|
| `IFNO` | integer | consecutive 1…`NSURFCELLS` |
| `CELLID` | integer(ncelldim) | surface cell |
| `MAX_POND_DEPTH` | double | static, `≥ 0`; `0` = no pond (default `0`) |
| `AUX` | double(naux) | optional |
| `BOUNDNAME` | string | optional |

#### 6.4 PERIOD

`IFNO` followed by `PRECIPITATION` (time-series capable, `[L/T]`, `≥ 0`) and/or `AUX` values.
Values persist until changed. `PRECIPITATION` defaults to 0.

#### 6.5 Validation (read time, once)

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

### 7. VSSB/VSSP overlap

VSSB owns the top face of its cells, including exfiltration there. A VSSP entry with
`IHC = 0` on the same node would double-count that face, so it is an error. VSSP entries with
`IHC = 1` (lateral faces) on a VSSB cell are allowed.

Implementation note: the check needs both packages; implement it at the model level after all
packages are defined, not inside either package's input reading. `IHC = 0` does not say whether
the face is up or down; the conservative rule above treats every vertical VSSP face on a VSSB
cell as a conflict.

### 8. Budget, observations, output

#### 8.1 Package budget (closed for each step)

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

#### 8.2 Observations

`precipitation`, `infiltration` (`GWF` term, positive into the aquifer), `runoff`,
`runoff-to-mvr`, `pond-storage` (storage rate), `pond-depth`. `pond-depth` is a state, not a
flow. Boundnames and cell-based ids follow the existing package conventions.

#### 8.3 Output

Optional binary pond-depth file, like VSF's `PHEAD`.

### 9. Time stepping and state

- `bnd_ad`: commit `d0 ← d` from the previous step. Pond depth starts at 0.
- `bnd_cf` (each outer iteration): compute `Cs`, `Q_inf`, `d`, `Q_ro` from the current `h` and
  `d0`; set `hcof`/`rhs`.
- Package convergence: the outer convergence check also includes the pond-depth change
  between iterations.
- `bnd_cq`/`bnd_bd`: compute final flows and the budget from the converged `h`.
- No pond state is saved or restored across runs in v1.

### 10. Tests (integration only, `autotest/`)

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

### 11. Implementation checklist

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
  that evaporation (from ponded water and bare soil) is simulated when a potential `EVAPORATION`
  rate is specified.
- Integration tests per §10; format Fortran with `fprettify` and Python with `ruff`.

### 12. Open details to settle during implementation

- Smoothing width and function for the `min()`.
- Treatment of `Cs` when `kr` is extremely small (floor on `Cs` for the closed form, and
  `β` overflow handling when `Cs` is large).
- Exact observation and budget label spellings.
- Whether `RUNOFF`/`RUNOFF-TO-MVR` accounting after MVR needs a one-iteration lag, as in UZF.
- Whether to expose pond depth through BMI in v1.

---

## Design note: soil and pond evaporation in VSSB

Status: implemented. The evaporation design in Sections 2-4 is implemented in the
VSSB Package and its surface-exchange kernel. The seepage-face slope extension in
Section 6 remains out of scope.

This note records the design for adding evaporation (from bare soil and from
ponded water) to the Variably Saturated Surface Boundary (VSSB) Package, the
relationship between VSSB and the seepage-face package (VSSP), and the
complexity of a possible future extension to seepage faces under a slope that is
not aligned with the grid. The slope extension is discussed for scoping only and
is explicitly out of scope for implementation.

### 1. Current architecture

The variably-saturated family is a core Richards flow package plus boundary
packages:

- `src/Model/GroundWaterFlow/gwf-vsf.f90` — the Richards (unsaturated) flow
  package. It flags Richards cells (`IVARSAT`) and publishes relative
  permeability `KREL` and its head derivative `DKRDH` to NPF, which the boundary
  packages read.
- `src/Model/GroundWaterFlow/gwf-vssb.f90` — the land-surface boundary:
  precipitation to infiltration / ponding / exfiltration / runoff on the **top
  face** of surface cells.
- `src/Model/GroundWaterFlow/gwf-vssp.f90` — a seepage-face boundary on an
  **arbitrary face** (`IHC=0` vertical, `IHC=1` horizontal) that discharges
  groundwater when head rises above the face.
- `src/Model/GroundWaterFlow/gwf-vset.f90` — variably-saturated **transpiration**
  (root-water uptake, Feddes). A pressure-head-dependent sink on a cell; it is
  the closest existing template for an evaporation sink, but it models
  transpiration only.

The physics heart of VSSB is the package-independent kernel
`surface_exchange` in `src/Model/ModelUtilities/VsfSurfaceExchange.f90`, called
from `vssb_evaluate`. It solves an implicit pond balance in closed form as a
smoothed minimum of three terms (supply-limited, pond-storing, pond-full) and
returns `qinf`, `pond`, `qro`, and the linearization (`hcof`, `rhs`, Newton
`deriv`). The transient pond balance it closes is:

```
(pond - pond0) * area / dt = qp - qinf - qro
```

This note originally recorded the design for the gap that the pre-evaporation
documentation named: *"Evaporation from ponded water is not simulated by the VSSB
Package or by any other package."* That gap is now closed — `doc/mf6io/gwf/vssb.tex`
documents the implemented evaporation behavior described below.

### 2. Evaporation demand and wiring (decision: Option A)

#### 2.1 Two separate partitioning levels

Only the first level is a VSSB-versus-VSET question:

```
Atmospheric PET
  |-- potential SOIL+POND evaporation  E_pot  -> VSSB   (new)
  |      |-- pond (open-water) evaporation       (internal to VSSB)
  |      '-- bare-soil evaporation               (internal to VSSB)
  '-- potential TRANSPIRATION           T_pot  -> VSET   (exists today, = rate)
```

- **Level 1 (E_pot vs T_pot):** the split of atmospheric demand into
  bare-soil/open-water evaporation (VSSB) and plant transpiration (VSET).
- **Level 2 (pond vs soil):** internal to VSSB. Open water evaporates first at
  the potential rate, and only the residual demand,
  `max(E_pot - E_pond, 0)`, is exposed to the drying soil. VSET is not involved.

#### 2.2 How VSET defines its demand today

VSET does not read a shared evapotranspiration demand. It reads a per-cell
`rate` in its own PERIOD block (`gwf-vset.dfn`), time-series-enabled, documented
as the potential transpiration flux rate for the cell. That value flows into
`vset_cf` as `this%rate(i)`, is multiplied by cell area, and is reduced to
actual uptake by the Feddes response. There is no "ET demand bus" to tap into.

#### 2.3 Decision: Option A — VSSB owns `E_pot`

VSSB defines its own potential evaporation rate per surface cell, mirroring
VSET's `rate` and VSSB's existing `PRECIPITATION`. Its documented meaning is
**potential soil + pond evaporation only**.

- There is **no code-level coupling** between VSSB and VSET. The two demands are
  additive sinks.
- The Level-1 split (PET into E_pot and T_pot, e.g. via a cover fraction /
  Beer's-law LAI relation) is the responsibility of the user or the preprocessor
  (flopy), supplied as two independent input series.
- This matches the MODFLOW 6 convention that each stress package reads its own
  time-series rate and packages stay loosely coupled.

Consequences to document and test:

- Nothing enforces `E_pot + T_pot <= PET`; a mis-split can let total actual ET
  exceed PET. This is the same trust model MODFLOW already uses for separate RCH
  / EVT packages, but it must be stated plainly in `doc/mf6io/gwf/vssb.tex`.
- **Co-located cells:** if a VSSB surface cell is also a VSET cell, both sinks
  act on the same top cell. They couple through head in the matrix, but each
  caps itself at its own potential, so the combined draw is bounded by
  `E_pot + T_pot` — correct provided the split was done upstream. This warrants
  an explicit autotest.

A future single-conserved-PET scheme (previously called Option B) can be layered
on later using the existing not-owned memory-pointer mechanism (the same pattern
VSSB already uses to read NPF `KREL`/`DKRDH` and VSF `IVARSAT`), without changing
VSSB's input surface. Choosing Option A now does not foreclose it.

#### 2.4 Decisions

- **Input:** `EVAPORATION` is a per-surface-cell, time-series-enabled potential
  evaporation rate (`L/T`) in the VSSB PERIOD block, always specified
  (positional, following `PRECIPITATION`), mirroring UZF's `PET`. A value of
  zero reproduces the current no-evaporation behavior exactly.
- **Soil drying reduction model:** reuse the VSF relative permeability `KREL`
  already published to NPF; actual bare-soil evaporation is `KREL * E_soil_pot`.
  No new input and no new option. (A pressure-head ramp was considered and
  rejected for now.)
- **Budget / observations:** a single combined `EVAPORATION` term (pond + soil),
  not split.

### 3. Evaporation formulation in VSSB

#### 3.1 The two sinks

1. **Pond (open-water) evaporation** — removes water from pond storage at (close
   to) the potential rate while the pond holds water; limited by available pond
   water. It is a sink on the pond balance only; it reaches the GWF matrix only
   indirectly, through its effect on pond depth (and hence on
   infiltration/exfiltration in the next iteration).
2. **Bare-soil evaporation** — removes water directly from the top cell (the
   aquifer), limited both by the potential rate and by how dry the soil is. It
   is a head-dependent sink on the GWF cell, structurally like VSET's uptake
   term (`hcof`/`rhs` plus a Newton term from the `KREL`/pressure-head
   dependence), but keyed to soil-surface drying rather than a root-zone Feddes
   curve.

#### 3.2 Demand partitioning (Level 2)

Potential evaporation `E_pot` is one demand. Open water evaporates first at the
potential rate; only the residual `max(E_pot - E_pond, 0)` is exposed to the
soil, where it is further reduced by soil dryness. This keeps total evaporation
`<= E_pot` and avoids double counting. Because both the pond term and the
residual are capacity-limited, they fit the kernel's existing smoothed-minimum
machinery. The cleanest implementation extends `surface_exchange` to carry the
evaporation demand and emit the extra sink terms, rather than bolting
evaporation on outside the kernel (which would break the exact pond-balance
closure and the Newton consistency).

#### 3.3 Linearization / convergence

The soil-evaporation sink needs `hcof`/`rhs` contributions and a Newton term for
its dependence on head through soil wetness, mirroring `vssb_fn` and the `DKRDH`
handling already in the kernel. Pond evaporation enters the pond update and the
runoff residual; its head dependence is weak, but the smoothed limiters must
stay C1 to preserve the current convergence behavior.

#### 3.4 Interaction with exfiltration

Exfiltrated groundwater lands in the pond, so pond evaporation can consume
exfiltrated water. This is physically correct and the budget still closes, but
it means the GWF and EVAPORATION terms interact; the test suite should assert
the closed water balance explicitly.

#### 3.5 Derived kernel equations (as implemented)

Let `E_pot` be the potential evaporation (volumetric rate, `>= 0`) for the cell,
`qp` the precipitation supply, `d0` the start-of-step pond depth, `krel`,
`dkrdh` the relative permeability and its head derivative, and `cs`, `hfull`,
etc. as in the existing kernel.

Pond evaporation is limited by the water available in the pond this step and is
independent of head, so it is a constant during the Newton iterations:

```
supply_evap = qp                       (steady, no pond storage)
supply_evap = qp + d0 * area / dt       (transient)
E_pond      = min(E_pot, max(supply_evap, 0))
qp_eff      = qp - E_pond
```

The existing infiltration / pond / runoff kernel is then run unchanged **with
`qp_eff` in place of `qp`** (in `supply`, `hs0`, the steady pond diagnostic, and
the runoff residual). Because the kernel conserves `qp_eff + d0*area/dt` among
infiltration, runoff and storage, total mass is conserved:
`E_pond + qinf + qro + storage = qp + d0*area/dt`. Since `supply_evap` is a
start-of-step constant, `qp_eff` is constant during the solve and the existing
`hcof`/`rhs`/`deriv` are structurally unchanged.

The residual atmospheric demand is exposed to the drying soil, reduced by `krel`
(the chosen soil-drying model):

```
E_soil_pot = max(E_pot - E_pond, 0)
E_soil     = krel * E_soil_pot
```

`E_soil` is a head-dependent cell sink. Using the VSSB/VSET sign convention
`q_into_cell = hcof*h - rhs`, the soil-evaporation flow into the cell is
`q_soil = -E_soil = -krel * E_soil_pot`, giving `hcof_soil = 0`,
`rhs_soil = krel * E_soil_pot` (frozen `krel`, Picard), and a Newton term
`deriv_soil = d(q_soil)/dh = -E_soil_pot * dkrdh`. These are added to the pond
terms:

```
res%hcof  = hcof_pond                       (soil adds 0)
res%deriv = deriv_pond - E_soil_pot * dkrdh
res%rhs   = res%hcof * h - (qinf_pond - E_soil)   (total flow into cell)
```

Total evaporation never exceeds potential:
`E_pond + E_soil = E_pond*(1 - krel) + krel*E_pot <= E_pot` for `krel <= 1`.

**Budget closure (combined `EVAPORATION` term).** `E_pond` comes from the pond
(surface water); `E_soil` is extracted from the aquifer, so the aquifer must be
credited for it. With flows into the package positive:

```
PRECIPITATION = qp
GWF           = -(qinf_pond - E_soil)   (includes the soil-evaporation extraction)
RUNOFF        = -(qro - qrotomvr)
STORAGE       = -qsto
EVAPORATION   = -(E_pond + E_soil)      (combined, to atmosphere)
```

Substituting the pond balance `qsto = qp - qinf_pond - qro - E_pond` shows the
package budget closes exactly, and the `GWF` term matches the matrix flow into
the cell (`qinf_pond - E_soil`). The `infiltration` observation remains the pond
exchange `qinf_pond` only.

**Inactive / lake-covered cells.** The existing branch routes all water to
runoff and performs no aquifer exchange; evaporation is left at zero there (the
atmospheric demand is simply unmet), consistent with that simplification.

### 4. Code touch points

- **Input definition:** add a period variable (`EVAPORATION`,
  time-series-enabled, like `PRECIPITATION`, always specified) in
  `doc/mf6io/mf6ivar/dfn/gwf-vssb.dfn`;
  regenerate the IDM Fortran and mf6io docs
  (`pixi run update-fortran-definitions`, `pixi run run-mf6ivar`).
- **Type + allocation:** new arrays (`evap` input pointer, `qevap_pond`,
  `qevap_soil`) in `VssbType` and `vssb_allocate_vars`.
- **Read / validate:** source the new period data in `vssb_rp`; add a
  non-negativity check in `vssb_ad` alongside the precipitation check.
- **Kernel:** extend `SurfaceExchangeResultType` and `surface_exchange` to
  accept `E_pot` and return the two evaporation flows and their linearization;
  thread it through `vssb_evaluate`, `vssb_cf`, `vssb_fn`, `vssb_cq`.
- **Budget:** add a single combined `EVAPORATION` term in `vssb_setup_budobj`
  and `vssb_fill_budobj`, and fold the soil-evaporation extraction into the
  `GWF` term; update the documented balance to
  `PRECIPITATION = GWF + RUNOFF + RUNOFF-TO-MVR + STORAGE + EVAPORATION`.
- **Observations:** add an `evaporation` obs type in the `bnd_*_obs` routines and
  `doc/Common/gwf-vssbobs.tex`.
- **Docs + tests:** update `doc/mf6io/gwf/vssb.tex` (remove the "not simulated"
  sentence, add the Option A wiring explanation), add a release note, and add
  autotests covering pond-only, soil-only, combined, co-located VSSB/VSET, and
  exfiltration-then-evaporate water-balance closure.

### 5. Relationship to VSSP (exfiltration)

The overlap between VSSB and VSSP is about **exfiltration / seepage**.

- VSSB already performs seepage out of the **top vertical face** of a surface
  cell: when cell head exceeds the pond head, `qinf` becomes negative
  (exfiltration into the pond), where the water can run off, be moved by MVR,
  or — once this work lands — evaporate. This is why the current rule forbids a
  VSSB cell from also being a VSSP cell with `IHC=0`
  (`doc/mf6io/gwf/vssb.tex`).
- VSSP's remaining, non-overlapping role is seepage on faces VSSB does not own:
  - **Horizontal connections (`IHC=1`)** — vertical faces, e.g. the walls of a
    hillslope, a trench, or a lateral domain boundary.
  - **Vertical faces on cells that are not VSSB surface cells** (not at the top
    of the active domain).

**Recommendation:** keep VSSP as the general seepage-face package; do not fold
general seepage faces into VSSB. The clean division of labor is VSSB = the top
land surface (precip / pond / infiltration / exfiltration / runoff, and now
evaporation); VSSP = seepage on arbitrary non-top faces. A smaller, optional
cleanup is to retire only the VSSP `IHC=0`-at-the-top use case in favor of VSSB,
since that is the genuinely duplicated path.

### 6. Future extension (out of scope): seepage faces under a slope

A seepage face that cuts obliquely through cells (a slope not aligned with the
grid) belongs conceptually with VSSP, not VSSB, and is considerably more
involved than adding evaporation. It is documented here only to scope the
complexity; it must not be implemented as part of the evaporation work.

1. **Geometry definition.** A grid-aligned seepage face is described by `IHC`,
   `CL1`, `HWVA`, `ANGLDEGX`. A sloped face needs intersection geometry: the
   sloped face area within each intersected cell, the outward normal, and the
   distance from the cell center to the oblique face. None of this is a native
   grid quantity; it must be supplied or computed from a surface/slope
   description — a substantial new input-geometry problem.
2. **Conductance along an oblique normal.** VSSP currently uses a LAK-style
   approximation that reads grid-aligned `K11`/`K33` directly and ignores
   `ANGLDEGX` and horizontal anisotropy (`effective_k`). A sloped face has a
   mixed horizontal/vertical normal, so a correct conductance must project the
   full anisotropic tensor onto the oblique normal (`n . K . n`), resolving
   `K11`/`K22`/`K33` plus `ANGLDEGX`/dip — a formulation upgrade, not just input.
3. **Seepage elevation datum.** The pressure-head switch (`psi = head - z`)
   currently uses a single face-centroid `z`. On a slope, the seepage elevation
   varies across the face within a cell; choosing a representative datum
   (centroid vs sub-face discretization) affects accuracy and Newton behavior.
4. **Where the discharged water goes.** Grid-aligned seepage currently just
   leaves (MVR is not yet supported — see the TODO in `vssp_fc`). Slope seepage
   realistically wants to route downslope to the next cell or into a routing
   network (MVR / overland flow), entangling this extension with finishing mover
   support and possibly surface routing.
5. **Multiple faces per cell and bookkeeping.** An oblique slope can intersect
   several cells and give a cell more than one seepage face; the current
   one-entry-per-cell list model and the "at most one of VSSB / VSSP(`IHC=0`)
   per cell" constraint would need rethinking.
6. **Numerics.** Each oblique face adds a head-dependent, switch-activated term
   (dry vs seeping) that needs smoothing for a continuous Jacobian; VSSP's
   Newton routine is presently a stub (`vssp_fn` just returns), so Newton
   support would have to be built out for robust convergence.

Keeping VSSP as a distinct package (Section 5) is what makes this future work
tractable, because it isolates the oblique-face geometry from VSSB's top-face
assumptions.

#### A note on the 3D geometry

The geometry is the hardest part of this extension, and it is genuinely complex
in three dimensions. In 2D (a cross-section) a slope is just a polyline cutting
through a stack of cells, and each intersection is a line segment — tractable to
compute and reason about. In 3D the seepage surface becomes a two-dimensional
surface (a terrain / triangulated mesh) that must be intersected with a
three-dimensional grid of cells, and the difficulty jumps:

- **Polygon clipping, not segment intersection.** Each cell–surface overlap is a
  general polygon (sometimes non-convex, sometimes several disconnected pieces)
  that must be clipped against the cell faces to get the wetted face area and
  centroid. This is full computational-geometry work, with the usual robustness
  hazards: slivers, near-tangent intersections, surfaces passing exactly through
  a cell edge or vertex, and floating-point degeneracies.
- **Arbitrary grid types.** It must work for DIS, DISV, and DISU. DISV/DISU
  cells are general polyhedra, so there is no shortcut from a rectangular cell
  shape; the clipping has to handle arbitrary cell geometries.
- **Per-piece conductance geometry.** The outward normal and the cell-center-to-
  face distance (`CL1`) vary across the surface and must be derived per
  intersected piece, then fed into the oblique-normal tensor projection
  (point 2). A single representative value per cell is an approximation whose
  error is hard to bound on a curved or steep surface.
- **Connectivity for routing.** To route seepage downslope (point 4) the
  intersected pieces must be ordered by elevation and linked neighbor-to-
  neighbor along the surface — a topology that the grid connectivity does not
  provide and that must be reconstructed from the intersection.

In short, 2D slope seepage is a modest extension, but the 3D case is a
non-trivial geometry engine (surface–grid intersection, polygon clipping,
robustness handling) that is realistically a project in its own right, separate
from and much larger than the evaporation work.

### 7. Potential-rate partitioning: basis and choices

This section records *why* the evaporation partitioning in Sections 2–3 is set
up the way it is, and which parts are forced by physics versus deliberate
modeling choices. It is the literature anchor for reviewers deciding whether to
keep the current scheme or extend it.

#### 7.1 The partitioning is two distinct problems, not one ordering

The atmospheric demand is split at two levels with **different physical bases**;
conflating them is the usual source of confusion.

```
Atmospheric PET
  ├── Level 1 (energy / canopy cover): split PET into
  │     ├── E_pot  — bare-soil + open-water evaporation   → VSSB
  │     └── T_pot  — plant transpiration                  → VSET (Feddes)
  └── Level 2 (area / supply, internal to VSSB): split E_pot into
        ├── pond (open-water) evaporation   — from pond storage
        └── bare-soil evaporation           — from the top cell
```

The key reframing for Level 2: **pond vs. soil is not a time ordering** ("pond
first, then soil"); it is a **spatial competition for the same evaporative
energy over the same cell footprint**. Where the surface is ponded it evaporates
as open water; where it is bare it evaporates as soil. The current
"pond-first, residual to soil" formulation (§3.2, §3.5) is the limiting case of
an area weighting in which the pond, while present, is treated as covering the
full cell (see §7.4).

#### 7.2 Literature by sub-problem

**Level-1 (E vs. T) partition — canopy / LAI.**

- Ritchie, J.T. (1972), *Model for predicting evaporation from a row crop with
  incomplete cover*, Water Resources Research 8(5):1204–1213. The canonical
  Beer's-law split `E_pot = PET·exp(−β·LAI)`, `T_pot = PET·(1 − exp(−β·LAI))`;
  used by SWAP, DSSAT, and most agro-hydrology models.
- Shuttleworth, W.J. & Wallace, J.S. (1985), *Evaporation from sparse crops — an
  energy combination theory*, Q.J.R. Meteorol. Soc. 111:839–855. The physically
  rigorous dual-source (soil + canopy) model with coupled resistances, if more
  than Beer's law is ever needed.
- Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998), *Crop
  evapotranspiration*, FAO Irrigation and Drainage Paper 56. The operational
  **dual crop coefficient** `Kc = Kcb (transpiration) + Ke (evaporation)`; maps
  directly onto the Option A choice (§2.3) of two independent input series.

**Level-2 bare-soil evaporation reduction (the soil-drying limiter).** Three
families exist:

| Family | Mechanism | Example codes / refs | Fit for a Richards solver |
|---|---|---|---|
| (a) Hydraulic / Darcy-limited | soil supplies only what the unsaturated flux allows, via `krel` | **VSSB (this work)** | native — reuses the flow kernel's `krel` |
| (b) Critical surface pressure head | cap `E` when surface `ψ` would drop below atmospheric equilibrium `h_crit` | HYDRUS (Šimůnek et al.); VS2DT (Lappala, Wallace & Healy, 1987) | compatible, head-based |
| (c) Empirical two-stage (cumulative) | energy-limited stage 1 → `√t` / falling-rate stage 2 | Black, Gardner & Thurtell (1969); Boesten & Stroosnijder (1986); Ritchie (1972) — the SWAP options | adds a second, calibrated drying law |

VSSB uses **(a)**: `E_soil = krel · E_soil_pot`. It is the physically consistent
choice because `krel` is the *same* relative permeability already driving
infiltration/exfiltration (NPF `KREL`/`DKRDH`), so no new empiricism and no new
input are introduced (§2.4). It is equivalent in spirit to the head-based cap in
HYDRUS/VS2DT.

**Open-water / pond evaporation.**

- Penman, H.L. (1948), *Natural evaporation from open water, bare soil and
  grass*, Proc. R. Soc. London A 193:120–145.
- Priestley, C.H.B. & Taylor, R.J. (1972), *On the assessment of surface heat
  flux and evaporation using large-scale parameters*, Mon. Weather Rev.
  100:81–92.

The robust fact these establish: **open water is energy-limited, not
supply-limited** — it evaporates at (near) the potential rate until it is gone.
This is what justifies giving the pond the potential rate before the soil
(§3.5), independent of the ordering question.

**MODFLOW-specific conventions (for user-expectation consistency).**

- Niswonger, R.G., Prudic, D.E. & Regan, R.S. (2006), *UZF1 Package*, USGS
  Techniques and Methods 6-A19 — `PET`, extinction depth, extinction water
  content.
- Shah, N., Nachabe, M. & Ross, M. (2007), *Extinction depth and
  evapotranspiration from ground water under selected land cover classes*,
  Ground Water 45(3):329–338 — bare-soil ET vs. water-table / extinction depth;
  the usual reference for groundwater-ET parameterization in the MODFLOW
  community.

#### 7.3 What is forced vs. what is a choice

- **Forced by physics (keep):** pond at the potential (energy-limited) rate;
  soil reduced by a drying limiter; total actual evaporation `≤ E_pot`;
  `krel`-based limiting is the internally consistent limiter for a Richards
  kernel. These are not ambiguous.
- **A deliberate choice (Option A, §2.3):** Level-1 (E/T) partition is pushed
  **upstream** to the user / preprocessor as two independent time series, with
  Ritchie (1972) / FAO-56 dual-`Kc` as the recommended basis. Nothing enforces
  `E_pot + T_pot ≤ PET`; this is the same trust model as separate RCH/EVT.
- **The one genuinely open Level-2 choice:** strict **pond-first** vs.
  **area-weighted** pond/soil evaporation (§7.4).

#### 7.4 The open choice: pond-first vs. area-weighted

Strict pond-first (current) implicitly assumes that, whenever a pond exists, it
covers enough of the cell to meet the demand at the potential rate *and* that
the residual still evaporates from bare soil on the same footprint — a mild
inconsistency, since a patch is either ponded or bare. The principled
generalization is an **area weighting** by a ponded-area fraction `f_pond(d)`:

```
E = f_pond · E_pot_openwater  +  (1 − f_pond) · krel · E_pot_soil
```

where `f_pond` is a smooth, C1 ramp of pond depth `d` toward 1 that → 0 as the
pond vanishes (e.g. `f_pond = smooth(d / d_scale)`). This removes the
"pond-covers-soil-but-soil-still-evaporates" artifact, preserves all invariants
(C1 in head, exact mass balance, total `≤ E_pot`), and **reduces to the current
pond-first behavior in the thin-but-present-pond limit** (`f_pond → 1`). It is
the Shuttleworth–Wallace / FAO-56 area idea reduced to a single tunable number.

**Recommendation.** Scope this as an application-driven decision rather than a
physics gap:

- For thin, transient ponding (most infiltration/runoff problems), **keep
  pond-first** and document it explicitly as the `f_pond → 1` limit of §7.4.
- Add `f_pond(d)` area weighting only if **partial-area ponding** matters
  (playa, intermittent wetland, long-residence ponds), where simultaneous
  open-water and bare-soil evaporation over one cell is physically real.

Either way, hold the invariants fixed: total actual evaporation `≤ E_pot`, C1 in
head, exact pond-balance closure, and the soil limiter is the same `krel` that
drives the flow kernel (`surface_exchange` in
`src/Model/ModelUtilities/VsfSurfaceExchange.f90`).

#### 7.5 Verification targets

Compare a 1-D, fine-vertical column (VSSB + VSET) against a code that resolves
these processes rigorously — SWAP (Ritchie partition + Black / Boesten–
Stroosnijder soil limiters) or HYDRUS (head-based limiter) — for:

1. a drying bare-soil column — does `krel`-limiting reproduce the stage-1 →
   stage-2 (energy- to supply-limited) falling rate? (Analytical/lysimeter
   two-stage benchmarks: Ritchie 1972; Idso et al. 1974.)
2. a ponded-then-evaporating case — the pond → soil handoff and water-balance
   closure (§3.4);
3. a shallow-water-table ET case — compare the extinction-depth behavior against
   Shah et al. (2007).

---

## VSSP (Seepage) Package — Status Handoff

_Last updated: 2026-10-20_

### Where we are
The **VSSP seepage boundary package** has been reworked from the old prototype
into a Richards-only seepage face that approximates conductance from NPF
conductivity × relative permeability. It **builds cleanly** and its **autotest
passes** (2 scenarios).

### Current VSSP behavior
- **Applies only to Richards (VSF) cells.** Requires the VSF package active;
  validated per-cell via VSF's `IVARSAT` array.
- **Input** (`stress_period_data`): `cellid  ihc  cl1  hwva  angldegx`
  - `ihc` — 0 = vertical face, 1 = horizontal face
  - `cl1` — distance from cell center to the face
  - `hwva` — seepage face area
  - `angldegx` — **reserved/unused** for now (kept in the definition)
- **Conductance:** `cond = krel(node) · K · hwva / cl1`, LAK-style:
  - vertical face (`ihc=0`) → `K33` (falls back to `K11` if `IK33=0`)
  - horizontal face (`ihc=1`) → `K11`
  - `krel` = NPF relative permeability (`KREL`), hardwired via `mem_setptr`
- **Seepage switch:** discharges (Robin term, `hcof=-cond`, `rhs=-cond·z`) while
  `head > z` (cell-center elevation); no-flow otherwise.
- The old penalty / ibound-toggle machinery and the ConductanceProvider
  abstraction were **removed**.

### Key files
- `src/Model/GroundWaterFlow/gwf-vssp.f90` — package: `effective_k`, `vssp_cf`,
  `vssp_ck` (Richards-cell check), pointer setup/teardown.
- `doc/mf6io/mf6ivar/dfn/gwf-vssp.dfn` — input definition.
- `src/Idm/gwf-vsspidm.f90` — generated IDM (via
  `utils/idmloader/scripts/dfn2f90.py`).
- `src/Model/GroundWaterFlow/gwf.f90` — model wiring.
- `autotest/test_gwf_vssp.py` — representative VSF+VSSP test (active / no-flow),
  asserts an analytic 2-node Robin solution.

### Also removed this session (by request)
- **SGD package** deleted entirely (source, IDM, dfn, test, build/registration).
- **ConductanceProvider** + **NpfConductanceProvider** modules deleted; NPF's
  now-unused `calc_eff_hy` removed too.

### How to build & test
```bash
# build
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
ninja -C builddir

# regenerate IDM after any dfn change
cd utils/idmloader/scripts && python dfn2f90.py --outdir ../../../src/Idm

# refresh flopy classes after any dfn change (for the test)
cd autotest && python update_flopy.py

# run the test (binary staged at bin/mf6)
cd autotest && python -m pytest test_gwf_vssp.py -q
```
Note: `bin/mf6` is a gitignored staged copy of `builddir/src/mf6`.

### Open items / next steps
1. **`DKRDH` Jacobian term** — `krel` is currently lagged; adding `dkr/dh`
   (NPF stores `DKRDH`) would make the Newton term consistent. _Pending._
2. **`angldegx` decision** — still reserved; decide whether to use it
   (horizontal anisotropy / full tensor) or drop it.
3. **Convergence near the switch** — the `head > z` on/off combined with
   head-dependent `krel` can chatter; watch in richer models.
4. **Partial-domain Richards** — VSF today requires all cells `ivarsat=1`, so the
   per-cell VSSP check mainly guards "VSF absent" until partial-domain support
   exists.

### Known caveats (physics/numerics)
- Grid-aligned `K11`/`K33` only — ignores rotated anisotropy (`ANGLE1/2/3`) and
  XT3D.
- Fixed user `hwva`/`cl1` — wetted area doesn't change as heads move.
- Uses the seepage cell's own `krel` (no upstream weighting).

---

## Uncapped river leakage for Richards flow

### Physical interpretation

A Richards cell can have negative water pressure without losing its ability
to receive water. Water moves according to total hydraulic head:

$$
h = z + \psi,
$$

where $z$ is elevation and $\psi$ is pressure head, which can be negative in
unsaturated soil.

For a wetted riverbed, the river provides head $H_r = STAGE$. Representing the
bed as a hydraulic resistance gives:

$$
Q = C(H_r - h),
$$

where $Q$ is positive into groundwater and $C$ is the effective
riverbed/interface conductance. Negative soil pressure increases the head
difference and draws water from the river. This represents capillary suction,
not an invalid head or an inactive cell.

### Why conventional RIV caps leakage

The conventional RIV relationship is:

$$
Q = C[H_r - \max(h, RBOT)].
$$

It prevents the downstream head seen by the bed from falling below `RBOT`.
This corresponds to assuming atmospheric pressure at the bed underside once
groundwater head falls below that elevation.

Richards flow does not require that assumption: water below the bed may remain
under suction rather than at atmospheric pressure. Removing the cap allows
that suction to influence exchange.

### Assumptions and limitations

Being a Richards cell does not, by itself, make constant-conductance uncapped
leakage physically valid. The approximation additionally assumes:

- A wetted river supplies water continuously.
- The bed maintains liquid continuity and approximately constant conductance.
- The cell head adequately represents the downstream interface head, or
  `COND` includes the unresolved resistance between them.

If the bed desaturates, its conductivity decreases. If an atmospheric-pressure
separation develops, suction is no longer transmitted through the interface
as assumed. Uncapped constant-conductance leakage can then overestimate
infiltration.

Richards conductivity inside the aquifer does not automatically correct a
missing or misrepresented interface resistance.

The defensible description is therefore:

> A connected, riverbed-resistance boundary that permits suction-driven
> infiltration into Richards soil.

It is not a universally correct river boundary for Richards flow.

### Minimal, formulation-independent design

Expose the choice as a boundary leakage law, rather than making RIV detect
whether VSF is active:

| Leakage law | Relationship | Interpretation |
|---|---|---|
| Capped (existing default) | $Q = C[H_r - \max(h, RBOT)]$ | Conventional disconnected-river approximation |
| Uncapped (opt-in) | $Q = C(H_r - h)$ | Connected, resistive riverbed with suction-driven infiltration |

Both laws use total hydraulic head. Neither requires RIV to inspect VSF soil
models, relative permeability, or cell formulation flags. The uncapped law
can also be meaningful for saturated flow.

For a first implementation, preserve the existing default and introduce an
explicit package-level leakage-law option. In uncapped mode, retain `RBOT`
for input compatibility and applicable riverbed geometry checks, but document
that it does not cap leakage.

Do not automatically multiply riverbed conductance by the receiving cell's
relative permeability: that would conflate riverbed resistance with aquifer
resistance. Explicit treatment of bed desaturation or aquifer half-cell
resistance would require a more detailed interface model.

This is a conceptual proposal only; no implementation is included.

---

## Feddes transpiration for MODFLOW 6 VSF

This is a parallel proposal to
[VS2D-compatible transpiration](#vs2d-compatible-transpiration-for-modflow-6-vsf), covering an
**uncompensated Feddes** root-water uptake model instead of the VS2D law. It is
the lower-effort, more MODFLOW-idiomatic option. The trade-off is that it is a
*different* model and will **not** reproduce the VS2D uptake law.

### 1. Physics: piecewise-linear pressure-head stress

Uptake from cell $i$ is local — a prescribed stress function times a normalized
root fraction times atmospheric demand:

$$
U_i = \alpha(\psi_i)\,\beta_i\,T_p\,A
$$

where:

- $\psi_i = h_i - z_i$: soil pressure head;
- $\alpha(\psi)\in[0,1]$: the Feddes stress response, piecewise-linear between
  four threshold heads;
- $\beta_i$: normalized root fraction in cell $i$, with $\sum_j \beta_j = 1$ over
  the rooted column;
- $T_p$: potential transpiration rate;
- $A$: column surface area.

The stress function uses four pressure-head thresholds:

$$
\alpha(\psi)=
\begin{cases}
0, & \psi \ge h_1 \ \text{(anoxia / ponded)}\\
\dfrac{\psi-h_1}{h_2-h_1}, & h_2 \le \psi < h_1\\
1, & h_3 \le \psi < h_2 \ \text{(optimal)}\\
\dfrac{\psi-h_4}{h_3-h_4}, & h_4 < \psi < h_3\\
0, & \psi \le h_4 \ \text{(wilting point)}
\end{cases}
$$

(Head convention: $h_1 > h_2 > h_3 > h_4$. The groundwater budget receives
$-U_i$.)

**Important consequences:**

- $\beta_i$ is a **normalized root fraction**, not an absolute conductance
  parameter — the opposite convention from VS2D's root activity $r(z)$.
- Uptake is purely local: it does **not** depend on soil conductivity
  $K(\psi)$, and there is no soil-to-root pressure-gradient transfer.
- There is **no column-sum cap and no inter-cell coupling**. Each cell's uptake
  depends only on its own pressure head.
- The oxygen-stress limb ($h_1,h_2$) is an explicit feature here, whereas VS2D
  has none.
- This is *uncompensated*. Compensation (shifting uptake to wetter layers)
  would reintroduce a column-sum rescaling and the inter-cell coupling that
  makes the VS2D law harder — so the simplicity holds only without it.

### 2. Architecture: a generalized EVT-style stress package

The structure is nearly identical to the VS2D proposal (a new boundary/stress
package on the [BndExtType pattern](src/Model/GroundWaterFlow/gwf-evt.f90),
requiring VSF), but the uptake kernel is far simpler and the coupling to VSF is
lighter.

Key differences from the VS2D design:

- **Reuse EVT as the template.** [gwf-evt.f90](src/Model/GroundWaterFlow/gwf-evt.f90)
  already implements piecewise-linear rate-vs-head with clean HCOF/RHS terms and
  a segment mechanism (`nseg`, `pxdp`, `petm`). Feddes is essentially "EVT keyed
  on pressure head instead of water-table depth," so most of the assembly logic
  has a direct analogue.
- **Minimal VSF coupling.** The package needs only the current pressure head
  (i.e. $h_i$ and cell-center elevation $z_i$); it does **not** call the VSF soil
  model for $K(\psi)$ or its derivative.
- **Diagonal Jacobian.** $\partial U_i/\partial\psi_j = 0$ for $i\ne j$, and
  $d\alpha/d\psi$ is a constant slope within each segment, so Newton terms are
  trivial and there is no Picard-vs-coupled-derivative decision.
- Evaluate pressure head from the **current nonlinear iterate**, as in the VS2D
  design — the published VSF pressure-head array
  ([vsf_cq](src/Model/GroundWaterFlow/gwf-vsf.f90)) is not guaranteed current
  during assembly.

Everything else carries over unchanged from the VS2D proposal: explicit rooted
columns (DISV/DISU-friendly), time series, bound names, cell-by-cell budgets,
coexistence rules with EVT/UZF to avoid double-counting, and the
[IDM workflow](IDM.md) for definitions, generated input code, docs, and FloPy.

### 2a. Distribution vs. reduction: keep column logic out of the solver

A recurring design question is whether the user should specify only the **top
cell** of a column and let MODFLOW walk the vertical stack, distribute demand by
a root-density-vs-depth profile, and (as a side benefit) make the package usable
with standard (non-Richards) NPF flow. The conclusion is to **keep the per-cell
list as the package contract** and treat column discovery + root profiling as a
*preprocessing* concern. The reasoning:

**The two concerns are separable.** The feature bundles two independent things:

1. **Distribution** — which cells receive how much atmospheric demand (column
   discovery + root-density-vs-depth profile).
2. **Reduction** — the local uptake function turning potential into actual
   uptake given the cell state.

The current design already isolates these: distribution is the per-cell `rate`
in the period list, and reduction is the pluggable
[EtUptakeModel](src/Model/ModelUtilities/EtUptakeModel.f90). "Top cell → walk
stack → apply profile" is purely a *distribution* feature; it touches neither
the reduction side nor the solver assembly.

**What moving distribution into MODFLOW costs.** Vertical-connectivity traversal
is trivial on DIS but needs vertical-connection logic on DISV and top/bot overlap
or `ihc==0` handling on DISU — the cell-stack machinery UZF already carries, with
its edge cases (idomain holes, pinch-outs, inactive cells). Root-profile
integration over each cell's vertical extent, partial cells at the root-zone
bottom, and time-varying rooting depth add VS2D's `VSRDF`-style logic. This
roughly doubles the package surface area and reintroduces the grid-topology edge
cases a flat list avoids (a list works on any discretization).

**The non-Richards angle largely overlaps EVT.** In standard NPF there is no
pressure head, so the reduction would key off depth-to-water-table + extinction
depth — essentially what [EVT](src/Model/GroundWaterFlow/gwf-evt.f90) already
does (including `fixed_cell`/highest-active handling). The genuinely new thing
VSE would add for non-Richards flow is the column distribution + root profile,
not the sink assembly (HCOF/RHS at a node is identical regardless of the state
variable). Supporting it is a clean generalization of the uptake abstraction:
`fc_terms(h, z, c)` currently assumes pressure head, so a depth/saturation-based
model would need its state input generalized — a small change, not a redesign.

**Reduction already tracks the solution; only the distribution is static.** It is
worth separating two solution-dependent quantities. The per-cell *reduction*
$\alpha(\psi_i)$ with $\psi_i = h_i - z_i$ is evaluated from the current
nonlinear iterate every assembly, so it is already fully dynamic — a cell that
floods ($\psi \ge h_1$, anoxia) or dries ($\psi \le h_4$, wilting) ramps its
uptake to zero on its own, and a cell near a moving water table tracks it
smoothly. What a preprocessed (or period-specified) list fixes is the
*distribution*: the root-zone cell set and the potential weights
$\beta_i$. Those are **plant properties, not water-state properties**, so holding
them static within a stress period is physically complete; seasonal change
(growing roots, changing density or demand) is naturally expressed by new values
in successive stress-period blocks.

**Where in-model dynamic logic genuinely earns its keep.** Given the above, for
the **uncompensated Richards (Feddes) model there is no within-period driver** for
a dynamic cell list or weights — every water-state dependence is in
$\alpha(\psi)$, and the potential distribution is a static plant property. Only
two features would reintroduce solution-dependent distribution, and neither is in
scope:

1. **Column compensation** — reallocating the *unmet* demand of stressed or
   flooded cells to less-stressed cells in the same column (VS2D-style capping,
   $\sum$ over the column $\le T_p A$, or compensated Feddes). This depends on
   the current $\psi$ of every cell in the column, so it is inherently a runtime
   column operation. It is the *only* genuine within-period redistribution
   driver, and it was deliberately deferred by choosing uncompensated Feddes.
2. **EVT-style relocation** — moving the sink to the cell that currently contains
   the water table. This is a **non-Richards concern**: with standard NPF the ET
   cell changes as the table drops ([EVT](src/Model/GroundWaterFlow/gwf-evt.f90)'s
   `fixed_cell`/highest-active logic). In the Richards case cells are never
   relocated — each root-zone cell is always present and governed by
   $\alpha(\psi)$ — so this does not apply.

**Recommendation.** Keep the per-cell list as the package contract. Deliver "top
cell + root profile" as a FloPy-side helper / documented recipe that emits the
list (inputs: top cellid, rooting depth, profile → per-cell `rate`), with no
Fortran change, working on all grids immediately. Keep the uptake function
pluggable; when a standard-NPF use case lands, add a depth/saturation-based
uptake model and minimally generalize the model input. Promote column discovery
into MODFLOW only when committing to column compensation (or, for non-Richards
flow, EVT-style relocation) — at which point it is a deliberate, UZF-like
feature with its own dimensions/period input, not a quiet change to this
package.

### 3. Implementation: small, verifiable increments

1. Stress function $\alpha(\psi)$, thresholds, and normalized root profile.
2. Package wiring, nonlinear assembly (diagonal HCOF/RHS), budgets, observations.
3. Time-varying demand and vegetation parameters.
4. Documentation and input-generation integration.

Because the kernel is local, phase 2 is substantially lighter than in the VS2D
plan (no column normalization, no coupled derivatives). This package can also
serve as a **phase-0 scaffold** — wiring, budgets, and observations shared with
a later VS2D kernel if both are eventually wanted.

Soil evaporation is a separate surface-boundary problem, handled by the VSSB
Package rather than by this transpiration package.

### 4. Verification: agree on these before coding

| Case | Required check |
|---|---|
| Prescribed pressure heads | Exact cell uptake $\alpha(\psi_i)\beta_i T_p A$ |
| Each stress limb | Correct $\alpha$ on anoxia, optimal, and dryness segments |
| At thresholds $h_1..h_4$ | Zero/unit/zero behavior exactly at breakpoints |
| Dry soil / zero demand / no roots | Zero extraction |
| Root fraction normalization | Column uptake = $\alpha$-weighted $T_p A$; handled when some cells stressed |
| Transient closed column | Storage loss equals integrated uptake |
| Solver variants | Picard/Newton agree within declared tolerances |
| Invalid inputs | Non-monotonic thresholds, unnormalized roots rejected clearly |

Note there is **no VS2D reference-column case** here — by construction this model
does not target VS2D fidelity.

### Choosing between the two

| | Uncompensated Feddes (this doc) | VS2D law ([other doc](#vs2d-compatible-transpiration-for-modflow-6-vsf)) |
|---|---|---|
| Implementation effort | Lower (≈ generalized EVT) | Higher |
| Jacobian | Diagonal, local | Column-coupled when cap active |
| Touches VSF soil model | Pressure head only | Needs $K(\psi)$ and its derivative |
| Inter-cell redistribution | Only via compensation (adds coupling) | Built in via conductivity + column cap |
| Reproduces the VS2D uptake law | **No** | Yes |
| Natural MODFLOW analogue | EVT (pressure-head-keyed) | New nonlinear sink |

**Decision pending:** if strict VS2D compatibility is a hard
requirement, use the VS2D law despite the extra complexity. If a defensible
transpiration model in VSF is sufficient, uncompensated Feddes is the simpler
path and could scaffold a later VS2D kernel.

---

## Alternatives (not yet implemented)

These documents describe alternative or superseded approaches for the VSF family that have not been implemented. They are retained here for reference and context.

### SPG Seepage Package — Notes and Findings

#### Overview

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

#### Budget sign convention

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

#### Testing

- `autotest/test_gwf_spg.py`: a confined two-cell column (CHD top, SPG bottom)
  with four cases — `{penalty, toggle} x {discharging, no-flow}`. Checks the
  held head equals the seepage elevation, the discharge balances the
  constant-head inflow, and the no-flow state carries zero flow.
- `autotest/test_gwf_vsf_wt.py`: the `wt-vsf-spg` case adds SPG seepage cells on
  the `+x` boundary of a transient 2D Richards/VSF model and verifies per-cell
  mass balance (FLOW-JA-FACE residuals).

#### Findings from the VSF integration

1. **Toggle conflicts with VSF.** VSF claims the entire grid for its own flow and
   storage formulation (`iformulation` is set for all cells). When SPG toggle
   switches a cell to a true constant head (`ibound < 0`), VSF still tries to
   apply its formulation to that cell. This produces repeated
   `VSF convertible cell error` messages and a small residual imbalance. The
   **penalty** method keeps the cell active (`ibound > 0`), so it composes
   cleanly with VSF, and is the method used in the `wt-vsf-spg` test.

2. **Seepage on/off is a discontinuous nonlinearity.** Switching a cell between
   the held and no-flow states is discontinuous: a tiny head change can flip the
   state and cause a jump in the boundary flow. As a result the model converges
   to a slightly looser per-cell FLOW-JA-FACE residual near the seepage cells
   (~`3.7e-6`) than the rest of the grid. Both the penalty and toggle methods
   converge to the same residual level, confirming this comes from the seepage
   switch itself rather than from the penalty conductance. The test therefore
   uses `res_atol = 1e-5` for the seepage case and keeps the strict `1e-6`
   tolerance for the non-seepage cases.

#### Future work

- **Continue the investigation with the toggle method and how to avoid the
  reported errors.** In particular, determine how to reconcile the toggle
  (true constant-head) method with packages such as VSF that claim the whole
  grid, so that toggling a cell to constant head does not trigger
  `VSF convertible cell error` or introduce a residual imbalance. Candidate
  directions to explore:
  - Have SPG coordinate with VSF (and similar whole-grid packages) to release or
    skip the toggled cell's formulation while it is held as constant head.
  - Investigate whether the seepage on/off switch can be smoothed (continuous
    transition) to improve convergence and tighten the residual.
  - Evaluate Newton-Raphson handling of the state switch under VSF.

### VS2D-compatible transpiration for MODFLOW 6 VSF

> **Status: planned later addition.** This was originally proposed as the
> first implementation. That decision was revised: **uncompensated Feddes is the
> day-one model** (see [Feddes transpiration proposal](VSF_DEVELOPMENT.md#feddes-transpiration-for-modflow-6-vsf)),
> because it is substantially simpler to implement and carries higher immediate
> demand. The VS2D law remains planned — it plugs into the same VSE package and
> uptake-model abstraction as an additional model — and this document is the
> reference for that work. It is kept design-complete but is **not** the next
> thing to build.

This document reproduces the VS2D root-water uptake law for eventual
VS2D compatibility. The four phases below still
fit; verification cases should be defined alongside the physics, rather than
after implementation.

#### 1. Physics: reproduce the VS2D uptake law

The [original VS2D manual](https://pubs.usgs.gov/publication/wri834099),
particularly equations 50-52, was checked against the released
[VS2DH 3.3 source](https://water.usgs.gov/water-resources/software/VS2DI/).

Using positive quantities for extraction, the uncapped uptake from cell $i$ is:

$$
U_i^* = V_i K_i(\psi_i)\,r_i(t)\,
        \max\!\left(\psi_i-\psi_{\mathrm{root}}(t),0\right)
$$

where:

- $\psi_i=h_i-z_i$: soil pressure head;
- $K_i(\psi_i)$: soil hydraulic conductivity;
- $r_i(t)$: root activity, with dimensions $L^{-2}$;
- $\psi_{\mathrm{root}}$: prescribed root pressure head;
- $V_i$: cell volume, accounting for partial root-zone penetration.

For each rooted column, potential transpiration limits total extraction:

$$
U_i = \lambda U_i^*, \qquad
\lambda =
\begin{cases}
1, & \sum_j U_j^*=0,\\
\min\!\left(1,\dfrac{A T_p}{\sum_j U_j^*}\right), & \text{otherwise}.
\end{cases}
$$

The groundwater budget receives $-U_i$.

**Important consequences:**

- Root activity is an absolute conductance-related parameter, **not a normalized
  root-density fraction**.
- Uptake stops when soil pressure head reaches or falls below root pressure head.
- Wetter, more conductive cells can supply more of the column demand. This
  provides redistribution without a separate empirical compensation factor.
- There is **no explicit Feddes-style oxygen-stress cutoff**. Adding one would
  change the reference physics.
- VS2D uses a depth-linear root-activity profile, with partial-cell treatment at
  the root-zone bottom.

For the first version, retain the reference inputs: potential transpiration,
root pressure head, rooting depth, and root activities at the top and bottom.

**What drives this model:** general **VS2D compatibility** is the sole driver.
The originally-cited Naranjo (2023) study is a canal-seepage / heat-tracer
application of VS2D that does **not** use transpiration, so it imposes no
constraint on this uptake model. The reference target is therefore VS2D's
documented uptake law itself ([manual](https://pubs.usgs.gov/publication/wri834099)
and [VS2DH 3.3 source](https://water.usgs.gov/water-resources/software/VS2DI/)),
not a specific published ET application.

#### 2. Architecture: fits the existing VSE package

The VS2D law is an additional uptake model inside the **existing VSE package**
(the [BndExtType pattern](src/Model/GroundWaterFlow/gwf-evt.f90), requiring VSF),
not a new package. The uptake-model abstraction
([EtUptakeModel](src/Model/ModelUtilities/EtUptakeModel.f90) + factory) already
built for Feddes is the extension point: VS2D becomes a new model registered in
the factory, alongside the existing FEDDES implementation.

That keeps responsibilities clear:

- **VSF:** retention, conductivity, flow, and storage.
- **VSE package:** vegetation parameters, atmospheric demand, uptake,
  observations, and budgets.

Specific design points:

- Evaluate pressure head from the **current nonlinear iterate**. The published
  VSF pressure-head array is updated in
  [vsf_cq](src/Model/GroundWaterFlow/gwf-vsf.f90), so it should not be assumed
  current during assembly.
- Reuse VSF's hydraulic-property definitions; do not introduce a second
  retention/conductivity implementation.
- Represent rooted columns explicitly, avoiding assumptions that would prevent
  later DISV/DISU support.
- Support standard time series, bound names, and cell-by-cell budgets.
- Define coexistence rules with EVT/UZF to prevent unintended double-counting.
- Follow the [IDM workflow](IDM.md) for definitions, generated input code,
  documentation, and FloPy integration.

**The main numerical design issue is the column cap:** normalization couples
uptake between cells. Its full Newton derivative is not merely diagonal. We
should deliberately choose and test a lagged/Picard treatment versus consistent
coupled derivatives.

> **Note — distribution vs. reduction.** Whether column discovery and
> root-depth distribution belong in the solver or in preprocessing is discussed
> in the Feddes proposal, [§2a Distribution vs. reduction](VSF_DEVELOPMENT.md#2a-distribution-vs-reduction-keep-column-logic-out-of-the-solver).
> That analysis is shared by both models; the column cap above is exactly the
> *column compensation* case identified there as one of the two reasons to pull
> stack-walking into MODFLOW.

#### 3. Implementation: small, verifiable increments

1. Uptake calculation and reference parameters.
2. Register the model in the uptake-model factory; nonlinear assembly with the
   column cap, budgets, and observations (the VSE package already provides the
   wiring).
3. Time-varying demand and vegetation parameters.
4. Documentation and input-generation integration.

Soil evaporation is a separate surface-boundary problem, handled by the VSSB
Package rather than by this transpiration package.

#### 4. Verification: agree on these before coding

| Case | Required check |
|---|---|
| Prescribed pressure heads | Exact cell uptake and column-cap scaling |
| Dry soil / zero demand / no roots | Zero extraction |
| Partial root-zone cell | Correct integrated root activity |
| Mixed wet and dry layers | Correct redistribution and demand limit |
| Transient closed column | Storage loss equals integrated uptake |
| Solver variants | Picard/Newton agree within declared tolerances |
| VS2D reference column | Pressure profiles and cumulative transpiration converge toward the reference |
| Invalid inputs | Explicit, informative rejection |

**Next milestone (when VS2D work is scheduled):** settle the VS2D equations,
conductivity convention, root-zone geometry, and paper-specific parameters
before extending the uptake-model interface. Note this likely requires
generalizing `fc_terms` beyond the current pressure-head-only inputs (it needs
$K(\psi)$ and its derivative) and introducing the column cap — the
*compensation* coupling discussed in §2's note.
