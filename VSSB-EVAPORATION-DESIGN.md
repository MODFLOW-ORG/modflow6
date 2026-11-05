# Design note: soil and pond evaporation in VSSB

Status: implemented. The evaporation design in Sections 2-4 is implemented in the
VSSB Package and its surface-exchange kernel. The seepage-face slope extension in
Section 6 remains out of scope.

This note records the design for adding evaporation (from bare soil and from
ponded water) to the Variably Saturated Surface Boundary (VSSB) Package, the
relationship between VSSB and the seepage-face package (VSSP), and the
complexity of a possible future extension to seepage faces under a slope that is
not aligned with the grid. The slope extension is discussed for scoping only and
is explicitly out of scope for implementation.

## 1. Current architecture

The variably-saturated family is a core Richards flow package plus boundary
packages:

- `src/Model/GroundWaterFlow/gwf-vsf.f90` — the Richards (unsaturated) flow
  package. It flags Richards cells (`IUNSAT`) and publishes relative
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

The current documentation already names the gap to be filled: *"Evaporation from
ponded water is not simulated by the VSSB Package or by any other package"*
(`doc/mf6io/gwf/vssb.tex`).

## 2. Evaporation demand and wiring (decision: Option A)

### 2.1 Two separate partitioning levels

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

### 2.2 How VSET defines its demand today

VSET does not read a shared evapotranspiration demand. It reads a per-cell
`rate` in its own PERIOD block (`gwf-vset.dfn`), time-series-enabled, documented
as the potential transpiration flux rate for the cell. That value flows into
`vset_cf` as `this%rate(i)`, is multiplied by cell area, and is reduced to
actual uptake by the Feddes response. There is no "ET demand bus" to tap into.

### 2.3 Decision: Option A — VSSB owns `E_pot`

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
VSSB already uses to read NPF `KREL`/`DKRDH` and VSF `IUNSAT`), without changing
VSSB's input surface. Choosing Option A now does not foreclose it.

### 2.4 Decisions

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

## 3. Evaporation formulation in VSSB

### 3.1 The two sinks

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

### 3.2 Demand partitioning (Level 2)

Potential evaporation `E_pot` is one demand. Open water evaporates first at the
potential rate; only the residual `max(E_pot - E_pond, 0)` is exposed to the
soil, where it is further reduced by soil dryness. This keeps total evaporation
`<= E_pot` and avoids double counting. Because both the pond term and the
residual are capacity-limited, they fit the kernel's existing smoothed-minimum
machinery. The cleanest implementation extends `surface_exchange` to carry the
evaporation demand and emit the extra sink terms, rather than bolting
evaporation on outside the kernel (which would break the exact pond-balance
closure and the Newton consistency).

### 3.3 Linearization / convergence

The soil-evaporation sink needs `hcof`/`rhs` contributions and a Newton term for
its dependence on head through soil wetness, mirroring `vssb_fn` and the `DKRDH`
handling already in the kernel. Pond evaporation enters the pond update and the
runoff residual; its head dependence is weak, but the smoothed limiters must
stay C1 to preserve the current convergence behavior.

### 3.4 Interaction with exfiltration

Exfiltrated groundwater lands in the pond, so pond evaporation can consume
exfiltrated water. This is physically correct and the budget still closes, but
it means the GWF and EVAPORATION terms interact; the test suite should assert
the closed water balance explicitly.

### 3.5 Derived kernel equations (as implemented)

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

## 4. Code touch points

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

## 5. Relationship to VSSP (exfiltration)

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

## 6. Future extension (out of scope): seepage faces under a slope

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

### A note on the 3D geometry

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
