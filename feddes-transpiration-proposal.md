# Feddes transpiration for MODFLOW 6 VSF

This is a parallel proposal to
[vs2d-transpiration-proposal.md](vs2d-transpiration-proposal.md), covering an
**uncompensated Feddes** root-water uptake model instead of the VS2D law. It is
the lower-effort, more MODFLOW-idiomatic option. The trade-off is that it is a
*different* model and will **not** reproduce the VS2D uptake law.

## 1. Physics: piecewise-linear pressure-head stress

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

## 2. Architecture: a generalized EVT-style stress package

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

## 2a. Distribution vs. reduction: keep column logic out of the solver

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

## 3. Implementation: small, verifiable increments

1. Stress function $\alpha(\psi)$, thresholds, and normalized root profile.
2. Package wiring, nonlinear assembly (diagonal HCOF/RHS), budgets, observations.
3. Time-varying demand and vegetation parameters.
4. Documentation and input-generation integration.

Because the kernel is local, phase 2 is substantially lighter than in the VS2D
plan (no column normalization, no coupled derivatives). This package can also
serve as a **phase-0 scaffold** — wiring, budgets, and observations shared with
a later VS2D kernel if both are eventually wanted.

Soil evaporation remains a later, separate surface-boundary problem.

## 4. Verification: agree on these before coding

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

## Choosing between the two

| | Uncompensated Feddes (this doc) | VS2D law ([other doc](vs2d-transpiration-proposal.md)) |
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
