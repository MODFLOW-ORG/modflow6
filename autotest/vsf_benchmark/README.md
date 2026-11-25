# VSF performance studies

How to make variably saturated flow (VSF) in MODFLOW 6 fast at the scale of
regional models (1 million cells and more), and what each step gains. This
folder holds the two scripts used and the results of two studies, 2026-11-23/24:

- [`vsf_big.py`](vsf_big.py): one regional benchmark model, used for
  VSF_1D_ONLY, the preconditioner and backtracking.
- [`vsf_study.py`](vsf_study.py): a parametrised study over soils, scenarios and
  grid sizes, used for everything else.

A shared, editable copy of the second study lives in Claude Docs:
<https://claude.ai/artifact/KkcCmcZR74oNy3trQCY34H>.

## Contents

- [Summary](#summary)
- [Recommended settings](#recommended-settings)
- [Background: where the time goes](#background-where-the-time-goes)
- [Test models](#test-models)
  - [Regional benchmark (vsf_big.py)](#regional-benchmark-vsf_bigpy)
  - [Parametrised study (vsf_study.py)](#parametrised-study-vsf_studypy)
- [Methods and results](#methods-and-results)
  - [Fused soil-function evaluations](#fused-soil-function-evaluations)
  - [Preconditioner](#preconditioner)
  - [Backtracking](#backtracking)
  - [Looser inner tolerance](#looser-inner-tolerance)
  - [VSF_1D_ONLY and condensing the 1D columns](#vsf_1d_only-and-condensing-the-1d-columns)
  - [Column ordering](#column-ordering)
  - [Limiting the Newton Jacobian term](#limiting-the-newton-jacobian-term)
  - [Under-relaxation](#under-relaxation)
  - [Saturation step limiter](#saturation-step-limiter)
  - [Smoothing kr near saturation](#smoothing-kr-near-saturation)
- [Final comparison at 1M cells](#final-comparison-at-1m-cells)
- [Open items](#open-items)
- [Reproducing the results](#reproducing-the-results)

## Summary

VSF is now about as fast as VS2DT on the comparison problem, and a
million-cell model runs 3–8× faster with MILU0 than with the COMPLEX defaults.
On top of that, a per-cell step limiter, smoothing of kr near saturation and
faster-recovering under-relaxation make the 1M-cell runs another 2.2× faster
(1.4–4.8× per case), and make cases converge that failed before.

| Method | Kind | Impact | Status |
| --- | --- | --- | --- |
| [Fused soil-function evaluations](#fused-soil-function-evaluations) | code | 4.5× faster on the VS2DT strip benchmark, on par with VS2DT | committed |
| [MILU0 preconditioner](#preconditioner) | input | 2.6–7.6× faster at 1M cells | recommended |
| [Backtracking off](#backtracking) | input | a third of the run saved; it never triggered | recommended |
| [INNER_DVCLOSE = OUTER_DVCLOSE](#looser-inner-tolerance) | input | 1.3–1.5× faster at 1M cells | recommended |
| [VSF_1D_ONLY](#vsf_1d_only-and-condensing-the-1d-columns) | code | at most 1.06–1.2×, can be slower | not merged (branch `vsf-1d-only`) |
| [Column ordering](#column-ordering) | input (DISU) | none | dropped |
| [Limiting the Newton Jacobian](#limiting-the-newton-jacobian-term) | code | slower in most cases | dropped |
| [DBD with fast recovery](#under-relaxation) | input | 20–45% fewer outer iterations; needs kr smoothing | recommended |
| [Saturation step limiter](#saturation-step-limiter) | code | 2.2× together with the rows above and below; sand converges | `NEWTON_DSAT_MAX`, default 0.2 |
| [Smoothing kr near saturation](#smoothing-kr-near-saturation) | code | fixes a 1915 s → 75 s failure; heads change ≤ 0.4 cm | `KR_SMOOTHING`, `HAEVGN`, off by default |

One steep-soil case (Haverkamp `deep`) still fails in every configuration, and
IMS crashes rather than retrying when BiCGSTAB overflows. Both are under
[open items](#open-items).

## Recommended settings

For VSF models with Newton. This is the "candidate" set-up of the
[final comparison](#final-comparison-at-1m-cells): 2.2× faster than MILU0 alone
at 1M cells. Without UNDER_RELAXATION_KAPPA 0.2 it is 1.8× faster.

```
# IMS
BEGIN NONLINEAR
  OUTER_DVCLOSE 0.01
  UNDER_RELAXATION DBD
  UNDER_RELAXATION_THETA 0.8
  UNDER_RELAXATION_KAPPA 0.2
  BACKTRACKING_NUMBER 0
END NONLINEAR
BEGIN LINEAR
  INNER_DVCLOSE 0.01
  LINEAR_ACCELERATION BICGSTAB
  PRECONDITIONER_LEVELS 0
  PRECONDITIONER_DROP_TOLERANCE 0.0
  RELAXATION_FACTOR 0.97
END LINEAR

# VSF (Van Genuchten soils)
BEGIN OPTIONS
  SOIL_MODEL VANGENUCHTEN
  KR_SMOOTHING 2.0        # cm; 0.02 in a model in metres
  # NEWTON_DSAT_MAX 0.2   # default, shown for reference
END OPTIONS
```

DVCLOSE values and `KR_SMOOTHING` are in model length units (cm here) and scale
with them.

## Background: where the time goes

VSF adds Richards-equation flow to the GWF model. Each unsaturated cell gets a
pressure-dependent relative permeability kr(ψ) and saturation S(ψ) from a soil
model (Van Genuchten or Haverkamp), so every outer iteration evaluates those
curves and their derivatives in every cell.

The trigger for this work was the VS2DT comparison in the richards-modflow6
repository (`examples/vs2d_comparison/05_convergence_and_performance.ipynb`). On
a small strip problem VS2DT ran in about 1.1 s, while VSF needed 5.2 s with
Newton and 5.4 s with Picard. Our production models are 1 million cells and up,
so both studies measure at that scale. A run's time splits into three parts that
respond to different fixes:

- **Formulate**: building the matrix, including the soil-function evaluations.
  40–60% of the run once the linear solve is fixed.
- **Linear solve**: the inner BiCGSTAB iterations with their preconditioner.
  Up to 70% with the COMPLEX default, 25–50% with MILU0.
- **Number of outer (Newton) iterations**: multiplies both. In hard cases this
  dominates everything else.

## Test models

### Regional benchmark (vsf_big.py)

One model, 183 × 183 × 30 = 1.0 million cells, in cm and s:

- **Layers.** 10 fine layers of 30 cm for the unsaturated zone, over 20 coarse
  aquifer layers of 5 m.
- **Conductivity.** A lognormal K field in the aquifer, an aquitard in coarse
  layers 7–8 (K = 1e-5 cm/s), and Celia/Haverkamp sand on top.
- **Stresses.** Patchy recharge of 20 mm/d for 2 days, then 1 mm/d for 3 days.
  Drains at 1.2 m depth; fixed heads hold the lateral boundary at 1.5 m depth.
- **Cell size** 100 m or 10 m (`--delr`). Newton, ATS with a maximum step of 6 h.

### Parametrised study (vsf_study.py)

A synthetic regional block with the same layout (10–20 soil layers of 30–50 cm
over 20 aquifer layers of 5 m, 10 m cells), varied over size, soil and scenario:

- n × n columns; n = 183 gives about 1 million cells, n = 58 about 100 000.
- Transient, 5 days in two stress periods, ATS from 30 s up to 6 h.

| Soil | Model | Ks (cm/s) | Character |
| --- | --- | --- | --- |
| Loam | Van Genuchten, n = 1.56 | 2.9e-4 | medium |
| Silt loam | Van Genuchten, n = 1.41 | 1.25e-4 | fine, gentle curves |
| Sand | Van Genuchten, n = 2.68 | 8.3e-3 | coarse, steep curves |
| Haverkamp sand | Haverkamp (Celia et al., 1990) | 9.4e-3 | steep kr, classic test soil |

| Scenario | What happens |
| --- | --- |
| `shallow` | water table at 1.5 m, 2 days of rain at 5% of Ks, then drainage to drains at 1.2 m |
| `rain20` | as shallow, with a fixed 20 mm/d rain then 1 mm/d (as in vsf_big.py) |
| `deep` | water table at 9 m; a wetting front travels down through dry soil |
| `layered` | as shallow, with alternating layers 10× less conductive |
| `wet` | rain above Ks on a shallow profile; surface drains take the excess |

What we report per run:

- **Wall time** of the MF6 run. 1M-cell runs were made one at a time. Some
  100k-cell comparisons ran 4–5 cases in parallel on a 16-core machine; there the
  iteration counts are the main metric.
- **Outer iterations** (including those of failed ATS attempts), **inner
  iterations**, and **failed attempts** (time steps retried with a smaller step).
- **Accuracy**: the largest head (cm) and saturation difference at the end of the
  run against the reference. Head differences below the outer closure of
  0.01 cm are noise.

Unless stated otherwise the solver is IMS COMPLEX with BICGSTAB, OUTER_DVCLOSE
0.01 cm, RCLOSE 10 STRICT, INNER_DVCLOSE 0.001, and Newton.

## Methods and results

### Fused soil-function evaluations

**Result: VSF Newton went from 4.5× slower than VS2DT to on par (5.15 s → 1.15 s
on the strip benchmark); Picard 5.39 s → 2.47 s. Committed.**

Profiling showed the run was dominated by `pow`, `exp` and `log` calls inside
the soil models. Three things caused that:

1. dkr/dψ for the Newton Jacobian was a central difference per connection: two
   extra kr evaluations, each of which called the saturation function again.
2. kr, S and the moisture capacity Cm were computed in separate calls that each
   recomputed the same powers.
3. The storage term re-evaluated S(ψ_old) at the start-of-step head in every
   outer iteration, although it only changes once per time step.

The changes:

- New type-bound procedures `kr_dkrdh` and `sat_capacity` on the soil model
  return kr with its analytical derivative, and S with Cm, from one shared set of
  powers. For Van Genuchten, with x = (α|ψ|)^n, both need only three
  exponentials and two logarithms. The Newton terms reuse the per-cell `dkrdh`.
- Close to saturation (ψ > −1e-6) the central difference is kept, because the
  exact dkr/dψ is singular there for n < 2.
- S(ψ_old) is cached per cell and refreshed only when a new time step starts.

The results are unchanged to solver precision. All later runs include this.

### Preconditioner

**Result: MILU0 is 2.6–7.6× faster at 1M cells than the COMPLEX default ILUT(5),
with no change in the answer. An input setting.**

On the 10 m regional benchmark (full VSF, Newton):

| Preconditioner | Wall time | Linear solve | Formulate | Outer / inner iterations | Linear solve per inner iteration | Max head difference vs ILUT(5) |
| --- | --- | --- | --- | --- | --- | --- |
| **MILU0 (relaxation 0.97)** | **177 s** | **49 s** | 62 s | 211 / 855 | 0.06 s | 0.003 cm |
| ILU0 | 212 s | 81 s | 63 s | 216 / 1508 | 0.05 s | 0.005 cm |
| ILUT(1, 1e-3) | 259 s | 109 s | 73 s | 252 / 2812 | 0.04 s | 0.04 cm |
| ILUT(5, 1e-4), the COMPLEX default | 455 s | 326 s | 62 s | 212 / 1090 | 0.30 s | |
| ILU0 with RCM reordering | failed after 8.6 h | 27,230 s | 704 s | 2294 outer | ≈12 s per solve | |

Across the parametrised study at 1M cells (MILU0 here also has backtracking off,
see the next section):

| Case (1M cells, Newton) | COMPLEX (s) | MILU0 (s) | Speed-up |
| --- | --- | --- | --- |
| Loam `shallow` | 478 | 126 | 3.8× |
| Loam `layered` | 1166, failed | 154 | — |
| Loam `deep` | 867 | 632 | 1.4× |
| Silt loam `shallow` | 349 | 90 | 3.9× |
| Silt loam `deep` | 552 | 212 | 2.6× |
| Silt loam `layered` | 537 | 99 | 5.4× |
| Haverkamp `rain20` | 445 | 115 | 3.9× |

- **ILUT(5) costs more than it saves.** Each inner iteration costs 5–6× more than
  with ILU0 or MILU0, and MILU0 even needs fewer inner iterations.
- **Why ILU0-type preconditioners suffice.** The soil layers are thin and their
  vertical conductances dominate, so the matrix is close to a set of tridiagonal
  columns. ILU0 in natural (layer-major) ordering already captures that
  coupling. With MILU0 the inner iterations per outer iteration stayed flat
  (about 2–4) from 100k to 1M cells.
- **RCM reordering breaks convergence** for this problem and should be avoided.
- **MILU0 is the linear setup of `COMPLEXITY MODERATE`.** To keep the COMPLEX
  nonlinear settings, set `PRECONDITIONER_LEVELS 0`,
  `PRECONDITIONER_DROP_TOLERANCE 0.0` and `RELAXATION_FACTOR 0.97` explicitly.
- Loam `deep` gains least because it is dominated by outer iterations (952),
  which the later methods address.

### Backtracking

**Result: backtracking never triggered, yet checking for it cost a third of the
MILU0 run. `BACKTRACKING_NUMBER 0` took the 10 m benchmark from 177 s to 125 s
with bit-identical heads.**

| 10 m benchmark, MILU0, Newton | Wall time | Solve | Formulate | Linear solve | Outer / inner iterations |
| --- | --- | --- | --- | --- | --- |
| COMPLEX backtracking (20 tries) | 177 s | 170 s | 62 s | 49 s | 211 / 855 |
| `BACKTRACKING_NUMBER 0` | 125 s | 118 s | 64 s | 53 s | 211 / 855 |

- **Settings.** COMPLEX sets `BACKTRACKING_NUMBER 20`, tolerance 1.05, reduction
  factor 0.1 and residual limit 0.002. SIMPLE and MODERATE switch it off.
- **How it works** (`sln_backtracking` in `NumericalSolution.f90`). Every outer
  iteration rebuilds the full system (Picard terms, `inewton=0`) only to compute
  the residual L2 norm. A step is cut back by 0.1, with a rebuild each time, only
  if the residual grew more than 5% and 0.1·max|Δh| ≥ `dvclose`.
- **Observed.** The backtrack count was 0 in all 211 outer iterations of the
  ILUT(5) and MILU0 runs.
- **Where it does help.** In the parametrised study backtracking is the only
  global setting that makes the sand case converge and it helps loam `deep`, but
  it makes Haverkamp `shallow` 3.5× slower (see
  [under-relaxation](#under-relaxation)). The
  [saturation step limiter](#saturation-step-limiter) now covers those cases
  locally, so backtracking can stay off.

### Looser inner tolerance

**Result: INNER_DVCLOSE equal to OUTER_DVCLOSE (0.01 cm) is 1.3–1.5× faster at
1M cells, with head differences up to 0.5 cm in the hardest case and below
0.05 cm elsewhere. An input setting.**

Each Newton step only needs to point in roughly the right direction, so solving
the linear system to 0.001 cm is more than necessary. Three looser variants with
MILU0 at 1M cells:

- `dv01`: INNER_DVCLOSE 0.01 (= outer closure).
- `rel1` / `rel2`: stop when the residual norm has dropped by a factor 10 or 100
  (RCLOSE 0.1 or 0.01 with RELATIVE_RCLOSE).

| Case (1M cells) | Reference (s) | `dv01` (s) | `rel1` (s) | `rel2` (s) | Largest head difference, `dv01` (cm) |
| --- | --- | --- | --- | --- | --- |
| Loam `shallow` | 126 | **87** | 87 | 86 | 0.002 |
| Loam `deep` | 632 | **443** | 537 | 448 | 0.52 |
| Silt loam `shallow` | 90 | **60** | 67 | 64 | 0.006 |
| Silt loam `deep` | 212 | **143** | 154 | 154 | <0.001 |
| Haverkamp `rain20` | 115 | **85** | 135 | 97 | 0.05 |

`dv01` also needed 20–25% fewer *outer* iterations, because of the STRICT
criterion: an outer iteration only counts as converged when the linear solve
converges in one inner iteration, which happens sooner with a looser inner
tolerance. The relative-residual variants bring the inner work down to about one
iteration per outer iteration, but `rel1` then needs more outer iterations.

An adaptive (Eisenstat–Walker) inner tolerance was considered but not built:
with only 1–4 inner iterations per outer iteration left, it could at best match
`dv01`.

### VSF_1D_ONLY and condensing the 1D columns

**Result: restricting Richards flow to vertical connections gives at most
1.06–1.2× and can be slower; condensing the columns would save 7–8% at best.
The option lives on branch `vsf-1d-only` and is not merged here.**

The option applies Richards flow only to vertical connections (ihc = 0). A
horizontal connection conducts flow only when both cells are fully saturated
(ψ ≥ 0 at the cell node), with the saturated conductance; partially saturated
cells get no horizontal flow. It is tested with Picard and Newton in
`autotest/test_gwf_vsf_1d.py` on that branch.

| Model | Scenario | Effect of VSF_1D_ONLY |
| --- | --- | --- |
| 16,000 cells | water table 60 cm deep, recharge 0.05·K | 1.1× (Picard), 1.2× (Newton) faster; unsaturated heads differ by up to 65 cm |
| 16,000 cells | the domain floods | 2–4× slower: cells flip on and off at ψ = 0 |
| 1M cells, 100 m | regional benchmark | 427 s → 404 s (1.06×); heads differ by 0.7 cm |
| 1M cells, 10 m | regional benchmark | 455 s → 488 s (7% slower, 239 instead of 212 outer); heads differ by 1.9 cm |

- **The linear solve did not get cheaper.** The ratio of vertical to horizontal
  conductance is C_v/C_h = (Kv/Kh)(Δx/Δz)², about 1e5 for 100 m cells with 30 cm
  layers and 1e3 for 10 m cells. Horizontal coupling in the unsaturated zone is
  already weak, so removing it changes little.
- **The saturation switch costs outer iterations.** Smoothing it over a band of
  0.5 or 2.5 cm did not help: in the flood case ε = 0.5 cm failed to converge and
  ε = 2.5 cm made Newton 20× slower.
- **Condensing the columns** (eliminating each unsaturated tridiagonal column
  with a Thomas sweep before the solve) would save about 7–8%: only about 15% of
  the cells are unsaturated, the preconditioner already handles the vertical
  coupling, and the set of condensed cells changes every outer iteration.

### Column ordering

**Result: no gain; inner iteration counts were identical. Dropped.**

Numbering the cells column by column instead of layer by layer should let the
incomplete factorisation keep more of the vertical coupling. The same 1M-cell
models were built as DISU grids in both orderings and solved with MILU0.

| Case (1M cells) | Layer order: time (s) / inner | Column order: time (s) / inner |
| --- | --- | --- |
| Loam `shallow` | 118 / 713 | 120 / 725 |
| Loam `deep` | 589 / 2688 | 659 / 2690 |
| Silt loam `shallow` | 100 / 421 | 105 / 414 |
| Silt loam `deep` | 210 / 605 | 210 / 606 |
| Haverkamp `rain20` | 135 / 911 | 137 / 914 |

On a 7-point stencil, ILU0 in layer order already treats each vertical
connection exactly; only fill from mixed horizontal–vertical products is
dropped, and that fill is the same in either order. A line or block-column
preconditioner would also only reduce inner iterations, of which there are
already few.

### Limiting the Newton Jacobian term

**Result: not a speed-up. Full Newton already gives the fewest outer iterations;
damping the term costs 0–25% more and near saturation it breaks convergence.
Dropped.**

With upstream weighting, the flow from an upstream neighbour m into cell n is
Q = kr(ψ_m) · C · (h_m − h_n). Its Newton derivative has an extra term
dkr/dψ · C · Δh, which for inflow from above lands on the *off-diagonal* of row
n. It can exceed the diagonal, so the matrix is no longer diagonally dominant.

The ratio |dkr/dψ · Δh| / kr (100k cells) is 3.6 in an easy case, 9–28 at a
wetting front entering dry soil (loam and silt loam `deep`), 35–62 close to
saturation (loam `wet`, ψ ≈ −0.005 cm, where dkr/dψ is singular for n < 2), and
unbounded in runs that diverge (sand, Haverkamp `deep`). Scaling the term (θ) or
capping it changes only the convergence path, not the answer:

| Case (100k cells, outer iterations) | Full Newton | θ = 0.5 | θ = 0 | cap 3 |
| --- | --- | --- | --- | --- |
| Loam `shallow` | **194** | 210 | 246 | 197 |
| Loam `deep` | 561 (2 failed) | 546 | 625 | 572 |
| Loam `wet` | **406** | 429 | timeout | 1784 |
| Haverkamp `shallow` | 400 | **383** | 771 | 390 |
| Sand `shallow` | fails | fails | 2469 | **1123** |

The large term is Newton doing its job; near saturation it is needed. The
diagnostics did show where runs go wrong, which led to the next three methods.

### Under-relaxation

**Result: the default delta-bar-delta (DBD) under-relaxation over-damps normal
cases by 30–45%; switching it off breaks the hard cases. DBD with fast recovery
(UNDER_RELAXATION_KAPPA 0.2) keeps most of the gain; at 1M cells it needs
[kr smoothing](#smoothing-kr-near-saturation) to stay robust. An input
setting.**

COMPLEX turns on DBD with θ = 0.8 and κ = 1e-4. DBD keeps a relaxation factor
per cell: on every sign change of the cell's Newton step it multiplies the factor
by θ, otherwise it adds κ. Near convergence Newton normally overshoots slightly
and corrects back, so almost every cell sees a sign change. With κ = 1e-4 the
factor hardly recovers, the rest of the time step runs at 0.8× steps, and
convergence becomes linear instead of quadratic. (UNDER_RELAXATION_GAMMA has no
effect unless momentum is used.)

| Case (100k cells, outer iterations) | DBD default | No under-relaxation | DBD κ = 0.2 | SIMPLE 0.7 | Backtracking 20 |
| --- | --- | --- | --- | --- | --- |
| Loam `shallow` | 194 | **132** | 153 | 301 | 194 |
| Silt loam `deep` | 242 | **130** | 161 | 357 | 242 |
| Haverkamp `rain20` | 221 | **148** | 166 | 305 | 221 |
| Haverkamp `shallow` | 400 | **247** | 268 | 449 | 1398 |
| Haverkamp `layered` | 494 | 701 | **299** | 459 | 997 |
| Loam `wet` | 406 | fails | **274** | timeout | fails |
| Loam `deep` | 561 | 675 (52k inner) | 417 (42k inner) | 392 | **323** |
| Sand `shallow` | fails | timeout | fails | fails | **1745** |

- **Normal cases want no damping**: 32–46% fewer outer iterations without it.
- **The loam `wet` failure is a 2-cycle**: a saturated top-layer cell flips
  ±28 cm every iteration, most likely the surface drain switching. DBD exists to
  damp this, and with fast recovery it does so without slowing everything else.
- **Wetting fronts into dry soil (loam `deep`, sand) need their steps limited.**
  No global setting handles all three, because the trouble is local; that is
  what the step limiter adds.

### Saturation step limiter

**Result: with fast-recovering DBD, 41% fewer outer iterations at 100k cells over
the cases that converged before, up to 3× fewer in loam `deep`, and sand
converges. VSF option `NEWTON_DSAT_MAX`, on by default at 0.2.**

After each linear solve, VSF checks the proposed head change of every
unsaturated cell. When the change in effective saturation would exceed
`NEWTON_DSAT_MAX`, the step of that cell is scaled down proportionally. Because
the retention curve is curved, the resulting saturation change is close to, but
not exactly, the limit (for example 0.299 instead of 0.2 for a Haverkamp step
from ψ = −30 to 0 cm); the option is a damping, not a hard bound. Cells that
stay saturated, and cells whose saturation hardly changes, are never touched,
so normal cases keep full Newton steps. The limiter uses the Newton
under-relaxation hook in IMS (the one NPF uses to keep heads above the cell
bottom), but does not raise the under-relaxation flag: that flag lets IMS accept
an outer iteration in which a package has not converged (see the upstream
issue on the Newton under-relaxation convergence shortcut). The regular
convergence check therefore applies unchanged. It costs two
saturation evaluations per cell per outer iteration, about 3% of the run. A
separate limit on drying steps made no measurable difference and was dropped.

| Case (100k cells, outer iterations) | DBD default (reference) | ΔSe 0.2, no under-relaxation | ΔSe 0.2 + DBD default | **ΔSe 0.2 + DBD κ = 0.2** |
| --- | --- | --- | --- | --- |
| Loam `shallow` | 194 | 132 | 194 | **153** |
| Silt loam `deep` | 242 | 129 | 238 | **159** |
| Haverkamp `rain20` | 221 | 149 | 223 | **166** |
| Haverkamp `shallow` | 400 | 253 | 406 | **267** |
| Haverkamp `layered` | 494 | 265 | 426 | **282** |
| Loam `deep` | 561 (2 failed) | 464 (3 failed) | 328 | **187** (0 failed) |
| Loam `wet` | 406 | fails (2-cycle) | 371 | **269** |
| Sand `shallow` | fails | 343 | 703 | **377** |
| Haverkamp `deep` | fails | fails | fails | fails |

The limiter alone gives the lowest counts in smooth cases but cannot damp the
drain 2-cycle in loam `wet`; with fast-recovering DBD every case except
Haverkamp `deep` converges. Head differences against the reference are at most
0.0013 cm, except loam `deep` (0.31 cm), where the reference itself needed two
failed attempts and took different time steps. Haverkamp `deep` diverges to heads
of about 1e26 in the first step in every configuration and then crashes in
BiCGSTAB (see [open items](#open-items)).

The limiter is tested in `autotest/test_gwf_vsf_dsat.py`: infiltration into a dry
Haverkamp column gives the same heads with the limit off, at the default and at
0.05, and the tight limit is shown to be active.

### Smoothing kr near saturation

**Result: fixes the one case that fast-recovering DBD broke (silt loam `layered`
at 1M cells: 1915 s with 40 failed attempts → 75 s), while changing heads by at
most a few tenths of a centimetre. Option `KR_SMOOTHING`, off by default. An
air-entry value (`HAEVGN`) fixes it too, but changes the soil model.**

**The failure.** The step limiter with fast DBD worked at 100k cells, but at 1M
cells silt loam `layered` took 1915 s instead of 99 s. The top-layer cells
oscillated by ±15–20 cm, and the oscillation hopped from cell to cell. The cause
is the Mualem–Van Genuchten kr curve for n < 2:

- To pass 5% of Ks through a layer with Ks/10, kr there must be about 0.5. For
  silt loam (n = 1.41) that is reached at ψ ≈ −2 cm, where Se is still 0.997.
- Within those first 2 cm below saturation kr falls from 1.0 to 0.54, and its
  slope goes to infinity at ψ = 0. The low-K layers sit right on this cusp.
- A Newton step of ±20 cm swings kr between about 1 and 0.12, so Newton
  overshoots back and forth. The saturation limiter cannot see it (ΔSe ≈ 0.003),
  and the default DBD only hid it by keeping these cells damped.

**Two remedies, both options of the Van Genuchten model:**

- **`KR_SMOOTHING <width>`** (VSF options, model length units). Between
  ψ = −width and 0, kr is replaced by a cubic that joins the original curve with
  the same value and slope at −width, and reaches kr = 1 with zero slope at
  saturation. For steep curves (large n or a wide interval) the slope at −width
  is capped at 3·(1 − kr(−width))/width, so the cubic stays monotone and never
  exceeds 1 (Fritsch and Carlson, 1980); the curve then has a small kink at
  −width. This cap is not active for the soils in this study. The retention
  curve and kr below −width are unchanged.
- **`HAEVGN`** (per-cell array): an air-entry pressure head h_s < 0 after Vogel
  et al. (2001). The soil is saturated for ψ ≥ h_s, and Se and kr are rescaled
  to reach 1 at h_s. This removes the cusp entirely, but it is a different soil
  description.

| Silt loam `layered`, 1M cells, fast DBD | Time | Outer iterations | Failed attempts |
| --- | --- | --- | --- |
| No regularisation | 1915 s | 4623 | 40 |
| `KR_SMOOTHING` 2 cm | 93 s | 170 | 0 |
| `HAEVGN` −2 cm | 89 s | 159 | 0 |

| Change of the answer vs plain Van Genuchten (100k cells, 8 cases) | Largest head change | Largest saturation change |
| --- | --- | --- |
| `KR_SMOOTHING` 2 cm | 0.23 cm | 3e-4 |
| `HAEVGN` −2 cm | 7–150 cm | 0.03–0.085 |

So `KR_SMOOTHING` works as a numerical regularisation, while `HAEVGN` is a
physical parameter that users would set and calibrate. Neither affects the
Haverkamp model, which has no cusp. `autotest/test_gwf_vsf_vgnearsat.py` checks
both curves directly: in a steady column with unit gradient, the simulated
pressure head must match the head where the modified kr equals q/Ks, within
1e-4 cm.

## Final comparison at 1M cells

**Result: the candidate set-up solves the seven cases the reference also solves
in 725 s instead of 1614 s (2.2× faster), with no failed time steps, and it
solves the sand case on which the reference crashes. The conservative set-up,
which keeps the default DBD, is 1.8× faster.**

Timings measured on different days varied by up to 20% on this machine, so the
three set-ups were run back to back, interleaved case by case, with the same
executable (study tag `final1m`):

- **Reference**: MILU0 without backtracking, default DBD; limiter and smoothing
  off. The best set-up available before the second study.
- **Conservative**: reference + INNER_DVCLOSE 0.01 + `NEWTON_DSAT_MAX` 0.2 (the
  default) + `KR_SMOOTHING` 2 cm.
- **Candidate**: conservative + UNDER_RELAXATION_KAPPA 0.2 (the
  [recommended settings](#recommended-settings)).

| Case (1M cells, Newton) | Reference: s / outer | Conservative: s / outer | Candidate: s / outer | Speed-up, candidate | Largest head difference, candidate |
| --- | --- | --- | --- | --- | --- |
| Loam `deep` | 699 / 952, 4 failed | 227 / 302 | **146 / 193** | 4.8× | 1.1 cm |
| Silt loam `deep` | 216 / 303 | 154 / 205 | **122 / 152** | 1.8× | <0.001 cm |
| Loam `layered` | 187 / 346 | 134 / 249 | **109 / 197** | 1.7× | 0.4 cm |
| Silt loam `layered` | 120 / 245 | 87 / 176 | **75 / 147** | 1.6× | 0.02 cm |
| Haverkamp `rain20` | 133 / 232 | 101 / 184 | **88 / 158** | 1.5× | 0.05 cm |
| Loam `shallow` | 158 / 247 | 119 / 191 | **111 / 163** | 1.4× | 0.14 cm |
| Silt loam `shallow` | 101 / 184 | 81 / 140 | **74 / 125** | 1.4× | 0.18 cm |
| Sand `shallow` | crashes after 841 s | 467 / 704 | **300 / 386** | — | — |
| **Total, first seven** | **1614 s** | **903 s (1.8×)** | **725 s (2.2×)** | | |

The head differences combine the smoothing of kr within 2 cm of saturation (up
to 0.4 cm here) and the different time steps ATS takes when there are no failed
attempts (loam `deep`). Saturation differences are at most 7e-4.

The candidate needs 34–80% fewer outer iterations in every case. The speed-up is
smaller than the iteration reduction because each outer iteration still costs
the same formulate work, now 40–60% of the run.

## Open items

- [ ] Decide whether `KR_SMOOTHING` should be on by default for Van Genuchten
  soils with n < 2. Its default must then be tied to the model length units
  (2 cm = 0.02 m).
- [ ] Document `NEWTON_DSAT_MAX`, `KR_SMOOTHING` and `HAEVGN` in the mf6io guide
  (now only in the `.dfn` descriptions), and echo `KR_SMOOTHING` in the listing
  file.
- [ ] Make IMS report non-convergence instead of crashing when BiCGSTAB overflows
  on a near-singular matrix (`ddot` in `ImsLinearBase.f90`). Today a bad Newton
  matrix aborts the run (SIGFPE) instead of letting ATS retry with a smaller step.
  This is an existing IMS weakness, not caused by VSF.
- [ ] Investigate Haverkamp `deep` (steep soil, dry profile): it diverges in the
  first step in every configuration.
- [ ] Reuse the saturation already computed by the storage term in the step
  limiter, to remove its ~3% cost.
- [ ] Take the backtracking residual from the Newton system instead of an extra
  `inewton=0` build, so backtracking becomes nearly free when it is wanted. In
  MF6's Newton formulation A·x − b at the current heads should equal the Picard
  residual; this has not been checked in detail.
- [ ] Decide what happens to VSF_1D_ONLY (branch `vsf-1d-only`): keep it as a
  development option or drop it.

**On hold, by decision:** faster matrix assembly (both directions of a
connection in one call, skipping the NPF default loops for VSF connections),
which is the largest remaining lever with formulate at 40–60% of the run; and AMG
or PETSc solvers.

**Worth knowing:** the Newton formulation always uses upstream kr between cells,
while Picard uses the KR_AVERAGING option. The two therefore solve slightly
different discrete problems, which matters when comparing them.

## Reproducing the results

Build MODFLOW 6 into `bin/` first (`pixi run build _builddir_Linux_gfortran_release`
or `pixi run build builddir`), then run from `autotest/vsf_benchmark`.

**Regional benchmark** (`vsf_big.py`). Each run takes 2–8 minutes on one core
and writes `stats.json` with step and iteration counts and profiler timings:

```
pixi run python vsf_big.py runs/d10m_ilut5 --delr 1000
pixi run python vsf_big.py runs/d10m_milu0 --delr 1000 --solver milu0
pixi run python vsf_big.py runs/d10m_milu0_nobt --delr 1000 --solver milu0_nobt
```

Omit `--delr` for the 100 m model. `--one-d` needs the VSF_1D_ONLY option from
the `vsf-1d-only` branch. Do not run the `ilu0rcm` variant without a time limit.

**Parametrised study** (`vsf_study.py`). Results go to
`study/<tag>/<case>/` and `study/<tag>/results.jsonl`:

```
pixi run python vsf_study.py list --n 58 --soil loam sand --scenario shallow deep --form newton --dt0 30
pixi run python vsf_study.py run --tag mytag --n 183 --soil loam --scenario deep \
    --form newton --dt0 30 --solver milu0 milu0_dbdk_dv01 --krsmooth 2 --timeout 1800

# compare two tags case by case (wall time, iterations, head and saturation differences)
pixi run python vsf_study.py compare base mytag
```

Options: `--n` (columns per side; 58 ≈ 100k cells, 183 ≈ 1M), `--soil`,
`--scenario`, `--form newton|picard`, `--solver` (the IMS variants named above,
defined in the `SOLVERS` dictionary), `--grid dis|disu_layer|disu_column`,
`--exe` (the MF6 binary) and `--env KEY=VALUE` (extra environment for the run).
Every run also writes `ims_outer.csv` with the per-iteration maximum change and
its cell. The new VSF options are set with:

| Argument | VSF input | Default |
| --- | --- | --- |
| `--dsat 0.2` | `NEWTON_DSAT_MAX` (0 = off) | 0.2 |
| `--krsmooth 2` | `KR_SMOOTHING` (length units) | off |
| `--hae -2` | `HAEVGN`, same value in every cell | 0 (no air entry) |

The final comparison, per soil and scenario:

```
B="pixi run python vsf_study.py run --tag final1m --n 183 --form newton --dt0 30 --timeout 1500"
$B --soil loam --scenario deep --solver milu0 --dsat 0                     # reference
$B --soil loam --scenario deep --solver milu0_dv01 --krsmooth 2            # conservative
$B --soil loam --scenario deep --solver milu0_dbdk_dv01 --krsmooth 2       # candidate
```

Run outputs (`runs/`, `study/`, `study_*.log`) are not in the repository. VS2DT
for the strip comparison was built from source in the richards-modflow6
examples folder.
