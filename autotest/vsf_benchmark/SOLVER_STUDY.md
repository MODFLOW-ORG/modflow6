# Speeding up VSF: solver study and results

Study report, 2026-11-24. The cases are built and run by
[`vsf_study.py`](vsf_study.py) in this directory. A shared, editable copy of this
report lives in Claude Docs:
<https://claude.ai/artifact/KkcCmcZR74oNy3trQCY34H>.

## Summary

VSF is now about as fast as VS2DT on the comparison problem, and at 1 million
cells the right IMS settings make it 3–8× faster than the COMPLEX defaults. On
top of that, a per-cell step limiter, smoothing of kr near saturation and
faster-recovering under-relaxation make the 1M-cell runs another 2.2× faster
overall (1.4–4.8× per case), and make cases converge that failed before.

We tested eight methods. They attack three different costs: building the matrix,
solving it, and the number of Newton (outer) iterations.

| Method | Kind of change | Impact | Status |
| --- | --- | --- | --- |
| 1. Fused soil-function evaluations | code | 4.5× faster on the strip benchmark (5.15 → 1.15 s), on par with VS2DT | committed |
| 2. MILU0 preconditioner, no backtracking | input | 2.6–7.6× faster at 1M cells | recommended |
| 3. INNER_DVCLOSE = OUTER_DVCLOSE | input | 1.3–1.5× faster at 1M cells | recommended |
| 4. Column ordering | input (DISU) | none | dropped |
| 5. Limiting the Newton Jacobian term | code | slower in most cases, only sand gains | dropped |
| 6. DBD under-relaxation with fast recovery (κ = 0.2) | input | 20–45% fewer outer iterations; safe only together with method 8 | recommended |
| 7. Per-cell saturation step limiter | code | with methods 3, 6 and 8: 2.2× faster at 1M cells, sand converges | option `NEWTON_DSAT_MAX`, default 0.2 |
| 8. Smoothing Van Genuchten kr near saturation (and optional air entry) | code | fixes the silt loam `layered` failure (1915 → 75 s); heads change ≤ 0.4 cm | options `KR_SMOOTHING` and `HAEVGN`, off by default |

The final comparison at 1M cells (wall time per run, the three set-ups run back
to back; see [Final comparison](#final-comparison-at-1m-cells)) is 1614 s → 725 s
over the seven cases the reference solves.

One steep-soil case (Haverkamp `deep`) still fails in every configuration, and
IMS crashes rather than retrying when BiCGSTAB overflows. Both are listed under
open issues.

## Background

VSF (variably saturated flow) adds Richards-equation flow to the MODFLOW 6 GWF
model. Each unsaturated cell gets a pressure-dependent relative permeability
kr(ψ) and saturation S(ψ) from a soil model (Van Genuchten or Haverkamp), so
every outer iteration has to evaluate those curves and their derivatives in
every cell.

The trigger for this work was the VS2DT comparison in the richards-modflow6
repository (`examples/vs2d_comparison/05_convergence_and_performance.ipynb`). On
the same small strip problem VS2DT ran in about 1.1 s, while VSF needed 5.2 s
with Newton and 5.4 s with Picard. VS2DT is a lean code: one 5-point stencil,
three soil evaluations per cell per iteration and one SIP sweep per iteration.

Those comparison models are small. Our production models are 1 million cells and
up, so the study below measures at that scale and across soils and wetting
conditions. A run's time splits into three parts that respond to different
fixes:

- **Formulate**: building the matrix, including the soil-function evaluations.
  Usually 40–60% of the run.
- **Linear solve**: the inner BiCGSTAB iterations with their preconditioner.
  Usually 25–50%.
- **Number of outer (Newton) iterations**: multiplies both of the above. In hard
  cases this dominates everything else.

The methods below are grouped by which of these they attack.

## Test setup

All numbers come from one parametrised study script, `vsf_study.py`, which
builds, runs and parses each case. The model is a synthetic regional block:

- n × n columns of 10 m × 10 m cells (units are cm and s), 10–20 thin soil
  layers (30–50 cm) on top of 20 aquifer layers of 5 m, with a heterogeneous
  lognormal aquifer and an aquitard.
- n = 183 gives about 1 million cells; n = 58 about 100 000 cells. Recharge
  varies in space; the edges hold a fixed water table.
- Transient, 5 days in two stress periods, adaptive time stepping (ATS) starting
  at 30 s with a 6 h maximum.

| Soil | Model | Ks (cm/s) | Character |
| --- | --- | --- | --- |
| Loam | Van Genuchten, n = 1.56 | 2.9e-4 | medium |
| Silt loam | Van Genuchten, n = 1.41 | 1.25e-4 | fine, gentle curves |
| Sand | Van Genuchten, n = 2.68 | 8.3e-3 | coarse, steep curves |
| Haverkamp sand | Haverkamp (Celia et al., 1990) | 9.4e-3 | steep kr, classic test soil |

| Scenario | What happens |
| --- | --- |
| `shallow` | water table at 1.5 m, 2 days of rain at 5% of Ks, then drainage to drains at 1.2 m |
| `rain20` | as shallow, with a fixed 20 mm/d rain then 1 mm/d (realistic rates for every soil) |
| `deep` | water table at 9 m; a wetting front travels down through dry soil |
| `layered` | as shallow, with alternating layers 10× less conductive |
| `wet` | rain above Ks on a shallow profile; surface drains take the excess |

What we report per run:

- **Wall time** of the MF6 run. Runs at 1M cells were made one at a time. Some
  100k-cell comparisons ran 4–5 cases in parallel on a 16-core machine; there
  wall times are only roughly comparable, so the iteration counts are the main
  metric.
- **Outer iterations** (Newton iterations, including those of failed ATS
  attempts), **inner iterations** (linear solver), **failed attempts** (time
  steps that did not converge and were retried with a smaller step).
- **Accuracy**: the largest head difference (cm) and saturation difference at the
  end of the run against the reference configuration. Head differences below the
  outer closure of 0.01 cm are noise.

Unless stated otherwise the solver is IMS with the COMPLEX preset, BICGSTAB,
OUTER_DVCLOSE 0.01 cm, RCLOSE 10 STRICT, INNER_DVCLOSE 0.001 and Newton.

## Method 1: fused soil-function evaluations

**Result: VSF Newton went from 4.5× slower than VS2DT to on par (5.15 s → 1.15 s
on the strip benchmark); Picard 5.39 s → 2.47 s.**

Profiling showed the run was dominated by `pow`, `exp` and `log` calls inside
the soil models, not by the solver. Three things caused that:

1. dkr/dψ for the Newton Jacobian was a central difference: two extra kr
   evaluations per cell, each of which called the saturation function again.
2. kr, S and the moisture capacity Cm were computed in separate calls that each
   recomputed the same powers.
3. The storage term re-evaluated S(ψ_old) at the start-of-step head in every
   outer iteration, although it only changes once per time step.

The changes:

- New type-bound procedures `kr_dkrdh` and `sat_capacity` on the soil model
  return kr with its analytical derivative, and S with Cm, from one shared set of
  powers. For Van Genuchten, with x = (α|ψ|)^n, both need only three
  exponentials and two logarithms.
- Close to saturation (ψ > −1e-6) the central difference is kept, because the
  exact dkr/dψ is singular there for n < 2.
- S(ψ_old) is cached per cell and refreshed only when a new time step starts.

The results are unchanged to solver precision. This change is in every run
below; it is the baseline for the rest of the report.

## Method 2: MILU0 preconditioner instead of the COMPLEX default

**Result: 2.6–7.6× faster at 1M cells, with no change in the answer. An input
setting.**

The COMPLEX preset uses ILUT (fill level 5, drop tolerance 1e-4) with
backtracking. We replaced the preconditioner with MILU0 (level 0, no drop
tolerance, relaxation 0.97) and switched backtracking off:

```
COMPLEXITY COMPLEX
LINEAR_ACCELERATION BICGSTAB
PRECONDITIONER_LEVELS 0
PRECONDITIONER_DROP_TOLERANCE 0.0
RELAXATION_FACTOR 0.97
BACKTRACKING_NUMBER 0
```

| Case (1M cells, Newton) | COMPLEX (s) | MILU0 (s) | Speed-up |
| --- | --- | --- | --- |
| Loam `shallow` | 478 | 126 | 3.8× |
| Loam `layered` | 1166, failed | 154 | — |
| Loam `deep` | 867 | 632 | 1.4× |
| Silt loam `shallow` | 349 | 90 | 3.9× |
| Silt loam `deep` | 552 | 212 | 2.6× |
| Silt loam `layered` | 537 | 99 | 5.4× |
| Haverkamp `rain20` | 445 | 115 | 3.9× |

Why it works: the soil layers are thin and their vertical conductances dominate,
so the matrix is close to a set of tridiagonal columns. ILU0 in natural
(layer-major) ordering already captures that coupling, and its cost per inner
iteration grows linearly with the cell count. ILUT buys little extra convergence
for a much larger factorisation, and its backtracking adds a full matrix rebuild
every outer iteration. With MILU0 the inner iteration count per outer iteration
stayed flat (about 2–4) from 100k to 1M cells.

Loam `deep` gains least because it is dominated by outer iterations (952 with
MILU0), which methods 6–8 address.

## Method 3: looser inner (linear) tolerance

**Result: setting INNER_DVCLOSE equal to OUTER_DVCLOSE (0.01 cm) is 1.3–1.5×
faster at 1M cells, with head differences up to 0.5 cm in the hardest case and
below 0.05 cm elsewhere. An input setting.**

Each Newton step only needs to point in roughly the right direction, so solving
the linear system to 0.001 cm is more than necessary. We compared three looser
variants with MILU0 at 1M cells:

- `dv01`: INNER_DVCLOSE 0.01 (= outer closure).
- `rel1` / `rel2`: stop the inner solve when the residual norm has dropped by a
  factor 10 or 100 (RCLOSE 0.1 or 0.01 with RELATIVE_RCLOSE).

| Case (1M cells) | Reference (s) | `dv01` (s) | `rel1` (s) | `rel2` (s) | Largest head difference, `dv01` (cm) |
| --- | --- | --- | --- | --- | --- |
| Loam `shallow` | 126 | **87** | 87 | 86 | 0.002 |
| Loam `deep` | 632 | **443** | 537 | 448 | 0.52 |
| Silt loam `shallow` | 90 | **60** | 67 | 64 | 0.006 |
| Silt loam `deep` | 212 | **143** | 154 | 154 | <0.001 |
| Haverkamp `rain20` | 115 | **85** | 135 | 97 | 0.05 |

`dv01` also needed 20–25% fewer *outer* iterations. The reason is the STRICT
criterion: an outer iteration only counts as converged when the linear solve
converges in one inner iteration, which happens sooner with a looser inner
tolerance. So part of the gain is IMS accepting convergence earlier. The
relative-residual variants bring the inner work down to about one iteration per
outer iteration, but `rel1` then needs more outer iterations (Haverkamp: 366
against 232).

We considered an adaptive (Eisenstat–Walker) inner tolerance in IMS but did not
build it: with only 1–4 inner iterations per outer iteration left, it could at
best match `dv01`.

## Method 4: column ordering of the unknowns

**Result: no gain. Inner iteration counts were identical, so this was dropped.**

The idea: number the cells column by column instead of layer by layer, so the
incomplete factorisation keeps more of the strong vertical coupling. We built the
same 1M-cell models as DISU grids in both orderings and solved them with MILU0.

| Case (1M cells) | Layer order: time (s) / inner | Column order: time (s) / inner |
| --- | --- | --- |
| Loam `shallow` | 118 / 713 | 120 / 725 |
| Loam `deep` | 589 / 2688 | 659 / 2690 |
| Silt loam `shallow` | 100 / 421 | 105 / 414 |
| Silt loam `deep` | 210 / 605 | 210 / 606 |
| Haverkamp `rain20` | 135 / 911 | 137 / 914 |

In hindsight this is expected: on a 7-point stencil, ILU0 in layer order already
treats each vertical connection exactly; only fill from mixed
horizontal–vertical products is dropped, and that fill is the same size in
either order. The DISU layer-order runs reproduced the DIS results exactly. A
line or block-column preconditioner would also only reduce inner iterations, of
which there are already few, so we did not pursue it.

## Method 5: limiting the Newton Jacobian term

**Result: not a speed-up. Full Newton already gives the fewest outer iterations;
damping the term costs 0–25% more, and near saturation it breaks convergence.
Only sand benefits. Dropped; the experiment is kept as
`study/method5_jacobian_experiment.patch` (not in the repository).**

With upstream weighting, the flow between cell n and an upstream neighbour m is
Q = kr(ψ_m) · C · (h_m − h_n). Its Newton derivative has an extra term
dkr/dψ · C · Δh. For an inflow from above this term lands on the *off-diagonal*
of row n. It can then exceed the diagonal, so the matrix is no longer diagonally
dominant, which is hard for BiCGSTAB/ILU and can make Newton overshoot.

First we measured how large the term is (100k cells; the ratio is
|dkr/dψ · Δh| / kr):

| Case | Typical ratio | Where |
| --- | --- | --- |
| Loam `shallow` (easy) | 3.6 | moderately wet soil |
| Silt loam / loam `deep` | 9–28 | wetting front entering dry soil (ψ ≈ −500 cm below −50 cm) |
| Loam `wet` | 35–62 | close to saturation (ψ ≈ −0.005 cm, where dkr/dψ is singular for n < 2) |
| Sand, Haverkamp `deep` | 37 → ∞ | runs that diverge (heads at −4000 m) |

The term exceeds kr in thousands of rows even in the easy case, which converges
in about 5 outer iterations per step anyway. We then tried scaling the term
(θ = 0.5 and 0) and capping it at 1× and 3× kr. The residual is untouched, so the
converged answer is the same; only the path changes.

| Case (100k cells, outer iterations) | Full Newton | θ = 0.5 | θ = 0 | cap 3 |
| --- | --- | --- | --- | --- |
| Loam `shallow` | **194** | 210 | 246 | 197 |
| Loam `deep` | 561 (2 failed) | 546 | 625 | 572 |
| Loam `wet` | **406** | 429 | timeout | 1784 |
| Haverkamp `shallow` | 400 | **383** | 771 | 390 |
| Sand `shallow` | fails | fails | 2469 | **1123** |

Capping makes the wetting-front linear systems easier (loam `deep`: 2750 → 1600
inner iterations) but does not reduce outer iterations. Near saturation the term
is needed: removing it makes loam `wet` 4× slower or fail. So the large Jacobian
term is not the problem; it is Newton doing its job. The diagnostics did show
where runs go wrong, which led to methods 6–8.

## Method 6: IMS under-relaxation and backtracking

**Result: the default delta-bar-delta (DBD) under-relaxation over-damps normal
cases by 30–45%; switching it off breaks the hard cases. Changing DBD to recover
quickly (UNDER_RELAXATION_KAPPA 0.2) keeps most of the gain; at 1M cells it needs
the kr smoothing of method 8 to stay robust. An input setting.**

COMPLEX turns on DBD with θ = 0.8 and κ = 1e-4. DBD keeps a relaxation factor per
cell. Every time the Newton step of a cell changes sign it multiplies the factor
by θ; when the sign stays the same it adds κ. Near convergence Newton normally
overshoots slightly and corrects back, so almost every cell sees a sign change.
With κ = 1e-4 the factor then hardly recovers, the rest of the time step runs at
0.8× steps, and convergence becomes linear instead of quadratic.
(UNDER_RELAXATION_GAMMA has no effect unless momentum is used; we confirmed it
changes nothing.)

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

What this shows:

- **Normal cases want no damping**: 32–46% fewer outer iterations without it.
- **The loam `wet` failure is a 2-cycle**: one saturated top-layer cell flips
  ±28 cm every iteration, most likely the surface drain switching on and off. DBD
  exists to damp exactly this, and with fast recovery it does so without slowing
  everything else.
- **Wetting fronts into dry soil (loam `deep`, sand) need their steps limited**;
  only backtracking manages that here, and it is very costly elsewhere
  (Haverkamp `shallow` 3.5× slower).

No global setting handles all three, because the trouble is local. That is what
method 7 adds.

## Method 7: per-cell saturation step limiter in VSF

**Result: combined with fast-recovering DBD, 41% fewer outer iterations at 100k
cells over the cases that converged before, up to 3× fewer in loam `deep`, and
sand now converges. VSF option `NEWTON_DSAT_MAX`, on by default at 0.2.**

After each linear solve, VSF checks the proposed head change of every
unsaturated cell. When the change in effective saturation would exceed
`NEWTON_DSAT_MAX` (0.2 works well), the step of that cell is scaled down
proportionally. Cells that stay saturated, and cells whose saturation hardly
changes, are never touched, so normal cases keep full Newton steps. The limiter
uses the existing Newton under-relaxation hook in IMS (the one NPF uses to keep
heads above the cell bottom), so convergence checking and reporting work
unchanged. It costs two saturation evaluations per cell per outer iteration
(about 3% of the run). A second, drying-specific limit was also tried; it made no
measurable difference and was dropped.

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

The limiter alone (column 2) gives the lowest counts in smooth cases but cannot
damp the drain 2-cycle in loam `wet`. With fast-recovering DBD (last column) every
case except Haverkamp `deep` converges, at 1.3–3.0× fewer outer iterations than
the reference.

**Accuracy.** Head differences against the reference are 0.0013 cm or less in
all cases except loam `deep` (0.31 cm). There the reference itself needed two
failed attempts and different time steps, so the difference reflects a different
ATS path rather than a less accurate answer.

**Still failing: Haverkamp `deep`.** In every configuration heads diverge to
about 1e26 in the first time step, and the run then crashes with a
floating-point overflow inside BiCGSTAB (see open issues).

## Method 8: smoothing Van Genuchten kr near saturation

**Result: fixes the one case that fast-recovering DBD broke (silt loam `layered`
at 1M cells: 1915 s with 40 failed attempts → 75 s), while changing heads by at
most a few tenths of a centimetre. New option `KR_SMOOTHING`, off by default. An
air-entry value (`HAEVGN`) fixes it too, but changes the soil model.**

**The failure.** Methods 6 and 7 together worked at 100k cells, but at 1M cells
silt loam `layered` took 1915 s instead of 99 s. The top-layer cells oscillated
by ±15–20 cm, and the oscillation hopped from cell to cell. The cause is the
shape of the Mualem–Van Genuchten kr curve for n < 2:

- To pass 5% of Ks through a layer with Ks/10, kr there must be about 0.5. For
  silt loam (n = 1.41) that is reached at ψ ≈ −2 cm, where Se is still 0.997.
- Within those first 2 cm below saturation kr falls from 1.0 to 0.54, and its
  slope goes to infinity at ψ = 0. The low-K layers sit right on this cusp.
- A Newton step of ±20 cm then swings kr between about 1 and 0.12, so Newton
  overshoots back and forth. The saturation limiter cannot see it (ΔSe ≈ 0.003),
  and the default DBD only hid it by keeping these cells damped.

The same cusp is behind the large Jacobian terms seen near saturation in
method 5. Any fine-textured soil carrying flux through a less permeable layer can
meet it.

**Two remedies, both input options for the Van Genuchten model:**

- **`KR_SMOOTHING <width>`** (VSF options, model length units). Between
  ψ = −width and 0, kr is replaced by a cubic that joins the original curve with
  the same value and slope at −width, and reaches kr = 1 with zero slope at
  saturation. The retention curve and kr below −width are unchanged, and the
  slope of kr is bounded everywhere.
- **`HAEVGN`** (per-cell array): an air-entry pressure head h_s < 0 after Vogel
  et al. (2001). The soil is saturated for ψ ≥ h_s, and Se and kr are rescaled so
  they reach 1 at h_s. This removes the cusp entirely, but it is a different soil
  description.

| Silt loam `layered`, 1M cells, fast DBD | Time | Outer iterations | Failed attempts |
| --- | --- | --- | --- |
| No regularisation | 1915 s | 4623 | 40 |
| `KR_SMOOTHING` 2 cm | 93 s | 170 | 0 |
| `HAEVGN` −2 cm | 89 s | 159 | 0 |

How much each changes the answer, against plain Van Genuchten (100k cells, 8
cases):

| | Largest head change | Largest saturation change |
| --- | --- | --- |
| `KR_SMOOTHING` 2 cm | 0.23 cm | 3e-4 |
| `HAEVGN` −2 cm | 7–150 cm | 0.03–0.085 |

So `KR_SMOOTHING` works as a numerical regularisation, while `HAEVGN` is a
physical parameter that users would set and calibrate. Neither changes the
Haverkamp model, which has no cusp. The test
`autotest/test_gwf_vsf_vgnearsat.py` checks both curves directly: in a steady
column with unit gradient, the simulated pressure head must match the head where
the modified kr equals q/Ks, within 1e-4 cm.

## Final comparison at 1M cells

**Result: the candidate set-up solves the seven cases the reference also solves
in 725 s instead of 1614 s (2.2× faster), with no failed time steps, and it
solves the sand case on which the reference crashes. The conservative set-up,
which keeps the default DBD, is 1.8× faster.**

Timings measured on different days varied by up to 20% on this machine, so the
three set-ups were run back to back, interleaved case by case, with the same
executable (study tag `final1m`):

- **Reference**: MILU0 (method 2) with the default DBD; limiter and smoothing
  off. This is the best set-up available before this work.
- **Conservative**: reference + INNER_DVCLOSE 0.01 + `NEWTON_DSAT_MAX` 0.2 (the
  default) + `KR_SMOOTHING` 2 cm.
- **Candidate**: conservative + fast-recovering DBD
  (UNDER_RELAXATION_KAPPA 0.2).

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

The head differences combine two effects: the smoothing of kr within 2 cm of
saturation, which changes the answer slightly (up to 0.4 cm here), and the
different time steps ATS takes when there are no failed attempts (loam `deep`).
Saturation differences are at most 7e-4.

The candidate needs 34–80% fewer outer iterations than the reference in every
case. The speed-up is smaller than the iteration reduction because each outer
iteration still costs the same formulate work, which is now 40–60% of the run.

## Recommendations and open issues

**Recommended input** for VSF models with Newton (the candidate set-up):

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

This is 2.2× faster than MILU0 alone at 1M cells. Without
UNDER_RELAXATION_KAPPA 0.2 (the conservative set-up) it is 1.8× faster. DVCLOSE
values and `KR_SMOOTHING` are in model length units (cm here) and scale with
them.

**Open items:**

- [ ] Decide whether `KR_SMOOTHING` should be on by default for Van Genuchten
  soils with n < 2. Its default must then be tied to the model length units
  (2 cm = 0.02 m).
- [ ] Document the three options in the mf6io guide (now only in the `.dfn`
  descriptions), and echo `KR_SMOOTHING` in the listing file.
- [ ] Make IMS report non-convergence instead of crashing when BiCGSTAB overflows
  on a near-singular matrix (`ddot` in `ImsLinearBase.f90`). Today a bad Newton
  matrix aborts the run (SIGFPE) instead of letting ATS retry with a smaller
  step. This is an existing IMS weakness, not caused by VSF.
- [ ] Investigate Haverkamp `deep` (steep soil, dry profile): it diverges in the
  first step in every configuration.
- [ ] Reuse the saturation already computed by the storage term in the limiter,
  to remove its ~3% cost.

**Not recommended:** column ordering (method 4) and limiting the Newton Jacobian
(method 5). An adaptive inner tolerance (Eisenstat–Walker) has little left to
gain after method 3.

**On hold, by decision:**

- Faster matrix assembly: filling both directions of a connection in one call,
  and skipping the NPF default loops for VSF connections. Formulate is still
  40–60% of the run, so this is the largest remaining lever.
- AMG or PETSc solvers.

**Worth knowing:** the Newton formulation always uses upstream kr between cells,
while Picard uses the KR_AVERAGING option (arithmetic in this study). The two
therefore solve slightly different discrete problems, which matters when
comparing them.

## Reproducing the results

Branch `vsf-dev-optimization` (from `vsf-dev`). Run from
`autotest/vsf_benchmark` after `pixi run build _builddir_Linux_gfortran_release`:

```
# list or run cases; results go to study/<tag>/<case>/ and study/<tag>/results.jsonl
pixi run python vsf_study.py list --n 58 --soil loam sand --scenario shallow deep --form newton --dt0 30
pixi run python vsf_study.py run --tag mytag --n 183 --soil loam --scenario deep \
    --form newton --dt0 30 --solver milu0 milu0_dbdk_dv01 --krsmooth 2 --timeout 1800

# compare two tags case by case (wall time, iterations, head and saturation differences)
pixi run python vsf_study.py compare base mytag
```

Options: `--n` (columns per side; 58 ≈ 100k cells, 183 ≈ 1M), `--soil`,
`--scenario`, `--form newton|picard`, `--solver` (the IMS variants named in this
report, defined in the `SOLVERS` dictionary), `--grid dis|disu_layer|disu_column`,
`--exe` (the MF6 binary) and `--env KEY=VALUE` (extra environment for the run).
Every run also writes `ims_outer.csv` with the per-iteration maximum change and
its cell. The new VSF options are set with:

| Argument | VSF input | Default |
| --- | --- | --- |
| `--dsat 0.2` | `NEWTON_DSAT_MAX` (method 7; 0 = off) | 0.2 |
| `--krsmooth 2` | `KR_SMOOTHING` (method 8, length units) | off |
| `--hae -2` | `HAEVGN`, same value in every cell (method 8) | 0 (no air entry) |

The final comparison is reproduced with, per soil and scenario:

```
B="pixi run python vsf_study.py run --tag final1m --n 183 --form newton --dt0 30 --timeout 1500"
$B --soil loam --scenario deep --solver milu0 --dsat 0                     # reference
$B --soil loam --scenario deep --solver milu0_dv01 --krsmooth 2            # conservative
$B --soil loam --scenario deep --solver milu0_dbdk_dv01 --krsmooth 2       # candidate
```

Study results are not in the repository (`study/` is ignored). VS2DT for the
strip comparison was built from source in the richards-modflow6 examples folder.
