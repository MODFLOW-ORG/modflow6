# VSF_1D_ONLY and solver performance for large VSF models

Benchmark report, 2026-11-23. The benchmark model is built by
[`vsf_big.py`](vsf_big.py) in this directory.

## Summary

Restricting Richards flow to vertical connections (VSF_1D_ONLY) gives at most
1.06–1.2× and can be slower. Switching the IMS preconditioner from ILUT(5) to
MILU0 makes a million-cell model 2.6× faster.

- **VSF_1D_ONLY is not worth it for performance.** On a 1-million-cell regional
  model it gave 1.06× with 100 m cells and was 7% slower with 10 m cells. The
  linear solve did not get cheaper.
- **The preconditioner is the big lever.** MILU0 (the MODERATE linear preset) cut
  the linear solve from 326 s to 49 s and wall time from 455 s to 177 s. Heads
  stayed within 0.003 cm.
- **Soil functions dominate after that.** The COMPLEX backtracking check never
  triggered but rebuilt the system every outer iteration. Turning it off brings
  the 10 m model to 125 s (3.6× in total). Formulate is then just over half the
  run. This supports optimising the soil parametrisations next.
- **Condensing the 1D columns is a dead end.** It would save about 7–8% at best.

## The VSF_1D_ONLY option

The new VSF option applies Richards flow only to vertical connections (ihc = 0).
A horizontal connection conducts flow only when both cells are fully saturated.

- **Saturation criterion.** A cell counts as fully saturated when its pressure
  head is non-negative (psi = h − z_node ≥ 0, z_node = ½(top + bot)). The test is
  at the cell node, not the cell top.
- **Horizontal conductance.** A horizontal connection uses the saturated
  conductance (`condsat`) when both cells are saturated, and zero otherwise.
  Partially saturated cells get no horizontal flow and no smoothing.
- **Where it lives.** The keyword is in `gwf-vsf.dfn` (IDM regenerated). It is
  read in `gwf-vsf.f90` as `ivsf1d` and passed to `VsfFlowType%initialize`.
  `VsfFlow.f90` routes horizontal connections to `calculate_coeffs_horizontal`
  in both the Picard and the Newton path.
- **Tests** (`autotest/test_gwf_vsf_1d.py`), with Picard and Newton:
  - `test_uniform_columns`: identical columns give the same heads as full VSF
    (atol 1e-6 cm).
  - `test_horizontal_flow_saturated_only`: an xz cross-section has no horizontal
    flowja where either cell has psi < 0, has saturated horizontal flow
    elsewhere, and closes its budget.
  - `test_benchmark` (marked slow): a parametrised timing comparison.

All 4 new tests and the 35 existing VSF tests pass.

## Small benchmark and smoothing

On a 16,000-cell test model, VSF_1D_ONLY gave 1.1–1.2× when the water table
stays deep, and was 2–4× slower when the model floods.

| Scenario | Picard | Newton | Notes |
| --- | --- | --- | --- |
| Realistic: water table 60 cm deep, recharge 0.05·K, 12 h | 1.1× faster | 1.2× faster | Unsaturated heads differ by up to 65 cm |
| Near-saturated: the domain floods | 2–4× slower | 2–4× slower | Cells flip on and off at psi = 0 |

The slowdown comes from the on/off switch at psi = 0. Each time a cell flips,
the horizontal conductance jumps, which costs extra outer iterations.

**Smoothing does not fix this.** We tested a smoothing band ε of 0.5 and 2.5 cm
on the saturated side. With ε = 2.5 cm, Newton recovered its outer iterations in
the realistic case. In the flood case, ε = 0.5 cm failed to converge and
ε = 2.5 cm made Newton 20× slower. Smoothing was removed, and the option uses
the hard criterion.

## Million-cell regional model

On a 1-million-cell model, VSF_1D_ONLY did not reduce the linear solve. It was
1.06× faster with 100 m cells and 7% slower with 10 m cells.

**Model setup** (cm and s, 183 × 183 × 30 = 1.0 million cells):

- **Layers.** 10 fine layers of 30 cm for the unsaturated zone, over 20 coarse
  aquifer layers of 5 m.
- **Conductivity.** A lognormal K field in the aquifer, an aquitard in coarse
  layers 7–8 (K = 1e-5 cm/s), and Celia/Haverkamp sand on top.
- **Stresses.** Patchy recharge of 20 mm/d for 2 days, then 1 mm/d for 3 days.
  Drains sit at 1.2 m depth, and fixed heads hold the lateral boundary at 1.5 m
  depth.
- **Solver.** Newton, ATS with a maximum step of 6 h, and IMS COMPLEX with
  BICGSTAB (ILUT level 5, drop tolerance 1e-4).

All runs use Newton with the ILUT(5) preconditioner:

| Cell size | Variant | Wall time | Linear solve | Outer / inner iterations | Max head difference |
| --- | --- | --- | --- | --- | --- |
| 100 m | Full VSF | 427 s | 214 s | | |
| 100 m | VSF_1D_ONLY | 404 s | ≈214 s | | 0.7 cm |
| 10 m | Full VSF | 455 s | 326 s | 212 / 1090 | |
| 10 m | VSF_1D_ONLY | 488 s | | 239 outer | 1.9 cm |

At 100 m, 14.7% of cells were unsaturated. The linear solve took 50% of the run
and formulate 22% (94 s).

**Cell size controls the horizontal coupling.** The ratio of vertical to
horizontal conductance is C_v/C_h = (Kv/Kh)(Δx/Δz)². That is about 1e5 for
100 m cells with 30 cm layers and about 1e3 for 10 m cells. At both sizes,
horizontal coupling in the unsaturated zone is already weak, so removing it
changes little. In the 10 m case, the saturation switch added outer iterations
that outweighed the savings.

## Condensing the 1D columns

Eliminating the unsaturated columns before the linear solve would save about
7–8% of the run, so we did not implement it.

With VSF_1D_ONLY, each unsaturated column is a tridiagonal chain attached to the
saturated system. Those cells could be eliminated with a Thomas sweep, giving a
Schur complement on the saturated cells. Only about 15% of the cells were
unsaturated, though, and the preconditioner already handles strong vertical
coupling well. The gain is therefore small, and the code would get more complex:
the set of condensed cells changes every outer iteration as the water table
moves.

## Preconditioner comparison

MILU0 is 2.6× faster overall and 6.6× faster in the linear solve than the
COMPLEX default ILUT(5). These runs use the 10 m model with full VSF and Newton.

| Preconditioner | Wall time | Linear solve | Formulate | Outer / inner iterations | Linear solve per inner iteration | Max head difference vs ILUT(5) |
| --- | --- | --- | --- | --- | --- | --- |
| **MILU0 (relaxation 0.97)** | **177 s** | **49 s** | 62 s | 211 / 855 | 0.06 s | 0.003 cm |
| ILU0 | 212 s | 81 s | 63 s | 216 / 1508 | 0.05 s | 0.005 cm |
| ILUT(1, 1e-3) | 259 s | 109 s | 73 s | 252 / 2812 | 0.04 s | 0.04 cm |
| ILUT(5, 1e-4), the COMPLEX default | 455 s | 326 s | 62 s | 212 / 1090 | 0.30 s | |
| ILU0 with RCM reordering | failed after 8.6 h | 27,230 s | 704 s | 2294 outer | ≈12 s per solve | |

- **ILUT(5) costs more than it saves.** Each inner iteration costs 5–6× more than
  with ILU0 or MILU0. MILU0 even needs fewer inner iterations.
- **RCM reordering breaks convergence** for this problem and should be avoided.
- **MILU0 is the default linear setup of `COMPLEXITY MODERATE`.** To keep the
  COMPLEX nonlinear settings, set the linear part explicitly:

```
PRECONDITIONER_LEVELS 0
PRECONDITIONER_DROP_TOLERANCE 0.0
RELAXATION_FACTOR 0.97
```

With MILU0 the run splits into linear solve 28% (49 s), formulate 35% (62 s),
and other work 37% (≈60 s). That last part is the backtracking residual check,
covered in the next section.

## Backtracking

Backtracking never triggered in these runs, yet checking for it cost a third of
the MILU0 run. Turning it off took the 10 m model from 177 s to 125 s, with
bit-identical heads.

| 10 m model, MILU0, Newton | Wall time | Solve | Formulate | Linear solve | Outer / inner iterations |
| --- | --- | --- | --- | --- | --- |
| COMPLEX backtracking (20 tries) | 177 s | 170 s | 62 s | 49 s | 211 / 855 |
| `BACKTRACKING_NUMBER 0` | 125 s | 118 s | 64 s | 53 s | 211 / 855 |

- **Settings.** COMPLEX sets `BACKTRACKING_NUMBER 20`, tolerance 1.05, reduction
  factor 0.1 and residual limit 0.002. SIMPLE and MODERATE switch backtracking
  off.
- **How it works** (`sln_backtracking` in `NumericalSolution.f90`). Every outer
  iteration rebuilds the full system (Picard terms, `inewton=0`) only to compute
  the residual L2 norm. A step is cut back by 0.1, with a rebuild each time, only
  if the residual grew more than 5% and 0.1·max|Δh| ≥ `dvclose`.
- **Observed.** The backtrack flag and count were 0 in all 211 outer iterations
  of the ILUT(5) and MILU0 runs. Newton under-relaxation (θ = 0.8) and ATS keep
  the iterations well behaved. Backtracking can still matter for sharp wetting
  fronts or larger time steps.
- **Combined effect.** MILU0 plus `BACKTRACKING_NUMBER 0` takes the 10 m model
  from 455 s to 125 s (3.6×). Formulate is then just over half the run.

**Possible improvement: make the check nearly free.** In MF6's Newton
formulation, the Newton terms are added to both the matrix and the right-hand
side. At the current heads, A·x − b should therefore equal the Picard residual.
If so, the residual can come from the Newton system that is built anyway, and
the extra build is needed only when a backtrack actually happens. That keeps
backtracking as a safety net at almost no cost. This has not yet been checked
against the code in detail.

## Next steps

The next gains are in evaluating the soil functions, now that MILU0 has made the
linear solve cheap.

- [ ] **Optimise the soil parametrisations.** Formulate is just over half the run
  with MILU0 and no backtracking.
- [ ] **Reuse `dkrdh` in the Newton terms.** `calculate_coeffs_nwt` recomputes
  dk_r/dh numerically per connection with two `krelative` calls. `uft_cf`
  already stores `dkrdh` per cell.
- [ ] **Take the backtracking residual from the Newton system.** Verify that the
  residuals are equal, then drop the extra `inewton=0` build unless a backtrack
  happens.
- [ ] **Recommend or default to MILU0 for VSF models.** Confirm it on the 100 m
  model first.
- [ ] **Decide what happens to VSF_1D_ONLY.** Keep it as a development-mode
  option or drop it.

## Reproducing

Build MODFLOW 6 into `bin/` first (`pixi run build builddir`). Each run takes
2–8 minutes on one core:

```
cd autotest/vsf_benchmark
pixi run python vsf_big.py runs/d10m_ilut5 --delr 1000
pixi run python vsf_big.py runs/d10m_milu0 --delr 1000 --solver milu0
pixi run python vsf_big.py runs/d10m_milu0_nobt --delr 1000 --solver milu0_nobt
```

Omit `--delr` for the 100 m model. `--one-d` needs the VSF_1D_ONLY option from
the `vsf-1d-only` branch. Each run writes `stats.json` with step and
iteration counts and profiler timings. Do not run the `ilu0rcm` variant without
a time limit.
