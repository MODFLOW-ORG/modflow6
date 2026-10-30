# VS2D-compatible transpiration for MODFLOW 6 VSF

> **Status: planned later addition.** This was originally proposed as the
> first implementation. That decision was revised: **uncompensated Feddes is the
> day-one model** (see [feddes-transpiration-proposal.md](feddes-transpiration-proposal.md)),
> because it is substantially simpler to implement and carries higher immediate
> demand. The VS2D law remains planned — it plugs into the same VSE package and
> uptake-model abstraction as an additional model — and this document is the
> reference for that work. It is kept design-complete but is **not** the next
> thing to build.

This document reproduces the VS2D root-water uptake law for eventual
VS2D compatibility. The four phases below still
fit; verification cases should be defined alongside the physics, rather than
after implementation.

## 1. Physics: reproduce the VS2D uptake law

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

## 2. Architecture: fits the existing VSE package

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
> in the Feddes proposal, [§2a Distribution vs. reduction](feddes-transpiration-proposal.md#2a-distribution-vs-reduction-keep-column-logic-out-of-the-solver).
> That analysis is shared by both models; the column cap above is exactly the
> *column compensation* case identified there as one of the two reasons to pull
> stack-walking into MODFLOW.

## 3. Implementation: small, verifiable increments

1. Uptake calculation and reference parameters.
2. Register the model in the uptake-model factory; nonlinear assembly with the
   column cap, budgets, and observations (the VSE package already provides the
   wiring).
3. Time-varying demand and vegetation parameters.
4. Documentation and input-generation integration.

Soil evaporation remains a later, separate surface-boundary problem, as selected.

## 4. Verification: agree on these before coding

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
