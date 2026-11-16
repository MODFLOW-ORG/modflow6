# Conductivity Provider

An abstract interface that gives boundary (and other) packages a cell's full
**hydraulic conductivity tensor** without depending on NPF internals.

## Why

Packages such as a seepage-face boundary need the hydraulic conductivity of a
cell to resolve flow across a face. In fully 3D anisotropic media the Darcy flux
is a tensor law `q = -K * grad(h)`, and the normal flux across a face,
`q_n = -n . K . grad(h)`, couples the face normal to the *full* gradient. A
single directional scalar (an "effective K" along one direction) cannot
represent that coupling when the gradient is not aligned with the face normal.
So the provider hands back the full tensor and lets the consumer do the physics.

That conductivity data (`k11/k22/k33` + rotation `angle1/2/3`) lives in NPF.
Reaching into NPF memory from each consumer duplicates physics and couples
packages to NPF internals. This interface centralizes the capability and is
injected by the model.

## Design

- `ConductivityProviderType` (abstract, `src/Model/ModelUtilities/ConductivityProvider.f90`)
  with deferred:
  - `conductivity(n) -> K(3,3)`: **saturated** hydraulic conductivity tensor of
    reduced node `n`, in model coordinates (anisotropy and rotation resolved).
    The consumer owns any projection or conductance assembly.
  - `krel(n) -> kr`: current relative permeability (1 when saturated, < 1 where
    an unsaturated formulation is active).
  - `dkrel_dh(n) -> dkrdh`: head derivative of `krel` (0 when saturated), for
    Newton-Raphson boundary terms.
- `GwfNpfType%calc_conductivity(n)` (`gwf-npf.f90`)
  - builds `K = R * diag(k11, k22, k33) * R^T`, where `R` is the rotation
    implied by `angle1/2/3` (same construction XT3D uses for its conductivity
    tensor). Returns the full 3x3 symmetric tensor.
  - `GwfNpfType%krel` / `%dkrdh` arrays hold the relative permeability and its
    head derivative. They default to 1 and 0 and stay inert until an unsaturated
    flow formulation (e.g. UZR) populates them, so saturated models are
    unaffected.
- `NpfConductivityProviderType` (`src/Model/ModelUtilities/NpfConductivityProvider.f90`)
  - concrete adapter holding a `GwfNpfType` pointer, delegating `conductivity`
    to `calc_conductivity` and `krel`/`dkrel_dh` to the NPF arrays. Needed
    because `GwfNpfType` already extends `NumericalPackageType` and Fortran is
    single-inheritance, so NPF cannot also extend the abstract provider directly.
- `GwfModelType%conductivity_provider` (`gwf.f90`)
  - created in `gwf_ar` (`create_npf_conductivity_provider`), deallocated in
    `gwf_da`.
- Injection hook on `BndType` (`BoundaryPackage.f90`)
  - `conductivity_provider` pointer + `set_conductivity_provider`; `gwf_ar`
    calls `packobj%set_conductivity_provider(this%conductivity_provider)` for
    every boundary package. The base implementation only stores the (possibly
    null) pointer, so there is no behavior change until a package consumes it.

```
consumer -> ConductivityProviderType%{conductivity,krel,dkrel_dh} -> NpfConductivityProvider -> GwfNpfType
```

A consumer assembles the Darcy flux it needs from these pieces, e.g. the normal
flux across a seepage face of normal `n` with gradient `grad_h`:

```
q_n = -krel(n) * dot_product(nrm, matmul(K, grad_h))
```

with `dkrel_dh` feeding the Newton terms (the tensor and face normal are
head-independent, so differentiation goes through `krel` and `grad_h`).

## Scope of this branch

Infrastructure only: the provider is created, held by the model, and injected
into every boundary package, but not yet consumed, so there is **no behavior
change** (all existing tests pass, library + `mf6` build clean). `krel` is 1 and
`dkrel_dh` is 0 everywhere until an unsaturated formulation populates the NPF
arrays. The first real consumer (a seepage-face boundary, initially for Richards
flow) lands on the UZR branch.

## Testing

`autotest/TestConductivityProvider.f90` (test-drive) covers abstract dispatch of
all three deferred procedures via a mock (including the full 3x3 tensor) and the
directional `hyeff` values that a consumer projection should reproduce. It is
registered but only built/run when the `test-drive` dependency is available.

## Future

Additional providers (e.g. non-NPF flow formulations) can implement
`ConductivityProviderType` and be injected the same way. BUY/VSC, which
currently pull NPF memory directly, could migrate later.
