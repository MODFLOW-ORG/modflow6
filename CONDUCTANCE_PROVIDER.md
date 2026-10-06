# Conductance Provider

An abstract interface that gives boundary (and other) packages a directional,
anisotropy-resolved hydraulic conductivity without depending on NPF internals.

## Why

Packages such as the seepage-face boundary need the effective K of a cell in a
specific direction (a face normal), which can be anisotropic. That logic and the
conductivity data live in NPF. Reaching into NPF memory from each consumer
duplicates physics and couples packages to NPF internals. This interface
centralizes the capability and is injected by the model.

## Design

- `ConductanceProviderType` (abstract, `src/Model/ModelUtilities/ConductanceProvider.f90`)
  - deferred `eff_hy(n, ihc, vg) -> hy`: effective K of reduced node `n` along
    unit vector `vg`; `ihc` selects the vertical/horizontal anisotropy branch.
- `GwfNpfType%calc_eff_hy(n, ihc, vg)` (`gwf-npf.f90`)
  - connection-agnostic public wrapper over the existing `hy_eff` (which already
    accepts an explicit direction and ignores the neighbor/ipos args).
- `NpfConductanceProviderType` (`src/Model/ModelUtilities/NpfConductanceProvider.f90`)
  - concrete adapter holding a `GwfNpfType` pointer, delegating `eff_hy` to
    `calc_eff_hy`. Needed because `GwfNpfType` already extends
    `NumericalPackageType` and Fortran is single-inheritance, so NPF cannot also
    extend the abstract provider directly.
- `GwfModelType%cond_provider` (`gwf.f90`)
  - created in `gwf_ar` (`create_npf_conductance_provider`), deallocated in
    `gwf_da`. Available for packages to consume.

```
consumer ──▶ ConductanceProviderType%eff_hy ──▶ NpfConductanceProvider ──▶ GwfNpfType%calc_eff_hy ──▶ hy_eff ──▶ hyeff
```

## Scope of this branch

Infrastructure only: the provider is created and held by the model but not yet
consumed, so there is **no behavior change** (all existing tests pass, library +
`mf6` build clean). Per-package injection (e.g. a `select type` handing the
provider to a boundary package in the `gwf_ar` loop) lands on the consuming
branch, e.g. the seepage-face (SPF) work.

## Testing

`autotest/TestConductanceProvider.f90` (test-drive) covers abstract dispatch via
a mock and the directional `hyeff` values. It is registered but only built/run
when the `test-drive` dependency is available; running it is deferred.

## Future

Additional providers (e.g. XT3D-specialized, or non-NPF flow formulations) can
implement `ConductanceProviderType` and be injected the same way. BUY/VSC, which
currently pull NPF memory directly, could migrate later.
