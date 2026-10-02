# Proposal: Centralize directional effective‑K and inject it into boundary packages

> Related: this proposal addresses action item #1 in
> [SPF_SEEPAGE_ACTION_ITEMS.md](SPF_SEEPAGE_ACTION_ITEMS.md) (the cached NPF
> pointers / non‑dispatched `bnd_ar` workaround).

> **Decision (2026‑10‑02):** going with the long‑term **Phase 2** (abstract
> provider) straight away. The provider infrastructure is being developed on a
> clean branch off `develop` — **`conductance-provider`** (worktree at
> `../modflow6-conductance-provider`) — to be merged first, independently of SPF.
> See that branch's `CONDUCTANCE_PROVIDER.md` for the as‑built design. The SPF
> package will consume it (replacing `effective_k` + the cached NPF pointers with
> `this%cond_provider%eff_hy(...)`) on this branch afterward.

## 1. Problem

`SpfType%effective_k` (in `src/Model/GroundWaterFlow/gwf-spf.f90`) computes the
anisotropy‑resolved hydraulic conductivity of a cell in the direction of a
seepage‑face normal. To do this it:

- **duplicates** the branch/flag logic that already lives in
  `GwfNpfType%hy_eff` (`src/Model/GroundWaterFlow/gwf-npf.f90`), and
- **reaches into NPF‑owned memory** (`K11`, `K22`, `K33`, `ANGLE1..3`, `IK22`,
  `IANGLE1..3`, `IAVGKEFF`) via `mem_setptr`, caching eleven pointers in SPF.

This is fragile: the physics can drift from NPF, and SPF knows far too much about
NPF internals. We want the logic in one central place, reusable by other packages
(current and future), with NPF remaining the owner of the conductivity data.

## 2. Goals / non‑goals

**Goals**
- Single source of truth for directional effective‑K (no duplicated flag logic).
- Boundary packages hold **one** handle, not a bag of NPF arrays/scalars.
- Reusable by future packages that need a face/segment conductance from the grid.
- No new circular module dependencies.

**Non‑goals (for now)**
- Changing the physics or the `hyeff` math itself.
- Supporting non‑NPF flow providers (kept as a later, optional extension).
- Reworking BUY/VSC, which already pull NPF memory by the same precedent.

## 3. What already exists

- `HGeoUtilModule%hyeff(k11,k22,k33,ang1,ang2,ang3, vg1,vg2,vg3, iavgmeth)` —
  the **pure, stateless** anisotropic‑K math. Already central. Keep as‑is.
- `GwfNpfType%hy_eff(n, m, ihc, ipos, vg)` — the **stateful wrapper** that knows
  which cell arrays/flags to read and calls `hyeff`. Public. Already does exactly
  what SPF needs; it even accepts an explicit `vg` (direction vector) and then
  ignores `m`/`ipos`.

So the physics is *already* centralized in NPF. The real gap is purely a
**plumbing** one: SPF has no handle to the NPF object, so it re‑implemented the
logic against raw memory instead.

## 4. Proposed design

### Part A — Expose a connection‑agnostic API on NPF

Add a thin public method that computes effective‑K for a cell in an arbitrary
direction, without pretending there is a neighbor connection:

```fortran
!> Effective hydraulic conductivity of cell n along unit direction vg.
!> ihc selects the vertical (0) or horizontal (1) anisotropy branch.
function calc_eff_hy(this, n, ihc, vg) result(hy)
  class(GwfNpfType) :: this
  integer(I4B), intent(in) :: n
  integer(I4B), intent(in) :: ihc
  real(DP), dimension(3), intent(in) :: vg
  real(DP) :: hy
  hy = this%hy_eff(n, n, ihc, vg=vg)   ! m/ipos unused when vg is supplied
end function
```

Optionally make `m`/`ipos` truly optional in `hy_eff` so the dummy `m = n`
disappears. A small companion helper builds the normal so callers don't repeat
trig:

```fortran
!> Unit normal of a face from its connection type and x-azimuth (degrees).
pure function face_normal(ihc, angldegx) result(vg)   ! HGeoUtil or NPF
  ...
  if (ihc == 0) then
    vg = [DZERO, DZERO, DONE]
  else
    vg = [cos(angldegx*DPIO180), sin(angldegx*DPIO180), DZERO]
  end if
end function
```

Net effect: SPF's `effective_k` collapses to one call
`this%npf%calc_eff_hy(n, ihc, face_normal(ihc, angldegx))`, and all the
`K*/ANGLE*/I*` pointers disappear from SPF.

### Part B — Inject the NPF object into the package (the pattern)

The blocker is that SPF is a generic `BndType` in the model's `bndlist`, created
by the package factory, with no route to the NPF object. Two patterns, phased:

#### Phase 1 (recommended now): model‑wired object pointer

Give SPF a single typed pointer and let the **model** wire it in `gwf_ar`, right
where it already special‑cases BUY/VSC per package:

```fortran
! gwf-spf.f90
use GwfNpfModule, only: GwfNpfType
type(GwfNpfType), pointer :: npf => null()
```

```fortran
! gwf.f90, inside the bnd loop after npf_ar (which runs first)
do ip = 1, this%bndlist%Count()
  packobj => GetBndFromList(this%bndlist, ip)
  call packobj%set_pointers(...)
  call packobj%bnd_ar()
  if (this%innpf > 0) then
    select type (packobj)
    type is (SpfType)
      packobj%npf => this%npf
    end select
  end if
  ...
end do
```

**Why model‑level wiring (not an NPF method that sets the package):** SPF already
`use`s `GwfNpfModule`. If NPF `use`d `SpfModule` to set the pointer, we would get a
**circular module dependency** (`SpfModule → GwfNpfModule → SpfModule`), which
Fortran cannot compile. The model (`gwf.f90`) already imports NPF and can import
SPF without a cycle, so it is the correct place for the `select type`. This mirrors
the existing `vsc_ar_bnd` `select type (packobj)` activation pattern.

This also sidesteps the known issue where SPF's `bnd_ar` override was never
dispatched — the model loop is guaranteed to run.

#### Phase 2 (optional, later): abstract provider interface

To fully decouple and allow non‑NPF providers (XT3D specializations, future flow
formulations, unit tests with a stub), define a small abstract type in a base
module that has **no** NPF dependency:

```fortran
! new low-level module, e.g. ConductanceProviderModule
type, abstract :: ConductanceProviderType
contains
  procedure(eff_hy_if), deferred :: eff_hy   ! (n, ihc, vg) -> K
end type
```

`GwfNpfType` `extends`/implements it; `BndType` can then hold a
`class(ConductanceProviderType), pointer` (generic, no NPF import), set by the
model. Packages call `this%provider%eff_hy(n, ihc, vg)` with zero knowledge of
NPF. This is the clean long‑term home and makes the capability trivially reusable
by any future package, but it is more machinery than needed to fix SPF today.

## 5. Migration steps (Phase 1)

1. NPF: add `calc_eff_hy` (and optionally make `hy_eff`'s `m`/`ipos` optional);
   add `face_normal` helper (NPF or `HGeoUtil`).
2. SPF: replace the eleven NPF pointers + `effective_k` + `set_npf_pointers` with
   a single `type(GwfNpfType), pointer :: npf` and a one‑line `calc_eff_hy` call.
3. `gwf.f90`: add the `select type (packobj)` wiring in the `gwf_ar` bnd loop
   (import `SpfType`).
4. Rebuild; `autotest/test_gwf_uzr_spf.py` must stay green (same numbers —
   `calc_eff_hy` is the same math as the current `effective_k`).

## 6. Tradeoffs

| Aspect | Phase 1 (object injection) | Phase 2 (abstract provider) |
|---|---|---|
| Effort | Small | Medium |
| Decoupling | SPF depends on concrete NPF type | SPF depends only on an interface |
| Reuse by other pkgs | Yes (same wiring per pkg) | Yes (cleanest) |
| Non‑NPF providers | No | Yes |
| Risk | Low | Low–medium (new abstraction) |

**Recommendation:** Do Phase 1 now (removes the duplication and the pointer soup
with minimal risk). Revisit Phase 2 when a second consumer or a non‑NPF provider
actually appears.

## 7. Notes / open questions

- `hy_eff` currently takes `(n, m, ihc, ipos, vg)`. Decide whether to add a clean
  `calc_eff_hy` wrapper (least churn) or refactor `hy_eff`'s signature.
- Where should `face_normal` live — `HGeoUtil` (pure geometry) or NPF? `HGeoUtil`
  keeps it reusable and dependency‑free.
- BUY and VSC use the same "pull NPF memory" approach; if Phase 2 lands, they
  could migrate to the provider too (out of scope here).
```
