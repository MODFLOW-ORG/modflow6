# NPF & STO flow-formulation refactor — migration guide for extension developers

Branch: `alt-formulations`

This document describes a breaking change to the NPF "flow formulation" and STO
"storage formulation" extension APIs in MODFLOW 6 core. If you maintain an
extension that implements a `GwfNpfFormulationType` (e.g. UZR flow, SWI flow, or
any custom flow term) **or** a `GwfStoFormulationType` (e.g. UZR/SWI storage),
your implementation **will not compile or behave correctly** against this branch
until it is updated. Feed this whole document to your AI coding assistant along
with your extension source to generate the required changes.

Both packages received the same conceptual change: dispatch moved from
**exclusive, one-formulation-per-element** to **additive, form-outer** (every
formulation runs over the whole grid and adds its terms; the standard term is now
itself a formulation that is always run first). Read sections 1–6 for NPF, then
section 7 for the STO specifics.

Four core files changed:

- `src/Model/ModelUtilities/GwfNpfExt.f90` — NPF extension contract (the abstract
  type `GwfNpfFormulationType` and its deferred/overridable methods).
- `src/Model/GroundWaterFlow/gwf-npf.f90` — the NPF package that drives the
  formulations.
- `src/Model/ModelUtilities/GwfStoExt.f90` — STO extension contract (the abstract
  type `GwfStoFormulationType`).
- `src/Model/GroundWaterFlow/gwf-sto.f90` — the STO package that drives the
  storage formulations.

---

## 1. What changed conceptually

### Before: exclusive, one-formulation-per-element dispatch

Each connection (and, for `cf`, each cell) was bound to exactly **one**
formulation via the `iformulation(:)` array (size `nja`). NPF owned the
cell/connection loop and, per element, dispatched to the single selected
formulation:

```fortran
! old npf_fc (and the same shape in npf_cf / npf_fn / npf_cq)
do n = 1, this%dis%nodes
  do ipos = ia(n)+1, ia(n+1)-1
    ...
    iform = this%iformulation(ipos)
    if (iform == DEFAULT_FLOW) then
      call this%fc_default_flow(n, m, ipos, ...)      ! hard-coded default
    else
      call this%flow_formulations(iform)%form%fc(n, m, ipos, ...)  ! one alternative
    end if
  end do
end do
```

Consequences of the old model:
- The formulation methods were **per-element** (`cf(kiter, n)`,
  `fc(n, m, ipos, ...)`), called from inside NPF's loop.
- A cell/face could have **only one** formulation. You could not have the
  standard conductance term on a face *plus* an additional term from another
  formulation on that same face.
- The default conductance path was **not** a formulation object; it was
  hard-coded as a special `DEFAULT_FLOW` branch.

### After: additive, form-outer dispatch

Now **every** formulation owns its own traversal of the grid and is called once
per phase. NPF simply loops over the active formulations and calls each; they
add their terms into the matrix/rhs/flowja additively:

```fortran
! new npf_fc (and the same shape in npf_cf / npf_fn / npf_cq)
call this%default_form%fc(kiter, matrix_sln, idxglo, rhs, hnew)
do iform = 1, MAX_EXT_FLOW_FORMS
  if (associated(this%flow_formulations(iform)%form)) then
    call this%flow_formulations(iform)%form%fc(kiter, matrix_sln, idxglo, rhs, hnew)
  end if
end do
```

Consequences of the new model:
- The formulation methods are now **whole-grid**: each method loops cells
  (`cf`) or connections (`fc`/`fn`/`cq`) **internally** and decides per element
  whether to contribute.
- Formulations **compose additively**. Multiple formulations can each add terms
  on the same cell/face.
- The standard conductance is now itself a formulation,
  `DefaultFlowFormulationType` (defined in `gwf-npf.f90`), and is always run
  first via `this%default_form`.
- Only `fc` is `deferred` (mandatory). `cf`, `fn`, and `cq` now have **no-op
  base implementations**, so a formulation overrides only the phases it
  contributes to.

### Claim mask (`iformulation`) semantics

`iformulation(:)` still exists (size `nja`, default value `DEFAULT_FLOW`) but its
meaning changed from "the one formulation for this face" to a **claim mask**: it
marks faces that an *exclusive* formulation owns. The default conductance
formulation **skips any face/cell whose `iformulation /= DEFAULT_FLOW`**:

```fortran
! inside DefaultFlowFormulationType routines
if (this%npf%iformulation(ipos) /= DEFAULT_FLOW) cycle
```

Use the claim mask only if your formulation **replaces** the standard
conductance on certain faces (set `iformulation(ipos)` to your form id so the
default skips them). If your formulation is purely **additive/supplemental**
(adds a term on top of the standard conductance), leave `iformulation` alone.

---

## 2. The new contract (`GwfNpfFormulationType`)

Full new definition in `src/Model/ModelUtilities/GwfNpfExt.f90`:

```fortran
type, abstract, public :: GwfNpfFormulationType
contains
  procedure(fc_if), deferred :: fc   ! mandatory
  procedure :: cf => cf_noop         ! override if you contribute to cf
  procedure :: fn => fn_noop         ! override if you contribute to newton
  procedure :: cq => cq_noop         ! override if you contribute to flows
end type GwfNpfFormulationType
```

### New method signatures (whole-grid)

```fortran
! Fill coefficients (DEFERRED — you MUST implement this)
subroutine fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
  class(<YourFormType>), intent(inout) :: this
  integer(I4B), intent(in) :: kiter
  class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
  integer(I4B), dimension(:), intent(in) :: idxglo
  real(DP), dimension(:), intent(inout) :: rhs
  real(DP), dimension(:), intent(inout) :: hnew

! Calculate coefficients (optional override; no-op by default)
subroutine cf(this, kiter)
  class(<YourFormType>), intent(inout) :: this
  integer(I4B), intent(in) :: kiter

! Fill newton terms (optional override; no-op by default)
subroutine fn(this, kiter, matrix_sln, idxglo, rhs, hnew)
  class(<YourFormType>), intent(inout) :: this
  integer(I4B), intent(in) :: kiter
  class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
  integer(I4B), dimension(:), intent(in) :: idxglo
  real(DP), dimension(:), intent(inout) :: rhs
  real(DP), dimension(:), intent(inout) :: hnew

! Calculate flows (optional override; no-op by default)
subroutine cq(this, hnew, flowja)
  class(<YourFormType>), intent(inout) :: this
  real(DP), dimension(:), intent(inout) :: hnew
  real(DP), dimension(:), intent(inout) :: flowja
```

### Signature change table (OLD → NEW)

| Phase | OLD (per-element)                                   | NEW (whole-grid, loops internally)                  | Status |
|-------|-----------------------------------------------------|-----------------------------------------------------|--------|
| `cf`  | `cf(this, kiter, n)`                                 | `cf(this, kiter)`                                   | no-op base; override optional |
| `fc`  | `fc(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)` | `fc(this, kiter, matrix_sln, idxglo, rhs, hnew)`   | **deferred; must implement** |
| `fn`  | `fn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)` | `fn(this, kiter, matrix_sln, idxglo, rhs, hnew)`   | no-op base; override optional |
| `cq`  | `cq(this, n, m, ipos, flowja, h_new)`               | `cq(this, hnew, flowja)`                            | no-op base; override optional |

Notes:
- `n`, `m`, `ipos` are **no longer passed in**. Your method now owns the loops
  and obtains them itself (see the loop template below).
- `fc` and `fn` gained `kiter`; lost `n, m, ipos`.
- `fc`/`fn` `hnew` is now `intent(inout)` (was `intent(in)` for `fc`).
- `cq` lost `n, m, ipos`; argument order is `(this, hnew, flowja)`.
- The old abstract interfaces `cf_if`, `fn_if`, `cq_if` were removed. Only
  `fc_if` remains as the single deferred interface.

---

## 3. What you must change in your extension

1. **Update every method signature** to match the table above.

2. **Move NPF's loop into each method.** Your method must iterate the grid
   itself. Use the model's connectivity through whatever back-pointer your
   formulation holds to the `GwfNpfType` / `DisBaseType` (the core default form
   keeps a `class(GwfNpfType), pointer :: npf` back-pointer — mirror that if you
   do not already have access to `dis`).

   Standard connection-loop template (for `fc`/`fn`/`cq`):

   ```fortran
   do n = 1, dis%nodes
     do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
       if (dis%con%mask(ipos) == 0) cycle            ! (cq historically does NOT mask)
       m = dis%con%ja(ipos)
       if (m < n) cycle                              ! upper triangle only
       ! ---- your applicability test here ----
       ! e.g. if (.not. applies_to_face(ipos)) cycle
       ! ---- add your terms for connection (n, m, ipos) ----
     end do
   end do
   ```

   Standard cell-loop template (for `cf`):

   ```fortran
   do n = 1, dis%nodes
     ! ---- your applicability test here ----
     ! ---- compute per-cell quantities ----
   end do
   ```

3. **Decide replace vs. additive:**
   - *Replace the standard conductance on a face:* set
     `npf%iformulation(ipos) = <your_form_id>` during setup so the default form
     skips that face, then supply the full term yourself in `fc`/`fn`/`cq`.
   - *Supplement (add on top):* leave `iformulation` as `DEFAULT_FLOW` and simply
     add your extra term; the default conductance still runs on that face.

4. **Only override the phases you need.** If your old formulation had a trivial
   `cf`, `fn`, or `cq`, delete the override entirely and inherit the no-op base.
   `fc` is mandatory (deferred).

5. **Registration is unchanged.** Keep calling:

   ```fortran
   call npf%add_flow_formulation(my_form_ptr, form_id)
   ```

   where `form_id` is in `1 .. MAX_EXT_FLOW_FORMS` and `my_form_ptr` is a
   `class(GwfNpfFormulationType), pointer`. The constants `DEFAULT_FLOW`,
   `UZR_FLOW`, `SWI_FLOW`, `MAX_EXT_FLOW_FORMS` are unchanged and still exported
   from `GwfNpfFormulationModule`.

---

## 4. Before/after example (illustrative)

### OLD extension (per-element, called from NPF's loop)

```fortran
subroutine my_fc(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
  class(MyFlowType), intent(inout) :: this
  integer(I4B), intent(in) :: n, m, ipos
  class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
  real(DP), dimension(:), intent(inout) :: rhs
  integer(I4B), dimension(:), intent(in) :: idxglo
  real(DP), dimension(:), intent(in) :: hnew
  ! ... add terms for this single connection (n, m, ipos) ...
end subroutine
```

### NEW extension (whole-grid, owns its loop)

```fortran
subroutine my_fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
  class(MyFlowType), intent(inout) :: this
  integer(I4B), intent(in) :: kiter
  class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
  integer(I4B), dimension(:), intent(in) :: idxglo
  real(DP), dimension(:), intent(inout) :: rhs
  real(DP), dimension(:), intent(inout) :: hnew
  ! local
  integer(I4B) :: n, m, ipos
  associate (dis => this%npf%dis)   ! however your form reaches dis
    do n = 1, dis%nodes
      do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
        if (dis%con%mask(ipos) == 0) cycle
        m = dis%con%ja(ipos)
        if (m < n) cycle
        if (.not. this%applies(ipos)) cycle
        ! ... add terms for connection (n, m, ipos) ...
      end do
    end do
  end associate
end subroutine
```

The type binding also changes — override only what you implement:

```fortran
type, extends(GwfNpfFormulationType) :: MyFlowType
  class(GwfNpfType), pointer :: npf => null()   ! add a back-pointer if needed
contains
  procedure :: fc => my_fc        ! required
  procedure :: fn => my_fn        ! only if you contribute newton terms
  procedure :: cq => my_cq        ! only if you contribute flows
  ! no cf binding needed if you don't implement cf (inherits no-op)
end type
```

---

## 5. How the core default formulation is wired (reference)

For a concrete, working example to mirror, see `DefaultFlowFormulationType` in
`src/Model/GroundWaterFlow/gwf-npf.f90`:

- Type definition (holds `class(GwfNpfType), pointer :: npf`) with
  `cf/fc/fn/cq` bound to `default_flow_cf/fc/fn/cq`.
- `default_flow_fc` / `default_flow_fn` / `default_flow_cq` implement the
  connection loop with the claim-mask skip
  (`if (this%npf%iformulation(ipos) /= DEFAULT_FLOW) cycle`).
- `default_flow_cf` implements the cell loop with the per-cell claim check at
  `idiag = this%npf%dis%con%ia(n)`.
- It is allocated in `allocate_arrays` (note the dummy argument is now
  `class(GwfNpftype), target :: this` so the back-pointer `form%npf => this` is
  valid) and deallocated in `npf_da`.

---

## 6. Behavior / correctness checklist (NPF)

- With no extension registered, output is bit-for-bit unchanged (the default
  form reproduces the original conductance, newton, and flow terms). Verified
  against the standard GWF regression tests (NPF, STO, TVK, newton paths).
- `xt3d` is still handled as an all-or-nothing early return inside
  `npf_cf/fc/fn/cq` and does not go through the formulation list. If your
  extension must coexist with `xt3d`, that integration point is unchanged and
  still needs separate handling.
- Guard against **double-counting**: if you replace the default on a face, make
  sure you set the claim mask so the default skips it; otherwise both terms are
  added.
- Keep `fc` and `fn` consistent per face (newton terms must match the residual
  terms your `fc` added).

---

## 7. STO (storage) package changes

The STO package was converted to the **same additive, form-outer model** as NPF.
The standard storage term is now a formulation, `DefaultStorageFormulationType`
(in `gwf-sto.f90`), always run first; registered extensions are then run and each
adds its terms.

### 7.1 What changed

- `sto_fc`, `sto_fn`, `sto_cq` no longer own the per-node loop and no longer
  select a single formulation via `iformulation(n)`. They now call the default
  storage formulation, then loop the registered external formulations
  (`sto_formulations(:)`, gated by the container's `is_active` logical) and call
  each. The steady-state (`iss`) / zero-`delt` guards and the `strgss`/`strgsy`
  reset remain in the `sto_*` drivers.
- The standard storage fill (`fc_default_sto` / `fn_default_sto` /
  `cq_default_sto`) is unchanged and is now invoked by
  `DefaultStorageFormulationType` via a back-pointer (`class(GwfStoType),
  pointer :: sto`).
- `iformulation(:)` (size `nodes`, default `DEFAULT_STORAGE`) is now a **claim
  mask**: the default storage formulation skips any node with
  `iformulation(n) /= DEFAULT_STORAGE`.

### 7.2 Contract changes in `GwfStoExtModule` (`GwfStoFormulationType`)

Only the three looping methods changed signature (they now own the node loop and
lost the `n` argument; `fc`/`fn` gained `kiter`). All methods remain **deferred**
(STO did not adopt no-op base methods — a registered storage formulation must
implement all six, exactly as before).

| Method       | OLD (per-element)                                      | NEW (whole-grid, loops nodes itself)                 | Changed? |
|--------------|--------------------------------------------------------|------------------------------------------------------|----------|
| `is_active`  | `is_active(this, n) -> logical`                        | *(unchanged)*                                        | no |
| `fc`         | `fc(this, n, matrix_sln, rhs, idxglo, h_old, h_new)`   | `fc(this, kiter, matrix_sln, rhs, idxglo, h_old, h_new)` | **yes** |
| `fn`         | `fn(this, n, matrix_sln, rhs, idxglo, h_old, h_new)`   | `fn(this, kiter, matrix_sln, rhs, idxglo, h_old, h_new)` | **yes** |
| `cq`         | `cq(this, n, flowja, h_new, h_old)`                    | `cq(this, flowja, h_new, h_old)`                     | **yes** |
| `bd`         | `bd(this, isuppress_output, model_budget)`             | *(unchanged)*                                        | no |
| `save_flows` | `save_flows(this, iprint, ibinun)`                     | *(unchanged)*                                        | no |

Important: gfortran requires the **overriding procedure's dummy-argument names to
match the interface exactly**. Use `h_old` and `h_new` (not `hold`/`hnew`) in
your `fc`/`fn`/`cq` implementations, or the override will be rejected and your
type will be treated as still-abstract.

### 7.3 What you must change in your STO extension

1. Update the `fc`, `fn`, `cq` signatures per the table (drop `n`, add `kiter` to
   `fc`/`fn`, keep `h_old`/`h_new` dummy names).

2. Move the node loop into each of `fc`, `fn`, `cq`, and use your `is_active(n)`
   as the per-node applicability test:

   ```fortran
   subroutine my_sto_fc(this, kiter, matrix_sln, rhs, idxglo, h_old, h_new)
     class(MyStoType), intent(inout) :: this
     integer(I4B), intent(in) :: kiter
     class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
     real(DP), dimension(:), intent(inout) :: rhs
     integer(I4B), dimension(:), intent(in) :: idxglo
     real(DP), dimension(:), intent(in) :: h_old
     real(DP), dimension(:), intent(in) :: h_new
     integer(I4B) :: n
     do n = 1, this%sto%dis%nodes          ! however your form reaches dis
       if (this%sto%ibound(n) <= 0) cycle
       if (.not. this%is_active(n)) cycle
       ! ... add your storage terms for node n ...
     end do
   end subroutine
   ```

3. Decide replace vs. additive (same rule as NPF): to **replace** standard
   storage on a node, set `sto%iformulation(n) = <your_form_id>` during setup so
   the default skips it; to **supplement**, leave it `DEFAULT_STORAGE`.

4. `bd` and `save_flows` are unchanged and are still called form-outer for every
   registered formulation (gated by the container `is_active` logical set in
   `add_sto_formulation`). Keep them as-is.

5. Registration is unchanged: `call sto%add_sto_formulation(my_form_ptr,
   form_id)` with `form_id` in `1 .. MAX_EXT_STO_FORMS`. The constants
   `DEFAULT_STORAGE`, `UZR_STORAGE`, `SWI_STORAGE`, `MAX_EXT_STO_FORMS` are
   unchanged.

### 7.4 Reference: core default storage formulation

See `DefaultStorageFormulationType` in `src/Model/GroundWaterFlow/gwf-sto.f90`:

- Holds `class(GwfStoType), pointer :: sto`; binds `is_active/fc/fn/cq/bd/
  save_flows` to `default_storage_*`.
- `default_storage_fc/fn/cq` loop all nodes, apply the `ibound` and claim-mask
  checks, and delegate to the existing `fc_default_sto/fn_default_sto/
  cq_default_sto`.
- `default_storage_is_active` returns `.false.` (the default uses the claim mask,
  not `is_active`); `default_storage_bd` and `default_storage_save_flows` are
  no-ops (the default storage budget/flows are handled directly by `sto_bd` /
  `sto_save_model_flows` through `strgss`/`strgsy`).
- Allocated in `allocate_arrays` (dummy is `class(GwfStoType), target :: this`)
  and deallocated in `sto_da`.

### 7.5 STO correctness checklist

- With no extension registered, STO output is unchanged. Verified against the
  STO regression tests (including newton and TVS paths).
- Guard against double-counting: if you replace standard storage on a node, set
  the claim mask so the default skips it.
- Keep `fc`/`fn` consistent per node (newton terms must match the residual terms
  your `fc` added).
- The `strgss`/`strgsy` arrays are reset once per `sto_cq` call before any
  formulation runs; your `cq` should only *add* to its own storage-rate arrays
  (external formulations maintain their own rate arrays and report them via `bd`
  / `save_flows`).
