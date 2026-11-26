# Plan: Memory-Manager Metadata for API-Modifiable Variables

## Motivation

The MODFLOW 6 API (BMI/XMI, via `libmf6`, `xmipy`, and `modflowapi`) exposes
**every** variable in the memory manager, with no indication of which variables a
user may meaningfully change. Today `get_input_var_names` and
`get_output_var_names` return the same list — "all of them"
(`srcbmi/mf6bmi.f90`). This is a foot-gun:

- Some variables are **derived** once and then used directly. Setting the raw
  input has no effect. The canonical example is intercell conductance: NPF
  computes `CONDSAT` from `K`/`K22`/`K33` during Allocate and Read, then uses it
  directly. To change conductance through the API a user must set `CONDSAT`, not
  `K` — which is surprising and undocumented.
- Some variables are **live input** read every formulate step (`STO/SS`,
  `STO/SY`) and can be set directly.
- Some are **period-reset** (re-read each stress period), **overwritten**
  (recomputed every step), **state**, or **structural/read-only**.

This plan proposes (1) adding metadata to the memory manager so each variable's
API role is explicit, queryable, and documentable; and (2) a mechanism so that
setting an *input* such as `K` can trigger recomputation of its derived
quantity (`CONDSAT`), removing the foot-gun. The companion user document is the
`doc/mf6api/` guide (the *Description of Variables Modifiable Through the API*).

## Goals

1. Make each memory-manager variable's API modifiability **explicit** in the
   code where it is allocated.
2. **Expose** that classification through the BMI so external tools (and
   `get_input_var_names`) are meaningful.
3. **Auto-generate** the per-package tables in the `doc/mf6api/` guide from the
   metadata, so docs cannot drift from the code.
4. Provide a **recompute-on-set** path for key derived variables, starting with
   `NPF/CONDSAT` from `K`, reusing existing infrastructure.
5. Enforce/guide at run time: optionally warn (or error) when a variable that is
   classified as having no effect (or read-only) is set through the API.

## Part 1 — Variable classification metadata

### 1.1 Taxonomy (mirrors the `doc/mf6api/` guide)

| Class | Meaning |
|---|---|
| `LIVE_INPUT` | read directly each formulate step; set freely (e.g. `STO/SS`) |
| `DERIVED` | computed from input and then used directly; set this to change behavior (e.g. `NPF/CONDSAT`) |
| `INERT_INPUT` | only used to compute a derived variable; setting has no effect unless a recompute is triggered (e.g. `NPF/K`) |
| `PERIOD_RESET` | re-read from input each stress period (e.g. `WEL/BOUND`, `GHB/COND`) |
| `OVERWRITTEN` | recomputed every step/iteration; setting is futile (e.g. `HCOF`, `RHS`, `SIMVALS`) |
| `STATE` | solver/model state; settable with care (e.g. `X`) |
| `STRUCTURAL` | grid/connectivity/sizes; read-only (e.g. `DIS/NODES`, `MAXBOUND`) |
| `UNSET` | default; not yet classified |

Note that a variable's class can be **option-dependent** (e.g. `NPF/CONDSAT` is
`DERIVED` normally, `OVERWRITTEN` when TVK/VSC is active, and absent under XT3D).
The metadata should be set to reflect the configured model where practical, or
documented as conditional.

### 1.2 Add an attribute to `MemoryType`

`MemoryType` (`src/Utilities/Memory/Memory.F90`) currently carries `name`,
`mastername`, `memtype`, `isize`, `element_size`, and `set_handler_idx`. Add:

```fortran
integer(I4B) :: api_class = MEMAPI_UNSET   !< API modifiability class
```

Define the class enum (and parameter names) in `ConstantsModule` or a small new
`MemoryApiClassModule`.

### 1.3 Set the class at allocation

Two options:

- **(a) New optional argument to `mem_allocate`** — `api_class=` defaulting to
  `UNSET`. Packages pass the class as they allocate. Most explicit, but touches
  every `mem_allocate` call we want to classify.
- **(b) A post-allocation tagging pass** — a helper
  `mem_set_api_class(var, path, class)` called once (e.g. in each package's
  `*_ar`/allocate routine) for the variables we want to classify. Less invasive;
  classification lives in one place per package.

Recommendation: start with **(b)** for the pilot packages (NPF, STO, the list
stress packages), since it is incremental and keeps the classification readable
and grouped. Migrate toward **(a)** if it proves cleaner.

### 1.4 Expose through the BMI

- Add a function (e.g. `get_var_api_class(c_var_address)`) returning the class
  string, mirroring `get_var_type`/`get_var_rank` in `srcbmi/mf6bmi.f90`.
- Make `get_input_var_names` / `get_output_var_names` filter by class
  (input-like vs output-like) instead of returning everything.

### 1.5 Auto-generate the `doc/mf6api/` tables

Add a small generator (analogous to `mf6ivar`) that runs a model (or introspects
the registered variables) and emits the per-package modifiability tables from the
`api_class` metadata, replacing the hand-written tables in `doc/mf6api/`. This
keeps the guide synchronized with the code and is the long-term payoff of the
metadata.

### 1.6 Caveat: the metadata is advisory, not enforceable

The classification cannot, by itself, *prevent* a user from making an
ineffective or unsafe change. A value written through a `get_value_ptr` pointer
goes straight to the variable's memory with no interception point (the same
reason a set-handler does not fire for pointer writes — see Part 2, Approach B).
There is therefore no hook at which the program could reject, warn about, or
react to such a write. The metadata is consequently **advisory**, useful for:

- **documentation** — driving the `doc/mf6api/` tables (§1.5);
- **discovery** — filtering `get_input_var_names` / `get_output_var_names` and
  answering `get_var_api_class` queries (§1.4);
- **best-effort run-time feedback** — warning (or erroring) when an
  `OVERWRITTEN`/`STRUCTURAL` variable is set, but **only** on the `set_value`
  path, since that is the only setter the program can observe.

It is **not** a hard guard: a user writing through a pointer can still set
anything, including read-only/structural variables, with no diagnostic. The
recompute mechanisms in Part 2 (especially the setter-agnostic trigger of
Approach A) are what actually keep derived quantities consistent; the metadata
informs the user, it does not enforce.

## Part 2 — Recompute-on-set for derived variables (`K` → `CONDSAT`)

### 2.1 The problem

`NPF/CONDSAT` is `DERIVED`: setting `K` through the API after Allocate and Read
has no effect, because the conductance was already computed. Users must set
`CONDSAT` directly, which requires understanding the symmetric connection
indexing and bypasses any consistency checks.

### 2.2 Existing infrastructure to reuse

Both pieces already exist:

1. **Recompute routine.** `NPF` has `calc_condsat(node, upperOnly)`
   (`src/Model/GroundWaterFlow/gwf-npf.f90`). For **time-varying K (TVK)**,
   `npf_ad` already recomputes `CONDSAT` for changed nodes:

   ```fortran
   ! npf_ad: if any K values changed this step, update CONDSAT (or XT3D coeffs)
   if (this%kchangeper == kper .and. this%kchangestp == kstp) then
     if (this%ixt3d == 0) then
       do n = 1, this%dis%nodes
         if (this%nodekchange(n) == 1) call this%calc_condsat(n, .false.)
       end do
     else
       ... recompute XT3D coefficients ...
     end if
   end if
   ```

   The flags `KCHANGEPER`, `KCHANGESTP`, and the `NODEKCHANGE(nodes)` array are
   the existing signalling mechanism; the VSC correction and the XT3D path are
   already handled here.

2. **Set-side-effect handler.** `src/Utilities/Memory/MemorySetHandler.f90`
   provides `mem_register_handler(var_name, mem_path, handler, ctx)` and
   `on_memory_set(var_name, mem_path, status)`. The BMI's `set_value_double`/
   `set_value_int` already call `on_memory_set` after copying a value
   (`srcbmi/mf6bmi.f90`). **No package registers a handler today**, so this
   infrastructure is dormant but ready.

### 2.3 Two implementation approaches

#### Approach A — explicit "recompute" trigger variable (discussed 2026-06-15)

Add a small integer flag to NPF, e.g. `IRECALCCONDSAT` (memory address
`MODEL/NPF/IRECALCCONDSAT`), default `0`. A user workflow is:

1. set `MODEL/NPF/K` (and `K22`/`K33`) to new values through the API;
2. set `MODEL/NPF/IRECALCCONDSAT = 1`.

NPF checks the flag at the start of the next formulate (or in `npf_ad`/a new
`npf_cf` guard): if `1`, recompute `CONDSAT` for all nodes (reusing the
`calc_condsat` loop, honoring the XT3D/VSC branches), then **reset the flag to
0**. This mirrors the user's description ("recalculate satcond when set to 1 and
then reset to 0").

Implementation notes:
- Reuse the existing recompute block. The simplest correct version sets
  `NODEKCHANGE(:) = 1`, `KCHANGEPER = kper`, `KCHANGESTP = kstp`, then runs the
  same loop the TVK path runs — i.e., funnel the trigger into the existing
  machinery rather than duplicating it.
- Pros: explicit, discoverable, easy to document and test; one trigger covers a
  batch of `K` edits; no dependence on which setter (`set_value` vs a pointer
  write) was used.
- Cons: a new public input-context/memory variable; two-step user workflow.

#### Approach B — register a set-handler on `K`/`K22`/`K33`

In `npf_ar`, call `mem_register_handler` for `K11`/`K22`/`K33` with a handler
that recomputes `CONDSAT`. Then a user simply does `set_value("MODEL/NPF/K", ...)`
and the recompute happens automatically (the BMI's `set_value` calls
`on_memory_set`).

> **Important limitation.** The handler fires **only for `set_value`**. A write
> through a pointer obtained from `get_value_ptr` goes directly to the variable's
> memory with no interception point, so `on_memory_set` is never called and the
> handler does **not** run. This matters because `get_value_ptr` + in-place write
> is the efficient, recommended access pattern (it is what the `doc/mf6api/`
> guide and the `modflowapi` examples use). A user following that pattern would
> set `K` through a pointer and see **no** `CONDSAT` update — silently wrong
> results. There is no way to intercept a raw pointer write, so a set-handler
> cannot be the *sole* mechanism; an explicit, setter-agnostic trigger
> (Approach A) is required for pointer-based workflows.

Implementation notes:
- Pros: no new variable; "set `K` and it just works"; uses the
  infrastructure that was clearly designed for this.
- Cons: only fires for `set_value` (a write through a `get_value_ptr` pointer
  bypasses `on_memory_set` — see the `doc/mf6api/` "pointers vs copies"
  section); recomputing on *every* `set_value` of `K` is wasteful if a user sets
  `K`, `K22`, `K33` in sequence (three recomputes); ordering/coalescing needs
  thought.

##### Example Fortran

The handler must match the `set_handler_iface` abstract interface in
`MemorySetHandlerModule` — `subroutine(owner, status)`, where `owner` is the
polymorphic context that was registered. Register the conductivity arrays in
`npf_ar` (after the arrays are allocated and `CONDSAT` first computed), passing
the NPF object itself as the context:

```fortran
subroutine npf_ar(this, ic, vsc, ibound, hnew)
  use MemorySetHandlerModule, only: mem_register_handler
  ! the passed-object must be a target so the context pointer can point at it;
  ! the NPF object lives for the whole simulation, so the pointer stays valid
  class(GwfNpfType), target :: this
  ! ... existing dummy arguments and body ...
  procedure(set_handler_iface), pointer :: k_handler
  class(*), pointer :: ctx

  ! ... existing Allocate and Read work, including the first CONDSAT calc ...

  ! register a recompute-on-set handler for the conductivity arrays, so that
  ! setting K (K22/K33) through the BMI set_value triggers a CONDSAT update
  k_handler => npf_recompute_condsat
  ctx => this
  call mem_register_handler('K11', this%memoryPath, k_handler, ctx)
  if (this%ik22 /= 0) &
    call mem_register_handler('K22', this%memoryPath, k_handler, ctx)
  if (this%ik33 /= 0) &
    call mem_register_handler('K33', this%memoryPath, k_handler, ctx)
end subroutine npf_ar
```

The handler (a module procedure) recovers the NPF object from the polymorphic
`owner` and funnels into the same recompute path that TVK uses in `npf_ad`:

```fortran
!> @brief Set-handler: recompute CONDSAT after K is set through the API.
!!
!! Matches set_handler_iface. The memory manager calls this (via on_memory_set,
!! which the BMI set_value invokes) after K/K22/K33 is set, passing back the NPF
!! object as the polymorphic owner/context registered in npf_ar.
!<
subroutine npf_recompute_condsat(owner, status)
  use TdisModule, only: kper, kstp
  class(*), pointer, intent(inout) :: owner  !< registered context (the NPF object)
  integer(I4B), intent(out) :: status        !< 0 = success, /= 0 = error
  integer(I4B) :: n

  status = 0
  select type (owner)
  class is (GwfNpfType)
    ! Debounce option: instead of recomputing here, set a flag (the
    ! IRECALCCONDSAT of Approach A) and let the next formulate do one recompute,
    ! so sequential K11/K22/K33 sets coalesce. Direct version shown here:
    owner%nodekchange(:) = 1
    owner%kchangeper = kper
    owner%kchangestp = kstp
    if (owner%ixt3d == 0) then
      do n = 1, owner%dis%nodes
        call owner%calc_condsat(n, .false.)
      end do
    else
      ! recompute XT3D coefficients for permanently confined connections,
      ! exactly as the TVK branch in npf_ad does
      if (owner%xt3d%lamatsaved .and. .not. owner%xt3d%ldispersion) then
        call owner%xt3d%xt3d_fcpc(owner%dis%nodes, .true.)
      end if
    end if
  class default
    status = -1  ! unexpected context registered for this handler
  end select
end subroutine npf_recompute_condsat
```

Notes:
- This reuses `calc_condsat` and the `NODEKCHANGE`/`KCHANGEPER`/`KCHANGESTP`
  machinery, so behavior matches the existing TVK per-step recompute (including
  the VSC-corrected K — apply the VSC correction first if VSC is active, as
  `npf_ad` does).
- The `class default` path returns a non-zero `status`; `on_memory_set`
  propagates it and the BMI reports the error, since the memory would otherwise
  be left inconsistent.
- For the debounce concern above, the handler can simply set
  `owner%irecalccondsat = 1` (Approach A's flag) and let the next formulate run
  one recompute — combining A and B with a single recompute path.

#### Recommendation

Implement **Approach A** (the trigger flag) as the primary mechanism. It is
explicit, batch-friendly, cleanest to document and test, and — decisively —
**setter-agnostic**: it works whether the user wrote `K` through `set_value` or
through a `get_value_ptr` pointer. Because the recommended, efficient access
pattern is an in-place pointer write, which **cannot** trigger a set-handler (see
the limitation under Approach B), Approach A is the only option that works for
the common workflow; a set-handler alone would silently fail to update `CONDSAT`
for pointer-based users.

Approach B can be added later as a convenience for `set_value` users, ideally
debounced so its handler just sets Approach A's `recompute` flag rather than
recomputing immediately (so multiple `K`/`K22`/`K33` edits coalesce into one
recompute at the next formulate). Both then route through the single
`calc_condsat`/TVK recompute path. If B is added, document clearly that it does
**not** fire for pointer writes, so users are not misled into relying on it.

### 2.4 Generalization

The same pattern generalizes to other derived quantities (e.g. anything else
computed once from input at AR). The trigger flag and/or set-handler approach,
funneling into the package's existing recompute routine, is the template.

## Phasing

1. **Metadata core**: `api_class` on `MemoryType`, the enum, and
   `mem_set_api_class`; classify NPF + STO as the pilot.
2. **BMI exposure**: `get_var_api_class`; make input/output var-name lists
   meaningful.
3. **`K` → `CONDSAT` recompute**: Approach A (`IRECALCCONDSAT`), reusing the TVK
   recompute path; classify `K*` as `INERT_INPUT` with a note pointing at the
   trigger, and `CONDSAT` as `DERIVED`.
4. **Doc generation**: emit `doc/mf6api/` tables from the metadata.
5. **Roll out** classification to the remaining GWF packages, then GWT/GWE/PRT.
6. **(Optional) enforcement**: warn when an `OVERWRITTEN`/`STRUCTURAL` variable
   is set through the API — best-effort only, on the `set_value` path (see §1.6;
   pointer writes cannot be observed).

## Testing

- Unit: a variable's `api_class` is set as expected after `*_ar`.
- API integration (extend `autotest/test_gwf_libmf6_apimod.py`):
  - setting `K` then `IRECALCCONDSAT=1` changes heads the same as directly
    scaling `CONDSAT`; the flag is reset to `0` after the recompute;
  - setting `K` *without* the trigger still has no effect (documents the
    contract);
  - XT3D-active and TVK-active models behave correctly (no double-recompute, no
    conflict with TVK's own per-step recompute).
- Regression: confirm TVK results are unchanged (the trigger reuses the same
  path).

## Level of effort and tooling

Scale of the codebase (for calibration): ~1,693 `mem_allocate` calls (572 in
GWF alone), 151 files in `src/Model`, `MemoryManager.f90` ~2,800 lines, and one
existing memory unit test (`autotest/TestMemoryStore.f90`).

The work splits into three buckets with very different character:

| Bucket | Scope | Effort (experienced MF6 dev) | Nature |
|---|---|---|---|
| A. Core infrastructure | `api_class` field on `MemoryType` + enum; `mem_set_api_class`; BMI `get_var_api_class`; filter `get_input/output_var_names` | 2–4 days | Mechanical |
| B. Classifying variables | assign a class to each meaningful variable | pilot (NPF/STO/stress/HFB) ~1 week; full GWF +1–2 weeks; GWT/GWE/PRT/utilities/solution +2–4 weeks | Judgment-heavy |
| C. Thorough testing | mechanism unit tests + per-package API behavior tests + regression | mechanism 1–2 days; behavioral coverage scales with B, ~2–4 weeks for major packages | Mostly mechanical, high volume |
| D. Doc generation (§1.5) | generator + build wiring | 3–5 days | Mechanical, deferrable |

**Rough totals:** a credible, tested **pilot** (the packages already explored —
NPF, STO, WEL/RIV/GHB/DRN, HFB) in **~1–2 weeks**; a **thorough, repo-wide,
well-tested** implementation in **~2–3 months** of expert effort (not
necessarily full-time). These are ranges, not commitments; the classification
breadth is open-ended (an audit that is never strictly "done"), so a first
deliverable should be anchored on the pilot scope.

The estimate is **dominated by bucket B, which is not mechanical**:

- Each variable's class requires tracing where it is written/read across the
  lifecycle (AR/RP/AD/CF/CQ).
- **Option-dependence** is the subtle trap: a class changes with active features
  (`CONDSAT` is derived normally, overwritten under TVK/VSC, absent under XT3D).
  Whether to store a configured class or a base class plus documented
  conditionals is a per-variable decision.
- The stakes are high: a wrong "settable" label causes silent wrong results, so
  each classification should be verified — ideally **empirically** (set it
  through the API and confirm the behavior).

Infrastructure note: the `mem_set_api_class` tagging-helper route (§1.3b) keeps
bucket A small; adding an optional argument to `mem_allocate` instead is more
work because of that routine's many type/rank overloads. The memory manager is
also foundational (parallel/PETSc, NetCDF export, serialization), so the field
addition needs regression care on those paths.

### Suitability for AI coding assistants

- **A (infra), C (tests), D (doc-gen): excellent fit.** Well-scoped,
  pattern-following, low-judgment — agents scaffold the field/enum/BMI function,
  generate the test matrix, and write the empirical probes reliably.
- **B (classification): strong accelerator, human-in-the-loop required.** An
  agent can read each package, reason about the lifecycle, draft a
  classification, *and verify it empirically* by running the model through the
  API. It still errs on subtle option-dependence, and because a wrong label is
  silently incorrect, an MF6 expert must review each package's output. This is
  "agent drafts and verifies, expert signs off," not "agent autonomous."

**Net:** AI assistance likely compresses wall-clock by ~2–3×, mostly by crushing
A/C/D and by drafting and probing B quickly; bucket B's review burden stays with
a human and is the long pole. The practical approach is to do A first, then work
**package by package** through B with an empirical API test written for each
classification (so the test *is* the verification), with expert review of each
package — the same loop used to classify NPF, STO, the list stress packages, and
HFB for the `doc/mf6api/` guide.

## Open questions

- Should the recompute trigger live on NPF only, or be a general
  "model needs reformulation" signal?
- Class is option-dependent (TVK/VSC/XT3D). Do we store the *configured* class,
  or a base class plus documented conditionals?
- Interaction with HFB, which also rewrites `CONDSAT` at each Read and Prepare
  (a user-set or recomputed `CONDSAT` on a barrier connection is overwritten by
  HFB) — the metadata/notes should capture this.
- For Approach B, how to coalesce multiple `K`/`K22`/`K33` sets into one
  recompute (debounce via the flag).
