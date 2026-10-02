# SPF Seepage Face (UZR) — Action Items

Status: Option A implemented, builds clean, all 5 cases in
`autotest/test_gwf_uzr_spf.py` pass, Fortran format check passes. Nothing committed.

## What was done
- `doc/mf6io/mf6ivar/dfn/gwf-spf.dfn`: period input is now exchange-style
  connection data — `cellid ihc cl1 hwva angldegx` (replaced `dist`/`area`).
- `src/Idm/gwf-spfidm.f90`: regenerated from the DFN.
- `src/Model/GroundWaterFlow/gwf-spf.f90`: removed the `condsat(node)` hack;
  `effective_k()` builds the face normal from `ihc` + `angldegx` and calls
  `hyeff(...)` for a directional, anisotropy-aware conductance (single cell, no
  two-cell averaging). Conductance = `effective_k * hwva / cl1`; seepage switch
  uses `z = (bot + top) / 2` as the face datum.
- `autotest/test_gwf_uzr_spf.py`: `spf_data` rows updated to the new columns.
- Regenerated flopy (`mfgwfspf.py`) and mf6ivar docs.

## Must-not-forget action items

### 1. Understand why the `bnd_ar` override was never dispatched (HIGH)
- The first implementation cached NPF pointers in a `bnd_ar => spf_ar` override.
  It was never called by the model's AR path, so `k11` stayed null and the model
  segfaulted the instant a seepage face first activated.
- Workaround in place: lazy `set_npf_pointers()` guarded by
  `associated(this%k11)`, called at the top of `spf_cf` (matches the file's
  original `mem_setptr`-in-`cf` idiom).
- Action: confirm *why* `bnd_ar` isn't invoked for this package before relying on
  an AR override anywhere else. May indicate a BndExt vs. Bnd AR-path nuance.
- Related: the centralization proposal in
  [SPF_EFFECTIVE_K_PROPOSAL.md](SPF_EFFECTIVE_K_PROPOSAL.md) would replace the
  cached NPF pointers with a single injected NPF handle, removing this concern.
- In progress: the abstract provider is being built on branch
  `conductance-provider` (worktree `../modflow6-conductance-provider`, see its
  `CONDUCTANCE_PROVIDER.md`), to merge first; SPF then consumes
  `this%cond_provider%eff_hy(...)` and drops `effective_k`/`set_npf_pointers`.

### 2. Face elevation datum `z` for `ihc = 0` (MEDIUM)
- `z = (bot + top) / 2` is correct for a full-height vertical face (`ihc = 1`),
  but only approximate for a true horizontal (top/bottom) face.
- Action: decide whether `ihc = 0` needs `z = top` (or `bot`), or an explicit
  datum, if horizontal seepage faces become a real use case.

### 3. Mover support (MEDIUM)
- `spf_fc` still hard-stops with an error when `imover == 1` (`TODO_UZR`).
- Action: implement seepage-to-mover, or document that it is unsupported.

### 4. Placeholder option `some_option` / `DEV_VS2D_BND` (LOW)
- Still a stub carried over from the original scaffold.
- Action: remove it or give it real meaning.

### 5. `angldegx` cannot express a sloped (non-axis-aligned vertical) normal (LOW)
- Known Option A limitation accepted during design: `angldegx` only gives the
  horizontal azimuth; a face normal with a vertical component is not expressible.
- Action: revisit only if sloped seepage faces are needed (would require a full
  normal vector — the Option B route).

### 6. Pre-commit housekeeping (before opening a PR)
- `codespell` (spelling), `ruff` format + lint on the Python test.
- Confirm `src/Idm/gwf-spfidm.f90`, flopy `mfgwfspf.py`, and mf6ivar docs are all
  regenerated and consistent with the DFN.
- Run the wider UZR test set, not just `test_gwf_uzr_spf.py`.
- Review the `.dfn` descriptions/longnames for `ihc`, `cl1`, `hwva`, `angldegx`.

## Key files
- `src/Model/GroundWaterFlow/gwf-spf.f90`
- `doc/mf6io/mf6ivar/dfn/gwf-spf.dfn`
- `src/Idm/gwf-spfidm.f90`
- `autotest/test_gwf_uzr_spf.py`
- `src/Utilities/HGeoUtil.f90` (`hyeff`) — the shared anisotropic-K primitive
