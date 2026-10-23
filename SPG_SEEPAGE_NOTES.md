# SPG (Seepage) Package — Status Handoff

_Last updated: 2026-10-20_

## Where we are
The **SPG seepage boundary package** has been reworked from the old prototype
into a Richards-only seepage face that approximates conductance from NPF
conductivity × relative permeability. It **builds cleanly** and its **autotest
passes** (2 scenarios).

## Current SPG behavior
- **Applies only to Richards (UZR) cells.** Requires the UZR package active;
  validated per-cell via UZR's `IUNSAT` array.
- **Input** (`stress_period_data`): `cellid  ihc  cl1  hwva  angldegx`
  - `ihc` — 0 = vertical face, 1 = horizontal face
  - `cl1` — distance from cell center to the face
  - `hwva` — seepage face area
  - `angldegx` — **reserved/unused** for now (kept in the definition)
- **Conductance:** `cond = krel(node) · K · hwva / cl1`, LAK-style:
  - vertical face (`ihc=0`) → `K33` (falls back to `K11` if `IK33=0`)
  - horizontal face (`ihc=1`) → `K11`
  - `krel` = NPF relative permeability (`KREL`), hardwired via `mem_setptr`
- **Seepage switch:** discharges (Robin term, `hcof=-cond`, `rhs=-cond·z`) while
  `head > z` (cell-center elevation); no-flow otherwise.
- The old penalty / ibound-toggle machinery and the ConductanceProvider
  abstraction were **removed**.

## Key files
- `src/Model/GroundWaterFlow/gwf-spg.f90` — package: `effective_k`, `spg_cf`,
  `spg_ck` (Richards-cell check), pointer setup/teardown.
- `doc/mf6io/mf6ivar/dfn/gwf-spg.dfn` — input definition.
- `src/Idm/gwf-spgidm.f90` — generated IDM (via
  `utils/idmloader/scripts/dfn2f90.py`).
- `src/Model/GroundWaterFlow/gwf.f90` — model wiring.
- `autotest/test_gwf_spg.py` — representative UZR+SPG test (active / no-flow),
  asserts an analytic 2-node Robin solution.

## Also removed this session (by request)
- **SGD package** deleted entirely (source, IDM, dfn, test, build/registration).
- **ConductanceProvider** + **NpfConductanceProvider** modules deleted; NPF's
  now-unused `calc_eff_hy` removed too.

## How to build & test
```bash
# build
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
ninja -C builddir

# regenerate IDM after any dfn change
cd utils/idmloader/scripts && python dfn2f90.py --outdir ../../../src/Idm

# refresh flopy classes after any dfn change (for the test)
cd autotest && python update_flopy.py

# run the test (binary staged at bin/mf6)
cd autotest && python -m pytest test_gwf_spg.py -q
```
Note: `bin/mf6` is a gitignored staged copy of `builddir/src/mf6`.

## Open items / next steps
1. **`DKRDH` Jacobian term** — `krel` is currently lagged; adding `dkr/dh`
   (NPF stores `DKRDH`) would make the Newton term consistent. _Pending._
2. **`angldegx` decision** — still reserved; decide whether to use it
   (horizontal anisotropy / full tensor) or drop it.
3. **Convergence near the switch** — the `head > z` on/off combined with
   head-dependent `krel` can chatter; watch in richer models.
4. **Partial-domain Richards** — UZR today requires all cells `iunsat=1`, so the
   per-cell SPG check mainly guards "UZR absent" until partial-domain support
   exists.

## Known caveats (physics/numerics)
- Grid-aligned `K11`/`K33` only — ignores rotated anisotropy (`ANGLE1/2/3`) and
  XT3D.
- Fixed user `hwva`/`cl1` — wetted area doesn't change as heads move.
- Uses the seepage cell's own `krel` (no upstream weighting).
