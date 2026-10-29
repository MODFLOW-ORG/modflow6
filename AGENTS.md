# AGENTS.md

Guidance for AI agents working in the MODFLOW 6 repository. This file is intentionally
short: it points to the authoritative docs rather than duplicating them, so routine,
well-defined tasks (building, testing, formatting, _etc._) are not re-derived. When a
procedure lives in a linked doc, follow that doc instead of guessing.

## Primary references

- **[DEVELOPER.md](./DEVELOPER.md)** — the canonical developer guide. Start here for
  environment setup, building, formatting, linting, testing, branching, and deprecation.
- **[CONTRIBUTING.md](./CONTRIBUTING.md)** — contribution and pull request process.
- **[EXTENDED.md](./EXTENDED.md)** — building and testing the extended version.
- **[IDM.md](./IDM.md)** — input definitions / the Input Data Model and generated files.

## Common tasks (follow the linked section; do not re-derive commands)

- **Environment** — Pixi is the canonical workflow and matches CI. Use
  `pixi run install` to set up the dev environment.
  See [DEVELOPER.md › Install the python environment](./DEVELOPER.md#install-the-python-environment).
  (Conda via `environment.yml` is the alternative.)
- **Building** → [DEVELOPER.md › Building](./DEVELOPER.md#building) (Pixi/Meson).
- **Testing** (unit + integration) → [DEVELOPER.md › Testing](./DEVELOPER.md#testing).
- **Formatting & linting** (fprettify, ruff, codespell) →
  [DEVELOPER.md › Formatting](./DEVELOPER.md#formatting).
- **Branching model** → [DEVELOPER.md › Branching model](./DEVELOPER.md#branching-model).
- **Deprecation policy** → [DEVELOPER.md › Deprecation policy](./DEVELOPER.md#deprecation-policy).

## Quick commands

Canonical Pixi invocations for the most common actions (run from the project root
unless noted). These mirror the linked DEVELOPER.md sections; see there for options,
the raw Meson/Pytest equivalents, and the Conda workflow.

- **Set up the dev environment:** `pixi run install`
- **Configure the build directory:** `pixi run setup builddir`
  (add `-Ddebug=true` for a debug build)
- **Build MODFLOW 6 + ZONEBUDGET** (installs to `bin/`): `pixi run build builddir`
- **Run all unit tests:** `pixi run test builddir`
- **Run a single unit-test module** (name from `autotest/tester.f90`):
  `pixi run test builddir --verbose ArrayHandlers`
- **Run the integration tests:** `pixi run autotest`
- **Format a single Fortran file:** `pixi run fix-style path/to/file.f90`
- **Format / lint Python:** `pixi run fix-python-format` and `pixi run check-python-lint`
- **Spell check:** `pixi run check-spelling`
- **Regenerate after a `.dfn` change** (IDM Fortran + mf6io input docs):
  `pixi run update-fortran-definitions` and `pixi run run-mf6ivar` (see [IDM.md](./IDM.md)).

## Repository orientation

- `src/` — Fortran source for MODFLOW 6 (models, exchanges, solution, IDM, utilities).
- `srcbmi/` — BMI/XMI interface sources.
- `autotest/` — Python test suite (unit and integration tests).
- `utils/` — utility programs (_e.g._, ZONEBUDGET, the mf5to6 converter).
- `doc/` — documentation sources, release notes, and the LaTeX manuals.
- `doc/mf6io/mf6ivar/dfn/` — input definition (`.dfn`) files; code and docs are generated
  from these.
- `distribution/`, `.github/` — release tooling and CI configuration.
- `builddir/` (and other `_builddir*` folders) — disposable Meson build output;
  never edit by hand, just rebuild.

## Conventions and guardrails

- **Generated files:** do not hand-edit generated artifacts. Edit the source of truth and
  regenerate. This includes IDM-generated Fortran and the doctoc-generated table of
  contents in `DEVELOPER.md`.
- **`.dfn` descriptions:** keep `description` entries concise. The detailed explanation
  belongs in the manual and the technical reference, not in the definition file.
- **Fortran:** run `fprettify` on changed Fortran before committing
  (see [DEVELOPER.md › Fortran formatting](./DEVELOPER.md#fortran-formatting)).
- **Python:** run `ruff` format and check on changed Python before committing
  (see [DEVELOPER.md › Python formatting](./DEVELOPER.md#python-formatting)).
- **Branching/PRs:** follow the repository
  [branching model](./DEVELOPER.md#branching-model) when creating branches and pull
  requests.
- **Testing scope:** scope tests to what changed, and rebuild first (`pixi run build
  builddir`) — integration tests run the installed binary. Run the affected Fortran
  unit-test module(s) by name, or the specific integration file(s). Anything broader —
  a smoke run (`pixi run autotest -S`) or the full suite — is slow, so ask the developer
  to confirm (or run it themselves) before starting it; escalate only when the change is
  cross-cutting (IDM, solver, I/O, shared utilities) or a targeted run surfaces failures
  you can't bound. See [DEVELOPER.md › Testing](./DEVELOPER.md#testing).
- **Commit messages:** keep them concise. See [the commit section in CONTRIBUTING.md](./CONTRIBUTING.md#commit-messages).
- **Tone:** reserve compliments for genuinely exceptional finds or ideas.
