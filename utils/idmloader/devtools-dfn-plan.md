# Plan: update modflow-devtools `dfn` branch for IDM code generation

When devtools v2 ships, `dfn2f90.py` should be updated to use the new `dfn`
package instead of `_load_v1_flat`. This document tracks the two small changes
needed in devtools to make that possible, plus the corresponding update to
`dfn2f90.py`.

---

## Context

The `dfn` branch restructures `modflow_devtools.dfn` into a package with a
proper `FieldV1` dataclass that already preserves every IDM-critical attribute
(`in_record`, `preserve_case`, `mf6internal`, `layered`, `block_variable`,
`developmode`). The `MapV1To2` stripping of those attributes is an explicit
opt-in for the flopy path — it does not affect a caller that uses schema v1.

Two small gaps prevent `dfn2f90.py` from using `FieldV1` directly today.

---

## Changes needed in devtools (`dfn` branch)

### 1. Add `time_series` to `FieldV1`

**File:** `modflow_devtools/dfn/schema/v1.py`

`time_series` is a DFN attribute used to mark IDM timeseries parameters. It is
currently silently dropped by `FieldV1.from_dict(strict=False)`.

```python
# before
@dataclass(kw_only=True)
class FieldV1(Field):
    valid: tuple[str, ...] | None = None
    reader: Reader = "urword"
    tagged: bool = False
    in_record: bool = False
    layered: bool | None = None
    preserve_case: bool = False
    numeric_index: bool = False
    deprecated: bool = False
    removed: bool = False
    mf6internal: str | None = None
    block_variable: bool = False
    just_data: bool = False

# after — add one line
    time_series: bool = False
```

### 2. Capture `# mf6 subpackage` comments in `parse_dfn`

**File:** `modflow_devtools/dfn/parse.py`

Six DFN files (`gwf-dis`, `gwf-disv`, `gwt-dis`, `gwt-disv`, `gwe-dis`,
`gwe-disv`) declare their subpackage using `# mf6 subpackage <abbr>` rather
than the flopy-style `# flopy subpackage <key> <abbr> <param> <val>`. The
current parser only captures `# flopy ...` comments so these are lost.

```python
# in parse_dfn(), in the comment-line handling block:

# before
if line.startswith("#"):
    _, sep, tail = line.partition("flopy")
    if sep == "flopy":
        if (
            "multi-package" in tail
            or "solution_package" in tail
            or "subpackage" in tail
            or "parent" in tail
        ):
            metadata.append(tail.strip())
    _, sep, tail = line.partition("package-type")
    if sep == "package-type":
        metadata.append(f"package-type {tail.strip()}")
    continue

# after — add mf6 subpackage capture
if line.startswith("#"):
    _, sep, tail = line.partition("flopy")
    if sep == "flopy":
        if (
            "multi-package" in tail
            or "solution_package" in tail
            or "subpackage" in tail
            or "parent" in tail
        ):
            metadata.append(tail.strip())
    _, sep, tail = line.partition("package-type")
    if sep == "package-type":
        metadata.append(f"package-type {tail.strip()}")
    _, sep, tail = line.partition("mf6 subpackage")
    if sep == "mf6 subpackage":
        metadata.append(f"mf6-subpackage {tail.strip()}")
    continue
```

Then add a helper alongside the existing ones in `parse.py`:

```python
def parse_mf6_subpackages(meta: list[str]) -> list[str]:
    """
    Return MF6 subpackage abbreviations declared via '# mf6 subpackage <abbr>'.
    These are distinct from flopy subpackages ('# flopy subpackage ...').
    """
    result = []
    for m in meta:
        if m.startswith("mf6-subpackage "):
            abbr = m.removeprefix("mf6-subpackage ").strip().upper()
            result.append(abbr)
    return result
```

---

## Corresponding update to `dfn2f90.py`

Once the devtools changes above are released, `parse_dfn` in `dfn2f90.py` can
be replaced with a thin wrapper over `modflow_devtools.dfn.load()`.

The IDM-specific logic that stays in `dfn2f90.py`:
- `_normalize_type()` / `_BASE_TYPE_MAP` — IDM Fortran type strings are an
  MF6 code-generation concern, not a devtools concern
- Shape processing (`NAUX NCPL`, `EXG CELLIDM` overrides, etc.)
- Longname wrapping for Fortran line length
- `Param` / `Block` / `DfnFile` dataclasses (or inline template access)
- Block-required inference logic

Sketch of the new `parse_dfn`:

```python
from modflow_devtools.dfn import load as load_dfn
from modflow_devtools.dfn.parse import parse_mf6_subpackages
from modflow_devtools.dfn.schema.v1 import FieldV1

def parse_dfn(dfnfspec: Path) -> DfnFile:
    component, subcomponent = dfnfspec.stem.upper().split("-")

    with dfnfspec.open(encoding="utf-8") as f:
        dfn = load_dfn(f, name=dfnfspec.stem, format="dfn")

    multi_package = dfn.multi
    mf6_subpkgs = parse_mf6_subpackages(...)  # needs meta exposed — see note below
    subpackages = [s.ljust(16) for s in mf6_subpkgs] or [" " * 16]

    block_names_ordered = []
    block_data = {}
    params = []

    for field in dfn.fields.values(multi=True):  # FieldV1 objects
        blockname_upper = (field.block or "").upper()
        if not blockname_upper:
            continue

        # ... same block tracking logic as today ...

        if field.block_variable:
            block_data[blockname_upper]["has_block_var"] = True
            continue

        vn = field.name.upper()
        mf6vn = (field.mf6internal or field.name).upper()
        shape = _process_shape(field.shape or "", component, vn, mf6vn)
        shapelist = shape.strip().split() if shape.strip() else []
        t = _normalize_type(field.type or "", " ".join(shapelist), len(shapelist),
                            (field.type or "").startswith("recarray"))

        # longname, required, etc. same as today but reading from FieldV1 attributes
        required = not field.optional
        in_record = field.in_record
        preserve_case = field.preserve_case
        layered = field.layered or False
        timeseries = field.time_series  # available after devtools change #1
        developmode = field.developmode

        params.append(Param(...))

    # ... rest of block building unchanged ...
```

**Note:** `load()` in the `dfn` branch does not currently expose the raw
metadata list after building the `Dfn`. To access `mf6-subpackage` metadata
(change #2 above), either:
- Call `parse_dfn(f)` directly for metadata, then `load()` for the `Dfn`, or
- Add a `metadata: list[str]` field to `Dfn` populated by `load()` —
  the cleaner option, worth suggesting in the devtools PR.

---

## Summary of devtools changes

| File | Change |
|------|--------|
| `modflow_devtools/dfn/schema/v1.py` | Add `time_series: bool = False` to `FieldV1` |
| `modflow_devtools/dfn/parse.py` | Capture `# mf6 subpackage` in metadata; add `parse_mf6_subpackages()` helper |
| `modflow_devtools/dfn/__init__.py` | (optional) Add `metadata: list[str]` to `Dfn` dataclass, populated by `load()` |
