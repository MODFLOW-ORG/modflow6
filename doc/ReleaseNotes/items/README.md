# Release notes

Add files to this directory to describe user-facing changes (bug fixes, new features, new examples, etc) to be included in the next release of MODFLOW 6.

To add a release note, create a `.toml` file with a short descriptive name, e.g. `wel-auto-flow-reduce-auxname.toml`:

```toml
section = "fixes"
subsection = "model"
description = "Plain-language description of the change..."
```

Valid values for `section` and `subsection` are defined in `doc/ReleaseNotes/schema.toml`.

Sections and subsections are rendered in the order listed in `schema.toml`. Within a subsection, items are ordered by file name.

Files are merged by the `mk_releasenotes.py` script underlying the `pixi run make-release-notes` task.
