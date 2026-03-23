# Distributing MODFLOW 6

This document describes how to release MODFLOW 6. This folder contains scripts used to create distributions.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Modes](#modes)
  - [Develop](#develop)
  - [Release](#release)
- [Versions](#versions)
  - [Patch](#patch)
  - [Minor](#minor)
- [Steps](#steps)
  - [Review features](#review-features)
  - [Review deprecations](#review-deprecations)
  - [Review release notes](#review-release-notes)
  - [Release examples repo](#release-examples-repo)
  - [Create a release branch](#create-a-release-branch)
  - [Build assets/distributions](#build-assetsdistributions)
  - [Merge release branch to master](#merge-release-branch-to-master)
  - [Publish the release](#publish-the-release)
  - [Reset the develop branch](#reset-the-develop-branch)
    - [Update version strings](#update-version-strings)
    - [Update release notes](#update-release-notes)
  - [Release downstream repos](#release-downstream-repos)
- [Scripts](#scripts)
  - [Updating version numbers](#updating-version-numbers)
  - [Regenerating build files](#regenerating-build-files)
  - [Collecting deprecations](#collecting-deprecations)
  - [Benchmarking examples](#benchmarking-examples)
  - [Building PDF documents](#building-pdf-documents)
  - [Building distributions](#building-distributions)
  - [Checking distributions](#checking-distributions)
  - [Testing scripts](#testing-scripts)
- [Workflows](#workflows)
  - [Workflow triggers](#workflow-triggers)
  - [Testing workflows](#testing-workflows)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Modes

There are two kinds of distribution: provisional development builds and standard distributions approved for release.

The release-related tooling has two corresponding modes: develop mode and release mode. This document uses the same terminology.

### Develop

Develop-mode distributions contain only input/output documentation, release notes, `code.json` metadata, and binaries:

- `mf6`: MODFLOW 6 executable
- `zbud6`: Zonebudget executable
- `mf5to6`: MODFLOW 5 to 6 converter executable
- `libmf6`: MODFLOW 6 dynamic library

Binaries are built in develop mode, with feature flags disabled, making prerelease features accessible.

Running `mf6 -d` shows language describing the preliminary state of the software.

### Release

Release-mode distributions contain everything in development builds, plus:

- Fortran source code
- Meson files, makefiles, MSVS project files
- MF6 example models
- extensive docs:
  - MF6 input/output guide
  - example model docs
  - release notes
  - supplementary technical information
  - docs for `mf5to6` and `zbud6`

Binaries are built in release mode, with feature flags enabled, making prerelease features inaccessible. Prerelease features are filtered out of the generated documentation.

Running `mf6 -d` shows language affirming that the software has been reviewed and approved for distribution.

## Versions

MF6 loosely follows a modified form of [semantic versioning](https://www.semver.org).

### Patch

Patch releases contain bugfixes and behavior changes related to correctness.

Patch releases may be made in one of two ways.

- Branch from `master` and cherry-pick commits. "Hotfix" release.
- Branch from `develop` and use feature flags: "Flagged" release.

Hotfixes are quick and suitable for urgent, narrowly-scoped patches. Flagged releases need more careful preparation, but [trunk-based development](https://www.atlassian.com/continuous-delivery/continuous-integration/trunk-based-development) allows incremental development of features on the mainline.

### Minor

Minor releases may contain bugfixes, changes and new features. Minor releases branch from `develop`.

There are no major releases as the MF6 major version number is constant. Breaking changes are made in minor releases.

## Steps

If building distributions locally, a development environment [must be configured](../DEVELOPER.md) including a Fortran compiler, a Python environment, a LaTeX environment, and the [`MODFLOW-ORG/usgslatex`](https://github.com/MODFLOW-ORG/usgslatex) styles.

If using GitHub Actions, nothing is needed besides `git`.  This section assumes releases are made with GitHub Actions.

To make a release,

0. freeze development
1. review features
2. review deprecations
3. review release notes
4. release examples repo
5. create a release branch
6. build assets/distributions
7. merge release branch to master
8. publish the release
9. reset the develop branch
10. release downstream repos

Complete steps 1-3 in consultation with the development team, then steps 4 and 5 once the release is greenlit. Step 5 triggers automation for steps 6-8. Step 6 happens automatically; steps 7 and 8 require review and manual sign-off. Steps 9 and 10 are performed manually.

It is typical to undergo several iterations of 5-6 as candidate distributions are reviewed and issues are identified and resolved. The MF6IO guide in particular should be carefully inspected. Some things to look for:

- Have input/output samples been substituted into the appropriate sections? Is the content correct?
- Are all expected packages mentioned and described in the relevant tables? Some tables are auto-generated, some aren't. Make sure manually managed tables have been modified if necessary.
- Are all expected package parameters present? Parameters may have been guarded with `developmode true` in DFN files during development, make sure this attribute has been removed if parameters are to be released.

### Freeze development

Pause merging of new features to the mainline for 1-2 weeks prior to the anticipated release date. Feature development can continue in the meantime, but keeping the `develop` branch stable helps to minimize churn and avoid last-minute surprises.

### Review features

Determine whether any new developments need feature flags. To help keep the trunk prepared for a prompt release, it may be convenient to guard code and docs with feature flags by convention until the new features are ready for release.

Code can be guarded with the `developmode` routine in `FeatureFlagsModule`. DFN variables may be guarded with an attribute `developmode true`. LaTeX sections may be wrapped in `\ifdevelopmode ... \fi`

### Review deprecations

Before proceeding with a release, check for deprecated DFN variables due for removal. See the [developer docs](../DEVELOPER.md#deprecation-policy) for more info on deprecation policy and how to search DFNs for deprecated variables. Deprecated/removed variables are automatically detected and inserted into a table in the release notes.

### Review release notes

Double-check release notes in `doc/ReleaseNotes/develop.toml` with the authors of any changes to be included in the release.

For hotfix releases, `develop.toml` must be trimmed manually on the release branch. For patch releases made from `develop`, release notes are automatically filtered to include only fixes.

**Note**: For all releases, add a line to the Release History section of `ReleaseNotes.tex` providing the version number, date and DOI of the release, e.g. `6.4.4 & February 13, 2024 & \url{https://doi.org/10.5066/P9FL1JCC}`. DOIs are updated with minor releases and remain the same for patch releases.

### Release examples repo

MODFLOW 6 [example models](https://github.com/MODFLOW-ORG/modflow6-examples) are bundled with official releases. The `release.yml` workflow attempts to download the latest release from the examples repository. An examples release is typically made before proceeding with an MF6 release (this may not be necessary if example models have not changed since the last release). See the [examples release instructions](https://github.com/MODFLOW-ORG/modflow6-examples/blob/develop/DEVELOPER.md#releasing-the-examples) for more information.

### Create a release branch

Create a release candidate branch from `develop` or `master`. The branch's name must begin with `v` followed by the version number. For a dry run in which a candidate distribution is built but the release does not proceed, append `rc` to the version string, e.g. `v6.4.0rc`. For an approved release, include *only* the version number.

```shell
git checkout develop
git switch -c v6.4.0
```

Push the branch to the repository. This triggers the release workflow. 

### Build assets/distributions

GitHub actions runs the release workflow and builds binaries, documentation and distribution archives for supported platforms.

If this is a dry run (branch name ends with `rc`) binaries are built with `IDEVELOPMODE` set to 1, and the workflow ends after uploading artifacts for review. If this is not a dry run, the workflow will continue, drafting a pull request against the `master` branch after the build completes.

### Merge release branch to master

Review the distributions. If they pass inspection, merge (don't squash) the pull request from the release branch into `master`. (Squashing will cause `master` and `develop` to diverge.) This triggers another job to tag the new tip of `master` with the release number, draft a release, and upload binaries and documentation as release assets.

### Publish the release

Inspect the draft release. The generated description should be in the form:

```
This is the approved USGS MODFLOW <semver> release.

<authors>, <release year>, MODFLOW 6 Modular Hydrologic Model version <semver>: U.S. Geological Survey Software Release, <release date>, <doi link>

Visit the USGS "MODFLOW and Related Programs" site for information on MODFLOW 6 and related software: https://doi.org/10.5066/F76Q1VQV
```

Update the DOI link in the citation if necessary. The DOI on the last line is the original and stays the same.

Publish the release.

### Reset the develop branch

Make a new branch from `master`:

```shell
git checkout master
git switch -c post-6.x.y-release-reset
```

#### Update version strings

Update the version number for the next development cycle:

```shell
pixi run update-version -v 6.x.y.dev0
```

This will substitute the new version number into the necessary files and set `IDEVELOPMODE` back to 1.

#### Update release notes

Generate a `develop.tex` file from `develop.toml`:

```shell
pixi run make-release-notes
```

Move/rename it to `doc/ReleaseNotes/previous/vx.y.z.tex` (where `x.y.z` is the version just released), then insert a new line `\input{./previous/vx.y.z.tex}` at the top of `doc/ReleaseNotes/appendixA.tex`.

If this was not a hotfix, trim `doc/ReleaseNotes/develop.toml` as necessary to remove items just released.

Create and merge (don't squash) a pull request from this branch into `develop`.

### Release downstream repos

MODFLOW 6 releases are typically followed by releases of [flopy](https://github.com/modflowpy/flopy), [pymake](https://github.com/modflowpy/pymake) and the combined [executables](https://github.com/MODFLOW-ORG/executables) distribution.

To release flopy, follow the steps in that repository's [developer documentation](https://github.com/modflowpy/flopy/blob/develop/docs/make_release.md).

To release pymake:

1. Update `usgsprograms.txt` with the path to the new MODFLOW 6 release.
2. Update other targets in `usgsprograms.txt` with the path to new releases.
3. Release a new version.

To trigger an executables release, follow the steps that repository's [developer documentation](https://github.com/MODFLOW-ORG/executables/blob/master/DEVELOPER.md#triggering-a-release).

## Scripts

This directory contains scripts for

- updating version numbers
- regenerating build files
- collecting deprecations
- benchmarking examples
- building PDF documents
- building distributions
- verifying distributions

### Updating version numbers

MODFLOW 6 version numbers follow the [semantic versioning](https://semver.org/) convention `major.minor.patch`. Release tags do not include an initial `v`.

The version string is stored in `version.txt` in the project root. The version string appears in several other files in the repository, as well as date and timestamp information.

The `update_version.py` script synchronizes updates to `version.txt` and other files containing version information.

```shell
pixi run update-version -v 6.4.1
python update_version.py -v 6.4.1 # or from the distribution/ folder
```

If a `--version` value is not provided, the version string will not be changed, just dates and timestamps. The `--version` value may contain trailing letters, e.g.

```shell
python update_version.py -v 6.4.2rc
```

These must begin immediately after the patch version number, and may contain numerics after an initial alphabetic character.

### Regenerating build files

Up-to-date makefiles must be generated for inclusion in a distribution. The `build_makefiles.py` script rewrites makefiles after Fortran source files have been added, removed, or renamed. 

```shell
pixi run build-makefiles
python build_makefiles.py # or from the distribution/ folder
```

MSVS project files are also included in the distribution. These currently cannot be auto-generated and must be updated/verified by hand.

### Collecting deprecations

Deprecated/removed variables are scraped from DFNs to create a table for insertion into the release notes:

```shell
pixi run collect-deprecations
python deprecations.py # or from the doc/mf6io/mf6ivar/ folder
```

### Benchmarking examples

The `benchmark.py` script benchmarks the current development version of MODFLOW 6 against the latest release rebuilt in development mode, using the models from the `MODFLOW-ORG/modflow6-examples` repository.

Paths to pre-built binaries for both versions can be provided via the `--current-bin-path` (short `-c`) and `--previous-bin-path` (short `-p`) command line options. If bin paths are not provided, executables are rebuilt in the default locations:

- `<project root>/bin`: current development version
- `<project root>/bin/rebuilt`: previous version

The examples repository must first be installed and prepared as described above. Its path may be explicitly provided with the `--examples-repo-path` (short `-e`) option. If no path is provided, the repository is assumed to be named `modflow6-examples` and live side-by-side with the `modflow6` repository on the filesystem.

```shell
pixi run benchmark
python benchmark.py # or from the distribution/ folder
```

### Building PDF documents

The `build_docs.py` script constructs documentation.

- regenerates LaTeX files from DFN files (via `mf6ivar.py`)
- downloads or reruns benchmarks (via `benchmark.py`)
- downloads publications hosted on the USGS website
- builds the MF6IO PDF document

```shell
pixi run build-docs
python build_docs.py # or from the distribution/ folder
```

The script is lazy &mdash; files are regenerated only if they do not already exist or the `--force` flag is provided. Likewise for installing example models and running benchmarks.

### Building distributions

The `build_dist.py` script constructs a complete distribution.

- builds binaries
- builds examples
- builds documentation
- collects release assets
- creates a distribution archive

```shell
pixi run build-dist -o $DISTDIR
python build_dist.py -o $DISTDIR # or from the distribution/ folder
```

The script is lazy &mdash; components of the distribution are recreated only if they do not exist or the `--force` flag is provided.

### Checking distributions

The `check_dist.py` script runs some checks on the distribution, e.g.:

- binaries are present
- source code is present
- PDF documents are present
- binaries successfully run example models 
- binaries emit expected version strings and other output
- build files are present and can successfully build the source code

```shell
pixi run check-dist --path $DISTDIR
pytest -v -s check_dist.py --path $DISTDIR # or from the distributions/ folder
```

### Testing scripts

Each script in `distribution/` contains its own tests. To run them, run `pytest` from the `distribution/` folder. (This happens in standard CI.) The tests will not be discovered if `pytest` is run from a different location, as the scripts in this folder are not named `test_*.py` and are only discoverable by virtue of the patterns provided in `distribution/pytest.ini`.  The tests use temporary directories where possible and revert modifications to tracked files on teardown.

There is a small additional suite of tests that can be used to validate a release distribution folder after it is built: `check_dist.py`. These tests are run as part of the release workflow &mdash; see below for more detail.

**Note:** the tests clean up after themselves by reverting changes to files in the following locations:

- `doc/`
- `make`
- `utils/**/make/`

Make sure you don't have any uncommitted changes in these locations before running the tests in your environment.

**Note:** to avoid contested file access, the tests will refuse to run in parallel with `pytest-xdist`.

## Workflows

The scripts in this directory are used by two GitHub Actions workflows:

- `.github/workflows/release.yml`: "callable" workflow for development or standard releases
- `.github/workflows/release_dispatch.yml`: dispatches a release on various triggers

### Workflow triggers

The `release.yml` workflow has no triggers of its own, and must be dispatched by `.github/workflows/release_dispatch.yml`, in one of two ways:

- Pushing a branch with a suitable name (e.g. `vx.y.z[rc]`) to the `MODFLOW-ORG/modflow6` repository. This is how releases are typically triggered.
- Triggering the workflow via GitHub CLI or web UI. Useful for testing release candidates or verifying the release automation before a final release is made.

The `release.yml` workflow is a callable function for producing distributions. It uses the scripts in this directory to build a distribution for each supported platform. Custom actions in `.github/actions/` are also used for extended builds.

This workflow is also used by the [nightly build repository](https://github.com/MODFLOW-ORG/modflow6-nightly-build).

### Testing workflows

The `workflow_dispatch` event is GitHub's mechanism for manually triggering workflows. This can be accomplished from the Actions tab in the GitHub UI. This is a convenient way to test the release procedure and evaluate release candidate distributions.

To dispatch the release workflow, navigate to the Actions tab of this repository. Select the release dispatch workflow. A `Run workflow` button should be visible in an alert at the top of the list of workflow runs. Click the `Run workflow` button, providing suitable values for the inputs.