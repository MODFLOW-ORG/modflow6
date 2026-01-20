# Contributing to MODFLOW 6

Contributions to MODFLOW 6 are welcome. We ask that contributors follow some guidelines.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Conduct](#conduct)
- [Submissions](#submissions)
  - [Bugs](#bugs)
  - [Questions](#questions)
  - [Requests](#requests)
- [Source code](#source-code)
  - [Format](#format)
  - [Style guide](#style-guide)
- [Commit messages](#commit-messages)
  - [Template](#template)
  - [Components](#components)
    - [Type](#type)
    - [Scope](#scope)
    - [Subject](#subject)
    - [Body](#body)
    - [Footer](#footer)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Conduct

Help us keep MODFLOW 6 open and inclusive. Please read and follow our [Code of Conduct](./CODE_OF_CONDUCT.md).

## Submissions

### Bugs

If you find a bug, you can help us by submitting an issue. Even better, you can submit a pull request with a fix.

Before submitting an issue, please search the issue tracker. Your issue may already have been reported, and the discussion might inform you of workarounds readily available.

We want to fix all issues as soon as possible, but we first need to reproduce and confirm them. We ask that all bug reports provide a minimal, complete, and verifiable example of the problem, including the following information:

- Operating system
- Version of MODFLOW 6
- A description of the problem contrasting actual with expected behavior
- A minimal recipe to reproduce the problem &mdash; e.g. Python/FloPy script or MF6 input files.

Please provide such an example. This saves maintainers time and ultimately allows us to fix more bugs. While we understand it might be difficult to isolate the essence of a problem, bugs can only be investigated or fixed if problems are accurately identified. Issues lacking sufficient information to reproduce may be closed.

File an issue by filling out our [issue form](https://github.com/MODFLOW-ORG/modflow6/issues/new).

### Questions

Please do not open issues for general support questions as we want to keep GitHub issues for bug reports and feature requests. Questions can be asked on the [discussions page](https://github.com/MODFLOW-ORG/modflow6/discussions) under the Q&A category or on [Stack Overflow](https://stackoverflow.com/questions/tagged/modflow6) with tag `modflow6`.

We will systematically close all issues that are requests for general support.

### Requests

You can request a new feature by submitting an issue to our GitHub Repository.

If you would like to implement a new feature:

- **Major** features should be discussed first. Please open an issue and outline your proposal. This will also allow us to coordinate our efforts, prevent duplication of work, and help you craft the change so that it can be accepted into the project.
- **Small** features can be submitted directly as a pull request.

To submit a pull request (PR):

1. To avoid duplicating effort, [search](https://github.com/MODFLOW-ORG/modflow6/pulls) for an open or closed PR that relates to your submission.
2. Fork the MODFLOW-ORG/modflow6 repo and make your changes in a new branch, following our style and commit message guidelines and [including appropriate test cases](./DEVELOPER.md#writing-tests).
3. [Check the spelling and formatting](./DEVELOPER.md#formatting) of any modified or new Fortran source files, python files definition files, markdown, and LaTeX files.
4. [Rebuild makefiles](./DEVELOPER.md#generating-makefiles) and update MSVS project files if you added, removed, or renamed any source files.
5. [Run the full test suite](./DEVELOPER.md#running-tests) and make sure all tests pass.
6. Push your branch to GitHub and create a pull request to the `develop` branch.
7. If we suggest changes:
  a. make the required updates
  b. make sure tests still pass
  c. rebase your branch if needed
  d. (force) push to update your PR

That's it! Thank you for your contribution!

If you have installed the pixi environment you can complete steps 3 and 4 using:

```shell
pixi run prepare-pull-request
```

## Source code

### Format

MODFLOW 6 source code is formatted automatically with the [fprettify formatter](https://github.com/pseewald/fprettify). Format rules are specified in the [fprettify configuration file](.fprettify.yaml). See the [developer docs](./DEVELOPER.md#formatting) for instructions on using `fprettify`.

### Style guide

Automated tools like `fprettify` cannot directly address all style concerns. Please review our [style guide](https://modflow6.readthedocs.io/en/latest/_dev/styleguide.html) before submitting your code.

## Commit messages

Adherence to commit message formatting rules results in consistent messages and a more readable project history.

Each commit message consists of three components: **header** (required), **body** (optional) and **footer** (optional). The header has three subcomponents: **type** (required), **scope** (optional) and **subject** (required):

See [the develop branch history](https://github.com/MODFLOW-ORG/modflow6/commits/develop) to get a sense for the commit style.

### Template

The general structure of a commit message is:

```text
<type>(<scope>): <subject>
<BLANK LINE>
<body>
<BLANK LINE>
<footer>
```

### Components

#### Type

Must be one of the following:

- **ci**: CI configuration files or scripts
- **docs**: Online or PDF documentation
- **feat**: New features
- **fix**: Bug fixes
- **perf**: Performance improvements
- **refactor**: Neither fixes a bug nor adds a feature, but may e.g. improve maintainability
- **style**: Non-functional format/style changes only
- **test**: New tests, corrections to existing tests, test refactors
- **revert**: Revert a previous commit &mdash; the commit message should identify the hash of the commit being reverted

#### Scope

The scope should identify the main target of the change. This might be, for instance:

- a filename
- an MF6 package name
- a Fortran module name
- a section of the documentation, e.g. "releasenotes" or "readme"

#### Subject

The subject contains a succinct description of the change.

Use the imperative, present tense: "change" not "changed" nor "changes". Avoid initial capitalization and do not end the line with a full stop.

#### Body

The body should include the motivation for the change and contrast previous and modified behavior.

Use standard capitalization and punctuation and break long explanations into suitable paragraphs. Single-backtick `code snippets` may be used.

Just as in the **subject**, use the imperative present tense: "change" not "changed" nor "changes".

#### Footer

If the commit closes an issue, the footer should contain a [closing reference](https://help.github.com/articles/closing-issues-via-commit-messages/).

If the commit includes breaking changes, the footer should note this prominently.
