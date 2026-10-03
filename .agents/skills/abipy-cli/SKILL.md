---
name: abipy-cli
description: Add or modify AbiPy command-line tools and subcommands, preserving parser, dispatch, output, exit-code, optional-dependency, safety, packaging, and subprocess-test conventions. Use for scripts under abipy/scripts and their shared CLI helpers; not for domain logic that belongs in library modules.
---

# Develop AbiPy command-line tools

Keep the CLI thin: parse user intent, validate boundaries, call reusable AbiPy APIs, format results, and return a
meaningful exit code. Put scientific algorithms, readers, input construction, and workflow mechanics in their domain
modules so they are usable without a subprocess.

Use `abipy-navigate` to locate the underlying API, `python-virtual-env` before running commands, and
`abipy-run-tests` for focused validation. Use `abipy-flowtk` when a command inspects or mutates a flow.

## Match the command's existing framework

Most modules in `abipy/scripts` use `argparse`, a `get_parser(with_epilog=False)` function, and
`@abipy.tools.cli_parsers.prof_main`. `abiml.py` uses Click and its decorator helpers. Extend the framework already
used by the command; do not mix Click and argparse inside one command merely for a new subcommand.

For argparse commands:

- keep reusable parent parsers for options shared by several subcommands;
- use `dest="command"` dispatch consistently with neighboring branches;
- provide actionable positional and option help, defaults, units, and constrained `choices` where appropriate;
- include examples in the raw-text epilog when syntax is not obvious;
- configure logging through `cli_parsers.set_loglevel` or the established local equivalent;
- reuse `EnumAction`, `range_from_str`, Panel-serving, plotting, and profiling helpers instead of duplicating them.

For Click commands, follow the established option decorators, context usage, configuration callbacks, and
Click-native path/type validation.

All top-level tools should support useful `--help` and `--version` behavior. Preserve profiling through `prof`,
`tuna`, and `snakeviz` when the script uses `@prof_main`.

## Packaging and entry points

AbiPy currently packages every Python file under `abipy/scripts` through the `scripts=` mechanism in `setup.py`, not
through `[project.scripts]` console entry points. A new top-level command therefore needs a script file with the
project shebang/main convention and an executable installation path. Do not add a second packaging mechanism for a
single command without an intentional project-wide migration.

Prefer adding a subcommand to the tool that already owns the domain. Create a new top-level script only when its
purpose and option space are distinct enough to be discoverable independently.

## Dispatch and exit codes

Keep parsing, dispatch, and command implementation visibly separated. A command handler should return an integer:

- `0` for success;
- nonzero for an expected failure the caller can act on;
- argparse/Click usage errors for invalid syntax.

`@prof_main` ultimately calls `sys.exit(main())`, so returned values become process exit codes. Do not print an error
and then return success. Avoid catching broad exceptions solely to turn every defect into an uninformative message;
when handling an expected exception, preserve enough context for diagnosis and return nonzero.

Write primary command results to stdout and warnings/errors to stderr when practical. Keep ordinary output stable
enough for shell use, and do not add decorative text to a branch intended to emit JSON, CSV, YAML, or another
machine-readable format. Respect existing no-color/no-logo options and avoid terminal control codes when output is
redirected.

## Paths, files, and AbiPy objects

Use `abilab.abiopen` or the appropriate domain loader rather than reproducing suffix dispatch. Close file-backed
objects deterministically, preferably with a context manager. Use `abipy.data` only for bundled examples/tests, not
as an implicit substitute for a missing user path.

Resolve and validate input paths before expensive work. For output paths:

- distinguish printing to stdout from writing a file;
- state the selected format and destination;
- do not overwrite an existing file or directory silently unless overwrite behavior is already explicit in that
  command;
- create parent directories only when the command contract promises it;
- report partial output if a later operation fails.

Treat pickles as trusted-input formats. Do not imply that arbitrary pickle files are safe to open.

## Optional and interactive features

Import optional dependencies inside the branch that needs them so `--help`, `--version`, and unrelated subcommands
remain usable. On a missing optional package, identify both the requested feature and dependency; do not mask an
unrelated import or runtime error as an installation problem.

Plotting, IPython, Jupyter, Panel servers, clipboard operations, web browsers, and GUI viewers require interactive
resources. Keep noninteractive commands headless, honor options such as `--no-browser`, and configure Matplotlib
before importing plotting backends when the existing helper expects that order.

Network-backed identifiers or downloads must be explicit in the command and should produce actionable failures
when credentials or connectivity are missing. Do not introduce network access into a formerly local operation as an
implicit fallback.

## Mutating and external operations

Classify each subcommand before implementing it:

- inspection/formatting;
- local file creation or modification;
- external process execution;
- flow or scheduler mutation;
- destructive cleanup.

The latter categories need clear help text and explicit flags. Before deleting, overwriting, submitting, cancelling,
restarting, or modifying a persisted flow, resolve the exact target and use the confirmation/dry-run convention of
the parent command. A `--yes` or force option is authorization for the documented target only, not broader cleanup.

Do not prompt from commands intended for pipelines unless the user selected an interactive mode. In noninteractive
contexts, a missing required confirmation should fail rather than assume yes. Keep analysis commands read-only even
when an underlying object exposes mutating methods.

## Reuse shared behavior

Place a helper in `abipy/tools/cli_parsers.py` only when multiple commands genuinely share it. Keep domain-specific
formatting or dispatch near its command or, preferably, in the domain library if useful programmatically.

Common reusable facilities already include:

- log-level setup;
- Panel server options and keyword conversion;
- Matplotlib/seaborn and figure-exposure options;
- enum parsing and range parsing;
- OpenMP thread normalization;
- profiling of `main`.

Do not broaden a shared helper with command-specific assumptions that could change other scripts.

## Tests

Add subprocess-level coverage in `abipy/scripts/tests/test_scripts.py` or a focused neighboring test module. The
existing `ScriptTest` pattern uses `scripttest.TestFileEnvironment`, an isolated working directory, and a headless
Matplotlib backend.

At minimum, test:

- `--help` and `--version` for a new top-level tool;
- the new subcommand's smallest successful invocation;
- return code and relevant stdout/stderr content;
- invalid syntax or a representative expected failure;
- output-file creation and overwrite behavior when applicable;
- missing optional dependency behavior if the branch can be isolated reliably.

Use small files from `abipy.data` and `abipy/test_files`. Do not open browsers, GUIs, notebooks, or Panel servers in
unit tests. Mock external services and scheduler mutations; skip only when the test genuinely requires an optional
executable or dependency.

Also test the underlying library function directly in its domain test module. CLI tests should prove parsing,
dispatch, formatting, and process semantics, not duplicate an entire scientific test suite.

Run the script through the same path users exercise, not only by calling `main()` in process. Report the exact
invocation, exit code, important output, and whether the command created files, launched a process, accessed the
network, or mutated a flow.
