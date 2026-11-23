---
name: abipy-run-tests
description: Select and run focused AbiPy pytest tests, distinguish unit from ABINIT-backed integration tests, and use AbiPy's test helpers and reference data. Use when validating AbiPy changes or diagnosing test failures.
---

# Run AbiPy tests

Activate and verify the project Python environment first. Run from the repository root and prefer
`python -m pytest` so the selected interpreter supplies pytest.

## Start focused

Run the closest test module or node first:

```bash
python -m pytest abipy/electrons/tests/test_gsr.py -q
python -m pytest abipy/electrons/tests/test_gsr.py::GsrFileTest::test_gsr -q
```

Then expand to the package tests affected by the change:

```bash
python -m pytest abipy/electrons/tests -q
```

Use `-x` while iterating, `-vv` for detailed collection/output, `-s` only when captured output hides relevant
diagnostics, and `--lf` only when the previous pytest cache is relevant to the current tree.

The repository-wide Invoke task runs coverage, doctests, and parallel workers while excluding integration tests,
reference data, scripts, examples, flows, and GUI code:

```bash
invoke pytest
```

Use it as a broader check, not as the first diagnostic command.

## Test categories

- Unit and file-reader tests live beside modules in `abipy/**/tests` and often derive from
  `abipy.core.testing.AbipyTest`.
- CLI tests live in `abipy/scripts/tests` and may exercise subprocess behavior.
- `abipy/integration_tests` can launch ABINIT workflows and depends on executables, managers, pseudopotentials, and
  machine configuration. Run these only when the change requires them and the environment is configured.
- Examples and bundled reference files are not automatically covered by the top-level Invoke task.

Inspect `abipy/core/testing.py` before recreating helpers. It provides dependency checks, temporary directories,
numerical assertions, and access patterns for data under `abipy/data`. Prefer `abipy.data.ref_file`, `cif_file`,
`pseudos`, and related helpers over paths tied to the checkout layout.

## Failure triage

Confirm that `abipy.__file__` points to this working tree. Separate missing optional dependencies, unavailable ABINIT
executables, absent external test data, and genuine assertion failures. Do not weaken tolerances, update binary
reference files, or turn failures into skips without establishing that the expected scientific behavior changed.

Report the exact test command, environment/interpreter, and pass/fail/skip summary.
