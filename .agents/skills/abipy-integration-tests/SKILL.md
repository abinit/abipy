---
name: abipy-integration-tests
description: Add, run, and diagnose AbiPy integration tests that execute ABINIT flows through TaskManager and PyFlowScheduler, including MPI parameterization, timeouts, cleanup, restarts, handlers, and scientific output checks. Use for abipy/integration_tests; use abipy-run-tests for ordinary unit tests.
---

# Develop AbiPy integration tests

Use an integration test only when correctness depends on a real ABINIT executable, task filesystem products,
dependency transfer, scheduler lifecycle, restart behavior, event handling, or multi-step scientific output. Keep
pure input, graph, reader, and transformation behavior in the faster domain unit tests.

Use `abipy-inputs` for input semantics, `abipy-flowtk` for workflow lifecycle, `python-virtual-env` before executing
pytest, and `abipy-reference-data` if the test requires bundled input or comparison data.

## Understand the harness

Tests live in `abipy/integration_tests` and are collected by its `pytest.ini`:

- files use `itest_*.py`;
- test functions use the `itest_` prefix;
- temporary results use `_integration_tests_` as the pytest base temporary directory;
- the default per-test timeout is 180 seconds when `pytest-timeout` is installed.

The shared `fwp` fixture provides:

- a pytest-managed temporary `workdir`;
- a `TaskManager` derived from `~/.abinit/abipy/manager.yml`;
- a configured `PyFlowScheduler` from the user scheduler configuration;
- `AbinitBuild` information for capability/version checks.

The `tvars` fixture parameterizes ABINIT variables such as `paral_kgb`. Manager configurations are parameterized as
well: fixed MPI execution is the default, controlled by `ABIPY_ITEST_MPI_PROCS` (default 2), while
`ABIPY_ITEST_AUTOPARAL=1` adds an autoparal configuration.

Use these fixtures rather than loading personal configuration or creating arbitrary work directories inside each
test. The autouse guard also:

- applies scheduler exception limits;
- adds ABINIT `autoparal 1` in the narrow supported case when AbiPy autoparal is disabled;
- tracks ShellAdapter processes and kills surviving process groups after failures/timeouts;
- appends task error-file tails and the flow workdir to failed pytest reports.

Do not bypass this harness without a concrete lifecycle requirement.

## Preflight before running

Verify the selected Python environment, working-tree import, pytest and `pytest-timeout`, the ABINIT executable, and
the manager/scheduler files. Inspect the active manager's queue-adapter type and resource limits before starting.

The suite can execute through the user's configured queue adapter. A command that reaches Slurm, PBS, or another
shared scheduler can submit and cancel real jobs; obtain explicit authorization for that external mutation. Do not
assume “run tests” means queue submission. Prefer a deliberately configured shell manager for local development
when it represents the behavior under test.

Confirm the intended MPI process count and ABINIT build capabilities. Do not silently replace the user's manager,
MPI launcher, executable path, or modules merely to make a failing environment run.

## Add a focused test

Build the smallest scientifically meaningful system and workflow. Use structures and pseudopotentials from
`abipy.data`, low but valid cutoffs and meshes, relaxed tolerances appropriate to the assertion, and as few tasks as
the behavior requires. Label deliberately non-production parameters in comments.

Accept `fwp` whenever the test creates a flow and `tvars` when the behavior should be checked across the suite's
parallel-variable matrix:

```python
def itest_feature(fwp, tvars):
    inputs = make_inputs(tvars)
    flow = build_flow(fwp.workdir, inputs, manager=fwp.manager)
    scheduler = run_flow(flow)
    assert scientific_invariant(flow)
```

Use `abipy.integration_tests.helpers.run_flow` for the normal scheduler-and-check path. Use `assert_flow_ok` when the
test needs custom scheduler steps or direct task execution. These helpers check `flow.all_ok`, scheduler exceptions,
and work finalization, and produce useful failure reports.

Do not duplicate scheduler polling loops or merely assert a zero scheduler return code. `Done` is not success;
verify `Completed`/`flow.all_ok` and finalization unless the test intentionally targets incomplete behavior.

## What to assert

Test the contract that requires integration, not every intermediate implementation detail. Depending on the case,
assert:

- graph dependencies and required products before execution;
- task states and bounded launch/restart/correction counts;
- existence and readability of essential output files;
- representative dimensions, units, values, conservation laws, or symmetry invariants;
- correct work finalization and post-processing outputs;
- event-handler corrections and preserved history for recovery tests;
- scheduler termination under deadlock or no-runnable-task conditions.

Use tolerances justified by the numerical setup and stable across supported platforms. Do not make a test pass by
loosening a tolerance without examining whether the difference is physical, parallel-order noise, or a regression.
Avoid exact text comparisons of ABINIT logs unless the text is itself the compatibility contract.

Assert exact `scheduler.nlaunch` only when the graph and restart policy make the count deterministic across the
parameterized managers. Otherwise assert the task/result behavior.

## Restarts, failures, and handlers

For a restart test, deliberately create a bounded and understood first failure, verify its status and restart
products, apply the documented correction, restart through the public API, and assert both the restart count and
final scientific result. Never use an unbounded retry to wait for an intermittent outcome.

For event-handler tests, trigger a specific ABINIT event with a minimal input and assert the correction's event type
and effect. If the triggering behavior changed in a newer ABINIT version, use a version/capability condition or remove
the obsolete test; do not leave unconditional `xfail` as permanent coverage.

Scheduler failure tests may use `abipy.flowtk.mocks` when the target is scheduler graph logic rather than a real
ABINIT failure. This keeps the failure deterministic while still exercising scheduler termination and deadlock
detection.

## Skips, xfails, and portability

Prefer skip conditions based on the actual missing capability: ABINIT version/build option, optional executable,
library, MPI feature, or scheduler support. Hostname lists are a last resort and must explain the concrete platform
limitation. Do not mark a scientifically failing test skipped merely because it is slow or inconvenient.

Use `xfail` only for a known, tracked incompatibility whose failure is still informative. Give a precise reason and
remove the marker when the underlying behavior is fixed. An unconditional xfail is not evidence that a handler or
workflow remains functional.

Override the timeout with `@pytest.mark.timeout(seconds)` only when the minimal valid calculation demonstrably needs
more than the suite default. Keep scheduler `frozen_timeout` below the pytest timeout so a stalled task is diagnosed
before pytest kills the test process.

## Run narrowly and preserve diagnostics

From the repository or integration-test directory, begin with one node and display parameter IDs:

```bash
python -m pytest abipy/integration_tests/itest_file.py::itest_feature -vv -s
```

Then expand to the affected integration module. Enable `ABIPY_ITEST_AUTOPARAL=1` only when the change concerns
autoparal or after the fixed-MPI case passes; it multiplies runtime and may exercise different decompositions.

Do not run the entire integration suite as the first check. It launches many ABINIT calculations and the ordinary
top-level `invoke pytest` intentionally excludes it.

On failure, preserve the reported flow workdir until the cause is understood. Start with the pytest report sections,
then inspect task status/history and the relevant `run.err`, `queue.qerr`, `__ABI_MPIABORTFILE__`, `run.log`, and
`run.abo`. Separate:

- invalid input or ABINIT scientific failure;
- missing dependency product;
- queue/launcher/environment failure;
- timeout or frozen output;
- parser/reader assertion failure after a successful calculation;
- nondeterminism across MPI decompositions.

Do not rerun repeatedly without changing the diagnosis. The teardown guard kills tracked local shell processes, but
after a timeout or harness failure also verify that no authorized external scheduler jobs remain.

## Cleanup and completion

Keep all produced files inside `fwp.workdir` or a pytest temporary directory. Avoid manual fixed paths and do not
write new outputs into `abipy/data`, the repository root, or a user's flow directory. Let pytest manage its base temp
tree; retain a failed workdir only long enough for diagnosis.

Report the exact pytest node, parameter IDs, interpreter, ABINIT version/build, manager adapter and MPI settings,
whether external jobs were submitted, elapsed time, pass/fail/skip summary, final flow status, and preserved failure
workdir when applicable.
