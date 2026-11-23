---
name: abipy-flowtk
description: Build, extend, inspect, persist, test, and troubleshoot AbiPy Flow/Work/Task workflows, including file dependencies, managers, schedulers, callbacks, task states, restarts, and event handling. Use for execution orchestration; use abipy-inputs for calculation input semantics alone.
---

# Work with AbiPy flowtk

Model execution as a graph of `Flow`, `Work`, and `Task` nodes. Inputs describe individual calculations;
`flowtk` owns their directories, file dependencies, execution policy, state transitions, recovery, and persistence.

Use `abipy-inputs` for constructing coherent inputs, `abipy-navigate` to find domain workflows, and
`python-virtual-env` plus `abipy-run-tests` before executing tests or workflows.

## Choose the highest existing abstraction

Search `abipy/flowtk/*_flows.py`, `*_works.py`, `flows.py`, and the examples before assembling nodes manually.

- Use an existing `Flow` constructor or specialized `Work` when it already represents the scientific workflow.
- Use `Flow.from_inputs` only for independent tasks; its API explicitly does not represent dependencies among them.
- Assemble a generic `Work` for a small, static group of related tasks.
- Add a specialized `Work` when coordination, post-processing, or `on_all_ok` behavior belongs to that group.
- Add a specialized `Task` only when executable invocation, output products, restart behavior, or event handling differs
  from existing task classes.
- Use a dynamic callback only when later topology genuinely depends on completed results and cannot be known when the
  flow is built.

Keep physical input policy in factories and workflow topology in `flowtk`. Do not make an input factory register
tasks or encode filesystem dependencies.

## Construct the graph

The usual static pattern is:

```python
flow = flowtk.Flow(workdir=workdir, manager=manager)
work = flowtk.Work()
scf_task = work.register_scf_task(scf_input)
nscf_task = work.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
flow.register_work(work)
flow.allocate()
```

Prefer typed registration methods such as `register_scf_task`, `register_nscf_task`, `register_phonon_task`, and
`register_eph_task`; they document intent and select the appropriate task class.

Dependencies map a producer node to one or more ABINIT products:

```python
deps={producer_task: "DEN WFK"}
```

Use only registered ABINIT extensions and request the minimum products required by the consumer. `Product` maps the
extension to the connecting input variables and expected output path. Do not pass a path as an informal dependency
when a producer node exists. File paths are appropriate for deliberate external inputs.

Use dependency getters such as `@structure` only for their implemented semantics. A dependency may combine products
and getters, but the consumer must remain scientifically valid after the getter changes its input.

Dependencies can be attached to a task or an entire work. Check parent/child relationships after `allocate()` and
run `flow.check_dependencies()` to detect invalid or cyclic topology.

## Allocation, building, and persistence

Registration creates the logical graph. `flow.allocate()` propagates work directories, managers, node positions,
and relationships. `flow.build()` creates runtime files and directories. `flow.build_and_pickle_dump()` builds and
persists the graph in `__AbinitFlow__.pickle`.

Do not treat these operations as interchangeable, and do not allocate repeatedly with a different work directory.
Use an explicit temporary directory in tests.

Flow state is persisted in its pickle database. For read-only inspection, load with the default spectator mode:

```python
flow = flowtk.Flow.pickle_load(workdir)
```

Spectator mode avoids signal-driven mutations. Removing a pickle lock is an exceptional recovery operation; first
establish that no scheduler or writer is active. Never delete or replace an existing flow work directory merely to
resolve a construction error.

Options such as `remove=True` and `flow_main --remove` delete the existing work directory. Use them only when the
user explicitly intends to discard that flow and after resolving the exact path.

## Managers and execution

`TaskManager` selects queue adapters, resource limits, launch commands, environment setup, and parallel policy. An
implicit manager is loaded first from `manager.yml` in the current directory and then from
`~/.abinit/abipy/manager.yml`. For reproducible tests and scripts, pass or construct the manager explicitly.

Inspect manager documentation with:

```bash
abidoc.py manager
abidoc.py manager slurm
```

Do not silently rewrite a user's manager configuration to make a workflow pass. Resource limits can be specialized
by task class; preserve that behavior when copying or replacing managers. Use `to_shell_manager` only when local
shell execution is intended and compatible with the executable environment.

The following operations can execute or submit jobs and are not mere graph validation:

- `flow.make_scheduler().start()`;
- `flow.rapidfire()` and `flow.single_shot()`;
- `PyLauncher` methods;
- task `start`, `restart`, `cancel`, and queue operations;
- scripts decorated with `flow_main` when invoked with `--scheduler`.

Obtain the required authorization before submitting, cancelling, or restarting external jobs. Use
`flow.abivalidate_inputs()` or `flow_main --abivalidate` when only ABINIT input validation is needed.

## Status and recovery

Interpret states precisely: `Done` means execution ended, not that results are correct; `Completed` is success.
`AbiCritical`, `QCritical`, `Unconverged`, and `Error` require different diagnosis. Call `check_status()` and inspect
task history, events, output, log, stderr, queue diagnostics, and corrections before proposing recovery.

Let task and event-handler APIs perform state transitions. Do not force status values or edit persisted state to
hide a failure. A restart must preserve the scientifically required products and update the dependency graph or
input only through established restart/correction mechanisms.

When adding recovery logic:

- distinguish ABINIT events from queue/submission failures;
- make corrections bounded and record them in node history;
- avoid retry loops without a clear stopping condition;
- confirm restart files exist before advertising restartability;
- preserve idempotence when status checks run repeatedly.

## Specialized works and callbacks

A specialized `Work` should register its internal tasks using the normal dependency API and keep group-level
finalization in lifecycle hooks such as `on_all_ok`. Hook code may be called after reload, so avoid assumptions based
only on transient process state.

For dynamic flows, register callbacks through the Flow API with serializable callback data and explicit dependency
products. A callback must be safe to evaluate after persistence/reload and must not create duplicate works if status
checking occurs more than once. Prefer a static graph whenever all tasks are knowable up front.

Expose new generally useful task/work/flow classes through `abipy.flowtk` following the existing import pattern.
Avoid expanding the public API for a one-off example.

## Flow scripts

For reusable example scripts, follow the established `@flowtk.flow_main` pattern and return a fully constructed
flow from `main(options)`. Respect the supplied `options.workdir` and `options.manager`; do not substitute personal
paths or configurations. Building a flow should not submit it unless the caller selects scheduler execution.

Keep examples small enough to communicate topology and use `abipy.data` structures and pseudos where appropriate.

## Tests

Put tests beside the affected layer in `abipy/flowtk/tests`. Separate three levels:

1. Pure graph tests with a temporary work directory and an explicit manager string/object.
2. Lifecycle tests using `abipy.flowtk.mocks` to simulate starts, statuses, and results.
3. Executable-backed integration tests only when real ABINIT or scheduler behavior is essential.

For graph tests, verify node/task classes, task counts, positions, dependencies and requested products, manager
propagation, `check_dependencies`, build output, and pickle reload. Use `Flow.temporary_flow()` where suitable and
clean up test-local state such as the configured user manager.

Do not submit real queue jobs from unit tests. Mock status changes through the provided helpers rather than editing
private fields when a helper exists. For executable-backed tests, state the ABINIT, manager, and scheduler
requirements and skip cleanly when they are unavailable.

## Completion information

Report the graph shape, dependency products, manager source, whether the flow was only allocated/built or actually
executed, the final task-status summary, and the focused tests used. If jobs were submitted or cancelled, report the
affected flow and scheduler identifiers.
