---
name: abipy-robots
description: Build, extend, and test AbiPy Robot classes for lifecycle-safe analysis, comparison, convergence, grouping, and plotting across collections of output files. Use for multi-file APIs; not for a single-file reader or generic plotting changes.
---

# Develop AbiPy Robots

A Robot owns a labeled collection of compatible AbiPy file objects and provides analysis that is meaningful across
files. Keep single-file decoding and physics on the file/domain object; add a Robot when comparison, aggregation, or
convergence is the actual abstraction. Use `abipy-file-reader` for a new format, `abipy-plotting` for figure
conventions, and `abipy-reference-data` plus `abipy-run-tests` for fixtures and validation.

## Start from the base contract

Subclass `abipy.abio.robots.Robot` and define the canonical `EXT` handled by the associated file class. Inspect a
mature Robot in the same scientific domain before designing dataframe columns, convergence methods, notebook cells,
or mixins such as `RobotWithEbands` and `RobotWithPhbands`.

Ensure `class_handles_filename` recognizes real produced basenames, especially formats with special names such as
`anaddb.nc` or exclusions such as `_DDB.nc`. A generally supported Robot must be discoverable through
`Robot.class_for_ext` and exported through the same `abilab` or package path as comparable classes.

Do not introduce a Robot merely to mirror one file's methods. A small collection with no meaningful cross-file
operation can use the base machinery until a stable comparison API emerges.

## Construct collections deliberately

Support the base constructors rather than inventing parallel discovery logic:

- `from_files` for explicit inputs and optional labels;
- `from_dir`, `from_dirs`, or `from_dir_glob` for extension-based scanning;
- `from_flow` or `from_work` for products already present in a workflow;
- JSON reconstruction only when path-based persistence is appropriate.

Labels are stable user-facing identifiers and dataframe indices, not necessarily filepaths. Require uniqueness,
preserve insertion order, and use `change_labels`, `remap_labels`, or `trim_paths` rather than mutating `_abifiles`
from client code. `abspath` controls labels; file objects still retain their actual paths.

When building from a Flow, honor `outdirs`, `nids`, `ext`, and `task_class` filtering and use node output-directory
helpers. Do not rerun tasks or infer missing products during collection construction.

## Preserve file ownership

Robots opened from path strings own those handles and must close them. File objects supplied by the caller are not
automatically owned. Maintain the base `_do_close` semantics when adding, filtering, removing, or popping entries.

Use Robots as context managers in scripts and tests:

```python
with abilab.GsrRobot.from_files(paths, labels=labels) as robot:
    frame = robot.get_dataframe()
```

If a candidate opened from a path is rejected by a filter or insertion fails, close the newly owned object. Do not
close caller-owned objects or reopen every file inside each analysis method. `remove()` deletes underlying files and
is materially different from `close()`; use it only when deletion is explicitly intended.

## Design comparison data

Provide a `get_dataframe` method when tabular comparison is useful. Use labels or normalized paths as a deterministic
index and make optional sections explicit, following nearby names such as `with_geo`, `with_params`, `with_spin`,
and `abspath`. Start from each file's semantic properties and `params`; do not reach through a low-level NetCDF
reader repeatedly when a high-level property exists.

Column names must identify quantities and units where ambiguity is possible. Preserve numerical dtypes, use
well-defined missing values for genuinely optional data, and keep row order aligned with Robot order. Accept custom
`funcs` only through the established `(name, value)` callback convention so failures are recorded in
`robot.exceptions` rather than corrupting unrelated columns.

Reuse `get_params_dataframe`, structure-dataframe helpers, `sortby`, and `group_and_sortby`. Their selectors may be
callables, dotted attributes, or keys in `abifile.params`; raise an informative error for an unresolved selector
instead of silently substituting labels or array indices.

## Implement convergence and plots

Put cross-file plots on the Robot and single-file plots on the file object. For convergence methods, use the base
sorting and hue grouping so callable selectors, dotted attributes, parameter keys, labels, and group order behave
consistently. Preserve input order when sorting is not requested and label axes with the selected quantity and units.

Separate dataframe or array construction from rendering when the transformation is substantial. Follow
`abipy-plotting` for `@add_fig_kwargs`, caller-supplied axes, `show=False` in nested calls, Plotly behavior, and
headless tests. Do not let plotting methods reopen files or mutate Robot labels.

`yield_figs` should expose a small, inexpensive set of useful automatic comparisons. Notebook generation should use
the base Robot code cells and reconstruct the intended labels/paths; it is optional unless nearby Robots support it.

## Serialization constraints

Base MSON serialization stores filepaths needed to reconstruct the Robot, not loaded scientific arrays or live file
handles. Keep this representation portable and deterministic. `trim_paths` changes labels for presentation; it does
not rewrite the stored file object's path.

When adding fields to serialization, verify a round trip and decide explicitly whether paths should be absolute or
relative to a documented base. Do not promise that a serialized Robot is self-contained when it still depends on
external output files.

## Tests

Put generic base behavior in `abipy/abio/tests/test_robots.py` and domain behavior in the owning package's tests. Use
small files from `abipy.data` and context managers. Cover the applicable invariants:

- filename handling and `class_for_ext` registration;
- construction from one file, several files, directories, and a Flow when relevant;
- stable unique labels, path trimming, sorting, and hue grouping;
- dataframe index, columns, units, row order, representative values, and optional branches;
- owned handles close while caller-owned file objects remain usable;
- JSON/MSON reconstruction when supported;
- convergence figures with display disabled;
- empty, duplicate-label, missing-selector, filtered-file, and heterogeneous-file behavior.

Do not delete reference outputs, scan large directory trees, start jobs, or require a graphical display. Test
scientific comparisons with meaningful invariants rather than only checking that a dataframe or figure is non-null.

Report the Robot class and `EXT`, supported construction paths, ownership behavior, dataframe/comparison APIs,
public registration, and focused tests used.
