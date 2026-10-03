---
name: abipy-inputs
description: Create, transform, validate, and test AbiPy AbinitInput and MultiDataset objects, including structures, pseudopotentials, k-point sampling, factories, decorators, DFPT inputs, and ABINIT variable semantics. Use for programmatic ABINIT input generation; not for merely parsing an existing .abi text file.
---

# Work with AbiPy inputs

Build inputs through AbiPy's domain objects and preserve the relationships among structure, pseudopotentials,
sampling, run level, tolerances, and file dependencies. Prefer an established factory or transformation from a
ground-state input when it represents the requested workflow.

Use `abipy-navigate` to find related implementations, `python-virtual-env` before executing Python, and
`abipy-run-tests` for test selection.

## Understand variables before setting them

Use AbiPy's variable database for names, defaults, dimensions, and documentation:

```bash
abidoc.py man ecut
abidoc.py apropos phonon
abidoc.py find eph
abidoc.py withdim natom
```

The same data is available through `abipy.abio.abivars_db`. For implementation details in the active ABINIT source,
run from its repository root:

```bash
./abisrc.py var VARIABLE
./abisrc.py var VARIABLE -v
```

`abisrc.py var` shows where a variable is declared, defaulted, read, checked, printed, documented, used, and set by
tests. Use it when behavior depends on the ABINIT version or the documentation does not explain a constraint.

Keep spell checking enabled. Do not disable it to accept an unknown variable unless the task explicitly targets a
new ABINIT variable that AbiPy's database does not contain yet; update the database through the established project
workflow instead of normalizing permanent spell-check bypasses.

## Choose the construction path

- Use `AbinitInput(structure, pseudos)` for a genuinely custom single calculation.
- Use a function in `abipy/abio/factories.py` for established SCF, bands, relaxation, GW, DFPT, phonon, conductivity,
  or related workflows. Inspect the factory signature and tests rather than duplicating its coupled defaults.
- Use transformations such as `*_from_gsinput` when a downstream calculation must inherit a validated ground-state
  setup. These functions handle details such as removing stale input-file-read variables and selecting compatible
  tolerances.
- Use `MultiDataset` for related datasets with the same pseudopotentials. Apply shared settings to the container,
  specialize individual entries, then call `split_datasets()` when registering tasks or returning independent inputs.

Do not mutate a caller-owned template unless the API promises in-place modification. Use `deepcopy`, `new_with_vars`,
`replicate`, or `MultiDataset.replicate_input` as appropriate.

## Structure and pseudopotential invariants

Pass any supported structure representation through `Structure`/`Structure.as_structure`. Do not set ABINIT geometry
variables such as `acell`, `rprim`, `xred`, `typat`, or `znucl` with `set_vars`; `AbinitInput` deliberately rejects
geometry variables. Use `set_structure` or construct a new input.

Let `PseudoTable` match and order pseudopotentials for the structure. Verify that every species has exactly one
compatible pseudo and do not depend on the order supplied by the caller. PAW and norm-conserving inputs have
different cutoff requirements; use pseudo hints or existing factory helpers rather than inventing `ecut` and
`pawecutdg` defaults.

Use `enforce_znucl` and `enforce_typat` only when compatibility with an external ABINIT file requires the original
type ordering. `typat` remains one-based in this interface and must match `natom`. Preserve these fields when copying,
serializing, or assembling a `MultiDataset`.

All entries in a `MultiDataset` must use the same pseudopotential set. Inputs may have different structures only
when their atomic types remain compatible with that set.

## Set coherent variables

Prefer semantic helpers over manually coordinating related variables:

- `set_kmesh`, `set_autokmesh`, `set_kpath`, and `set_qpath` for sampling;
- `set_cutoffs_for_accuracy` when pseudo hints are the intended source;
- `add_abiobjects` for objects exposing `to_abivars`;
- domain factories for perturbations and response-function datasets.

`set_vars` updates certain associated dimensions such as `nshiftk`, but do not assume it resolves every ABINIT
constraint. Assign arrays in the shapes expected by the helper or variable documentation.

An input should contain one appropriate SCF tolerance. Assigning a member of the SCF tolerance family replaces the
previous one. Choose the tolerance for the run level rather than layering `tolvrs`, `toldfe`, `tolwfr`, and similar
variables indiscriminately.

Use `set_vars_ifnotin` for defaults that must respect caller overrides. Use `pop_vars` for optional cleanup and
`remove_vars(..., strict=True)` when absence indicates a programming error. Preserve comments, tags, decorators,
and enforced type ordering across transformations when they remain semantically applicable.

## Executable-backed operations

Methods beginning with `abiget_` and `abivalidate` invoke ABINIT through a task and may create temporary work files.
They require a verified Python environment, ABINIT executable, and usable `TaskManager`. Do not present these as
pure in-memory checks.

Validate at the level justified by the change:

1. Inspect `to_string()` and object invariants without invoking ABINIT.
2. Test serialization or transformations when relevant.
3. Use `abivalidate()` for parser and cross-variable consistency when ABINIT is configured.
4. Run a real calculation only when scientific runtime behavior must be established.

On validation failure, inspect the returned log and stderr. Do not weaken variables or skip validation merely to get
a zero return code.

## Implementing or changing factories

Keep shared physical policy in the factory and return ordinary `AbinitInput`/`MultiDataset` objects. Reuse helper
functions for cutoffs, shifts, band counts, electronic settings, and stopping criteria. Make caller overrides
explicit in the signature and avoid silently overriding them later.

For a transformation from a ground-state input:

- copy the source input;
- remove variables that read files inappropriate for the new task;
- replace sampling, run-level, band-count, output, and tolerance variables coherently;
- preserve structure, pseudos, spin configuration, and intentional metadata;
- document assumptions about prerequisite files and task dependencies.

Do not encode Flow/Task dependency wiring inside an input factory; that belongs in `flowtk`.

## Tests

Add focused coverage in `abipy/abio/tests/test_inputs.py`, `test_factories.py`, or `test_decorators.py` according to
the changed layer. Use structures and pseudos from `abipy.data` rather than machine-specific paths.

Test observable invariants: returned type and dataset count, required variables and shapes, preservation of the
source input, pseudo/type ordering, serialization when supported, and expected failure for invalid combinations.
Use `AbipyTest.abivalidate_input` or `abivalidate_multi` only for tests that are already executable-backed and skip
cleanly when their ABINIT requirement is unavailable.

Report whether validation was in-memory, dry-run through ABINIT, or a real calculation, together with the interpreter
and ABINIT version when they affected the result.
