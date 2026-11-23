---
name: abipy-navigate
description: Locate AbiPy file readers, domain objects, flows, command-line entry points, reference data, and tests. Use when tracing behavior or deciding where an AbiPy change belongs.
---

# Navigate AbiPy

Use `rg` and the maps below to trace behavior from public entry points to implementation and tests. Activate the
project environment first when importing modules or running Python.

## Main entry points

- `abipy/abilab.py`: public convenience imports and `abiopen` file-extension dispatch. Start here when determining
  which class opens an ABINIT output file.
- `abipy/<domain>/`: scientific objects and file readers. Major domains include `electrons`, `eph`, `dfpt`, `waves`,
  `dynamics`, `wannier90`, and `lumi`.
- `abipy/iotools/`: common readers, NetCDF helpers, visualizers, and serialization utilities.
- `abipy/abio/`: ABINIT inputs, outputs, variable database, factories, and timers.
- `abipy/flowtk/`: flows, works, tasks, managers, schedulers, events, wrappers, and ABINIT execution machinery.
- `abipy/scripts/`: command implementations. Their tests are in `abipy/scripts/tests`.
- `abipy/panels/`: Panel-based user interfaces.
- `abipy/data/`: bundled structures, pseudopotentials, and reference outputs exposed through `abipy.data` helpers.
- `abipy/core/testing.py`: `AbipyTest`, dependency checks, numerical assertions, and test helpers.

## Search recipes

Find a public symbol, its implementation, and tests:

```bash
rg -n 'SymbolName' abipy/abilab.py abipy
rg -n 'class SymbolName|def symbol_name' abipy
rg -n 'SymbolName|symbol_name' abipy --glob '**/tests/**'
```

Find a file reader or NetCDF variable:

```bash
rg -n 'EXT\.nc|class .*File|abiext2ncfile' abipy/abilab.py abipy
rg -n 'variable_name|read_value|read_variable' abipy
```

Find CLI wiring and behavior:

```bash
rg -n 'def main|argparse|click|command' abipy/scripts abipy/flowtk
rg -n 'COMMAND_NAME' abipy/scripts/tests
```

Find bundled data through its API before hard-coding a path:

```bash
rg -n 'def (ref_file|cif_file|pseudos|structure_from_ucell)' abipy/data/__init__.py
rg -n 'BASENAME' abipy/data abipy
```

## Change placement

Keep domain behavior in its domain module, shared file-format mechanics in `iotools`, execution orchestration in
`flowtk`, and convenience exposure in `abilab.py`. Match an existing sibling reader or object before introducing a
new abstraction. Put unit tests in the nearest `tests/` directory and reuse `abipy.data` reference files when
possible.

Avoid navigating generated caches, build outputs, `_integration_tests_`, or submodule contents unless the task
specifically concerns them. Preserve lazy/optional dependency behavior: inspect nearby import patterns before adding
top-level imports for optional packages.
