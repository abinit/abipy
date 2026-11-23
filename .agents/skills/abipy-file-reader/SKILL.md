---
name: abipy-file-reader
description: Add or extend AbiPy readers and high-level objects for ABINIT output files, using abisrc.py to derive the format from the Fortran writer and integrating the result with abiopen, reference data, tests, and optional Robot or plotting APIs.
---

# Implement an AbiPy file reader

Treat the ABINIT writer and a representative output file as the format specification. Trace the complete path from
the Fortran definitions and writes to AbiPy's low-level reader, semantic file object, public registration, and tests.

Use the `netcdf` skill for detailed axis-order or value discrepancies, `abipy-navigate` to locate related Python
code, and `abipy-run-tests` when executing tests.

## Establish the schema from ABINIT

Locate the active ABINIT checkout rather than assuming the sibling repository is the right version. When it is
`../abinit`, run `abisrc.py` from that repository root so its source paths and caches are correct.

For a NetCDF format, start with:

```bash
cd ../abinit
./abisrc.py nc_explain
./abisrc.py nc_explain WRITER_PROCEDURE --follow 2
./abisrc.py nc_explain WRITER_PROCEDURE --follow 2 --check /path/to/example_FILE.nc
```

Without a procedure name, `nc_explain` lists routines that define NetCDF content. For a selected writer it reports
dimensions, variables, Fortran and NetCDF dimension order, types, groups, conditions, attributes, units, write
slices (`start`/`count`), fill values, and nested writer calls. `--check` compares the inferred schema with an actual
file; a mismatch can mean the file came from a different ABINIT version.

Use the other structural queries when names or semantics remain unclear:

```bash
./abisrc.py where NAME
./abisrc.py source NAME --decls --doc
./abisrc.py context PROCEDURE --no-tests
./abisrc.py callers PROCEDURE -d 2
./abisrc.py callees PROCEDURE -d 2
./abisrc.py callsites PROCEDURE
```

Prefer these focused queries to reading whole Fortran files. Fall back to `rg -in` for dynamically constructed
names, macros, or relationships the parser cannot resolve. Fortran is case-insensitive.

Before implementing, record:

- file suffix and the routine that creates it;
- dimensions and their physical meanings;
- variable names, types, units, groups, and conditional presence;
- NetCDF order as observed by Python;
- partial writes and fill-value semantics;
- one-based indices that must become zero-based in Python;
- complex-number representation and any version-dependent names.

Do not infer a transpose from reversed dimension listings alone. The NetCDF Fortran interface commonly exposes
Fortran `x(n1,n2,n3)` as `(n3,n2,n1)` to Python. Check the exact write call and slice.

## Match an existing AbiPy reader

Choose a nearby reader in the same scientific domain before designing a new abstraction. Typical layers are:

1. A specialized reader derived from `ETSF_Reader`, `ElectronsReader`, `BaseEphReader`, or another domain reader.
2. A high-level file object derived from `AbinitNcFile` plus applicable mixins such as `Has_Header`,
   `Has_Structure`, `Has_ElectronBands`, and `NotebookWriter`.
3. Optionally, a `Robot` for meaningful comparison or convergence analysis across multiple files.

The file object should own and close its reader. Existing code may expect both names, so follow the closest mature
implementation when assigning `self.r` and `self.reader`. Use `cached_property` for expensive immutable data, but
do not eagerly load large arrays merely to construct the object.

Keep low-level schema decoding in the reader and physical interpretation in well-named methods or domain objects.
Make unit conversions explicit with AbiPy or pymatgen unit helpers. Document array shapes in Python order and state
whether indices are zero- or one-based. Preserve masked/fill values when they carry meaning.

For compatibility across ABINIT versions, distinguish required variables from optional or renamed ones. Check
`rootgrp.variables`, groups, or use a supported `default` argument only when absence has a defined meaning. Do not
hide malformed required data behind broad fallbacks.

## Public integration

For a generally supported file format:

- import the file class in `abipy/abilab.py`;
- register its suffix in `abiext2ncfile` for NetCDF or `ext2file` for other formats;
- expose intentionally public classes through the module's established `__all__`/package pattern;
- verify `abilab.abiopen(path)` returns the new class;
- implement informative `to_string`/`__str__` output and deterministic `close` behavior.

Filename matching order matters when suffixes overlap. Follow the existing registry convention and test the actual
basename used by ABINIT.

Add these only when they provide real value:

- `Robot` with an `EXT`, dataframe, or comparison API for collections of files;
- `plot_*` methods using AbiPy plotting helpers and `@add_fig_kwargs`;
- `yield_figs` for automatic `abiopen -e` output;
- `write_notebook` or Panel support for interactive exploration.

Do not require every format to implement the entire presentation stack.

## Tests and reference data

Put tests in the nearest domain `tests/` directory. Prefer an existing file from `abipy/data/refs`, accessed via
`abipy.data.ref_file`, over a hard-coded repository path. Before adding a binary fixture, check its size, provenance,
and whether an existing file exercises the same schema.

Cover the behavior that proves the integration:

- construction and context-manager closing;
- representative dimensions, scalars, array shapes, units, and index conversion;
- conditional or version-dependent variables when supported;
- derived scientific quantities with known values or invariants;
- `abiopen` dispatch for a newly registered suffix;
- Robot or plotting behavior only when those features were added.

Use a real reference file for schema integration and small synthetic objects or mocks only for isolated transforms
and edge cases. Do not update numerical expectations or binary references merely to silence a failure.

## Completion check

Run the narrow test module first, then the affected package tests. Report the ABINIT source/version used to infer the
schema, the reference file inspected, supported compatibility behavior, and any schema fields intentionally left
unexposed.
