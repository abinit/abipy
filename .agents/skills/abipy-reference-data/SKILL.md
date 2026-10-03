---
name: abipy-reference-data
description: Find, add, regenerate, package, and validate AbiPy structures, pseudopotentials, scientific reference outputs, and test fixtures under abipy/data. Use when tests, examples, or readers need bundled data; avoid for transient outputs or large datasets better stored externally.
---

# Manage AbiPy reference data

Treat bundled data as maintained scientific test inputs, not as a convenient dump of calculation outputs. Reuse an
existing fixture whenever it exercises the required schema and values; add data only when it enables meaningful
coverage that cannot be obtained cheaply with a small synthetic object.

Use `abipy-file-reader` when the fixture supports a new file format, `netcdf` to inspect NetCDF content, and
`abipy-run-tests` when validating consumers.

## Find existing data first

Search both filenames and consumers:

```bash
find abipy/data -type f | sort
rg -n 'BASENAME|ref_file\(|cif_file\(|pseudo\(' abipy
```

Use the public helpers rather than checkout-relative paths:

- `abipy.data.cif_file` / `cif_files` for structures in `abipy/data/cifs`;
- `abipy.data.pseudo` / `pseudos` for bundled pseudopotentials;
- `abipy.data.ref_file` / `ref_files` for scientific outputs and other data-relative paths;
- `abipy.data.structure_from_ucell` for built-in unit cells;
- `abipy.data.pyscript` for examples.

`ref_file` indexes every `.nc` below `abipy/data/refs` by basename. NetCDF basenames should therefore be unique;
duplicates produce an ambiguous-data warning. For non-NetCDF files and explicit nested paths, pass the path relative
to `abipy/data`, for example `refs/case/file_DDB`.

Do not expose a new helper or module-level file list unless multiple consumers need a stable semantic API. A test
that needs one fixture can call `ref_file` directly.

## Choose the location

- Put reusable crystal structures in `abipy/data/cifs`.
- Put small redistributable pseudopotentials in the established pseudo directories only after checking format,
  species, provenance, and redistribution terms.
- Put related ABINIT outputs in a descriptive subdirectory of `abipy/data/refs` when several files form one
  calculation or comparison set.
- Put a single broadly reused NetCDF output at the appropriate existing location and keep its basename unique.
- Use `abipy/test_files` for CLI/parser fixtures that are not part of the public `abipy.data` scientific API.
- Keep transient work directories, scheduler artifacts, caches, plots, and generated scratch files out of the data
  tree.

Match the nearest existing dataset layout. Do not create a new directory hierarchy for one tiny file when an
established related reference set exists.

## Minimize and document scientific content

Before adding a binary file, inspect its size and contents:

```bash
du -h /path/to/file
ncdump -h /path/to/file.nc
```

Retain only variables and companion files required to reproduce or understand the tested behavior, but do not
rewrite an ABINIT NetCDF file with an ad hoc tool if doing so would make it cease to represent real ABINIT output.
Prefer a smaller physical calculation over post-hoc binary surgery.

For generated scientific outputs, preserve enough provenance in an existing input, log, generator script, or nearby
metadata convention to identify:

- producing program and version/commit;
- input variables and pseudopotentials;
- why this case is needed and which behavior it exercises;
- any intentional corruption, truncation, legacy schema, or post-processing.

Do not replace a legacy fixture merely because a new ABINIT version writes different metadata; compatibility tests
may depend on the old schema. Add or regenerate data only after deciding whether the expected behavior is backward
compatibility or adoption of the new format.

Avoid committing large outputs casually. If a fixture is several megabytes, estimate its repository and package
cost, look for a smaller calculation, and confirm that the test cannot isolate the behavior with mocks or arrays.
Datasets used for benchmarks, tutorials, or broad research analysis may belong in an external versioned dataset
rather than the Python package.

## Generate or update fixtures

Use the established `FilesGenerator`, `AbinitFilesGenerator`, or `AnaddbFilesGenerator` pattern when it fits. These
operations execute scientific codes and clean/rename outputs; run them only in a dedicated work directory with the
intended executables and environment.

Never generate directly inside `abipy/data/refs`. Generate in a temporary or dedicated calculation directory,
inspect results, then copy only the resolved files into the repository. Do not overwrite an existing fixture until
its consumers and compatibility purpose have been identified.

When changing a fixture, compare old and new schemas and representative values. Review changes in dimensions,
units, indexing, fill values, attributes, and numerical tolerances rather than accepting a new binary solely because
the producing calculation completed.

## Packaging

Being present in the source tree does not guarantee inclusion in both sdist and wheel. `MANIFEST.in` includes only
selected extensions recursively, while `setup.py` has explicit `package_data` patterns and a finite list of nested
reference directories.

When adding a new file type or nested reference directory:

1. update the narrowest appropriate packaging rule;
2. avoid broad wildcards that pull scratch files or large unrelated data into distributions;
3. build both sdist and wheel when packaging behavior changed;
4. inspect the archives and verify the installed package can resolve the fixture through `abipy.data`.

Do not rely on an editable checkout test as proof that packaged data is available.

## Tests and review

Add or update the nearest scientific consumer test. Add a focused assertion in `abipy/data/tests/test_data.py` only
when the public data helper itself changed.

Verify:

- the helper returns an existing regular file;
- the basename is unambiguous where basename lookup is used;
- the relevant reader opens and closes it;
- a small set of stable schema/value invariants proves it is the intended fixture;
- all known consumers still pass;
- packaging includes it when the file is part of distributions.

Avoid checksums as the only scientific test: they detect any byte change but do not explain correctness. Checksums
can supplement schema and value assertions when exact identity is required.

Report files added or replaced, byte sizes, provenance, producing-code version, affected consumers, package archive
verification, and whether an old fixture remains available for compatibility coverage.
