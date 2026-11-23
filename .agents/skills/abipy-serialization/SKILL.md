---
name: abipy-serialization
description: Add, evolve, and test AbiPy object serialization with Monty/MSON JSON, pickle helpers, slot state, Flow persistence, and generated notebooks. Use for persistence and round-trip compatibility; not for NetCDF scientific file readers.
---

# Develop AbiPy serialization

Choose the representation from the persistence contract, not convenience. JSON/MSON is inspectable and suitable for
portable object state; pickle preserves richer Python graphs but is version-sensitive and unsafe for untrusted
input; notebooks are executable reconstructions, not serialized data. Use `netcdf` or `abipy-file-reader` for ABINIT
scientific file formats and `abipy-flowtk` for workflow semantics beyond persistence.

## Select the established mechanism

Inspect `abipy/tools/serialization.py`, `abipy/core/mixins.py`, and a nearby object with the same lifetime.

- Use Monty's `MontyEncoder`/`MontyDecoder` and `mjson_write`, `mjson_load`, or `mjson_loads` for MSONable JSON.
- Derive small result containers needing both JSON and pickle file helpers from `Serializable`.
- Use `HasPickleIO` when the established API persists an object under a work directory and conventional basename.
- Use `SlotPickleMixin` only for classes with `__slots__` whose slot values form their complete pickle state.
- Use the specialized `Flow.pickle_*` implementation for Flow graphs; do not replace it with a generic mixin.
- Use `NotebookWriter` only when an object can generate a useful, reproducible analysis notebook.

Do not create another wrapper when one of these contracts fits. Avoid mixing YAML configuration, NetCDF output, and
object persistence merely because all are file I/O.

## Implement MSON round trips

An MSONable object provides `as_dict` and a compatible `from_dict`. Decorate `as_dict` with `@pmg_serialize` when the
method does not otherwise add `@module` and `@class`; the decoder needs those exact importable identifiers.

Serialize constructor-level semantic state rather than caches, readers, open handles, callbacks, loggers, or derived
arrays that can be rebuilt cheaply. Use nested objects' own `as_dict` representation and let Monty recursively encode
NumPy, datetime, pathlib, and other supported values. Do not manually stringify a rich object unless the string is
the actual stable contract.

`from_dict` must tolerate Monty having already decoded nested objects. Reconstruct through public constructors where
possible so validation and normalization still occur. Do not mutate the caller's dictionary while popping metadata
or compatibility fields; copy it first.

Treat the serialized keys as a compatibility surface. When renaming or adding state:

- give new optional fields a documented default when old files have an unambiguous meaning;
- accept legacy keys narrowly and normalize them in one place;
- reject missing required or malformed state with a useful error;
- never silently substitute scientifically different defaults.

Include an explicit schema/version field only when real migrations require it. Do not add version machinery without
a compatibility policy and tests for an older representation.

## Use pickle with the right safety model

Never load a pickle from an untrusted or unexplained source; unpickling can execute arbitrary code. JSON is not a
drop-in safe format either, but ordinary parsing does not carry pickle's execution semantics.

Pickles generally require importable class paths and compatible Python/library definitions. Do not advertise them as
long-term or cross-version interchange. Exclude or reconstruct resources that cannot be pickled reliably, such as
open NetCDF handles, subprocesses, locks, weak references, and live scheduler connections.

For `__slots__`, ensure `SlotPickleMixin.__getstate__` covers every slot that determines behavior and that
`__setstate__` can restore a valid instance. If inheritance contributes slots or extra state, implement the combined
state deliberately rather than assuming the mixin discovers it.

Flow persistence is transactional and concurrency-sensitive. Preserve `FileLock`, `AtomicFile`, the configured
pickle protocol, `PmgPickler` handling of external pymatgen objects, Flow version checks, and spectator-mode loading.
Removing a `.lock` file is a recovery action: first establish that no scheduler or writer is active. Loading a Flow
may call `check_status`; it is not equivalent to inert JSON parsing.

## Paths and external resources

State clearly whether a serialized object embeds data or only refers to files. Robot and Flow JSON commonly stores
paths needed to reconstruct an object; moving or deleting those files breaks reconstruction.

Normalize path-like values consistently, but do not replace absolute paths with relative paths unless a stable base
directory is part of the contract. Never persist machine-specific temporary paths as a supposedly portable fixture.
When a reader is reconstructed from a path, keep ownership and closing behavior consistent with the original class.

Write persistent state atomically when partial writes could corrupt a reusable database. Use AbiPy's established
atomic/file-lock pattern rather than open-coded temporary-file renames around Flow state.

## Generate notebooks deliberately

`write_notebook` should call `get_nbformat_nbv_nb`, append focused Markdown and code cells, and finish with
`_write_nb_nbpath`. Generated code must import from public AbiPy APIs and reconstruct the object from stable inputs.

Do not embed a live object through an opaque temporary pickle unless the existing object explicitly uses that
short-lived pattern. Quote paths safely, avoid machine-specific environment assumptions, and do not execute the
notebook during generation. Notebook creation must not open a browser or call `expose`.

## Tests

Use a pytest temporary directory and test the actual public round trip:

- `MontyEncoder` to `MontyDecoder`, plus `mjson_write`/`mjson_load` when file helpers are part of the API;
- pickle dump/load for supported objects, using the correct specialized loader;
- equality of semantic values, array dtype/shape, units, optional fields, and nested types;
- legacy dictionaries when backward compatibility is intentional;
- type checks and useful failures for incompatible or malformed state;
- slot state, path behavior, and exclusion/reconstruction of transient resources;
- notebook JSON validity, key cells, and explicit output paths without executing a kernel.

Do not merely assert that decoding returns the same class. Compare the state that determines scientific behavior and
verify the reconstructed object remains usable. Avoid committing generated pickle files: their bytes and protocol
are not stable fixtures. A JSON fixture is appropriate only when it represents a maintained compatibility contract.

For Flow persistence, use a temporary Flow and verify graph identity, dependencies, manager/state behavior, version
handling, and spectator mode without submitting jobs. Use `abipy-integration-tests` only when scheduler execution or
real task products are essential.

Report the serialization format and public helpers used, embedded versus referenced state, compatibility behavior,
security or version limitations, and round-trip tests run.
