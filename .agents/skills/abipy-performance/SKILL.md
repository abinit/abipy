---
name: abipy-performance
description: Measure, diagnose, optimize, and regression-test AbiPy and ABINIT performance, including Python CPU, memory and import profiling, large-file access, FFT results, and BenchmarkFlow parallel scaling. Use when runtime or resource behavior is an explicit goal; not for ordinary correctness testing.
---

# Work on AbiPy performance

Start with a reproducible measurement and optimize the dominant cost while preserving scientific results. Separate
Python performance from ABINIT executable scaling: they need different profilers, environments, and evidence. Use
`python-virtual-env` before Python measurements, `netcdf` for large-file layout issues, `abipy-flowtk` for execution
mechanics, and `abipy-plotting` for presentation of benchmark results.

## Define the measurement

Record the operation, representative input, dataset size and shape, environment, dependency versions, machine and
CPU count, thread settings, cold/warm-cache state, and metric. Decide whether the target is elapsed time, CPU time,
peak/resident memory, allocations, import latency, I/O volume, or parallel efficiency.

Use the smallest workload that retains the bottleneck, but confirm the result on a realistic size before claiming an
improvement. Measure a baseline and candidate in the same process model and environment. Repeat short operations,
include warm-up where imports, JITs, filesystem cache, or lazy properties matter, and report a distribution or robust
summary rather than one best run.

Do not mix correctness setup, queue wait time, plotting/display time, and the code under investigation unless they
are deliberately part of the metric. Avoid timing with debug logging or an interactive backend enabled.

## Profile Python and imports

Use a coarse timer first, then a profiler that matches the question. `cProfile` and `pstats` identify cumulative CPU
cost; `python -X importtime -c "import abipy"` plus the repository's `invoke tuna` workflow diagnoses imports;
allocation or RSS tools are needed for memory. Do not infer a memory leak from peak RSS alone.

Profile the actual public call with representative data. Inspect cumulative time, call count, and callers before
editing. A hot function may be expensive because its caller invokes it redundantly, not because its inner loop is
slow. Compare profiles after the change using the same command and input.

Keep optional dependencies lazy and avoid moving heavyweight imports into AbiPy's common import path. When reducing
import time, verify both `import abipy` and the first use of the moved feature so initialization cost has not merely
become a surprising runtime spike.

## Optimize scientific Python safely

Preserve units, indexing, shapes, dtypes, masks, tolerances, ordering, and lazy-loading behavior. Add or retain a
correctness test before changing vectorization, caching, chunking, parallelism, or numerical libraries.

For NetCDF-backed objects, read only needed variables and slices, avoid repeated decoding in loops, and distinguish
disk I/O from transformations. Cache an immutable expensive property only when the object's source and parameters
cannot change; bound caches whose keys can grow. Avoid eager materialization of a large array to save a small number
of calls.

Vectorization is useful when it reduces Python overhead without creating much larger temporaries. Estimate temporary
array sizes and preserve precision. Chunk work when full broadcasting increases peak memory beyond realistic files.
Do not add multiprocessing or threads until the profile shows parallel work large enough to repay startup, copying,
serialization, and oversubscription costs.

## Benchmark ABINIT scaling

Use scripts under `abipy/benchmarks` and `BenchmarkFlow` for ABINIT MPI/OpenMP studies. A benchmark script should
build a scientifically fixed workload while varying only the intended parallel configuration. Reuse `bench_main`,
`build_bench_main_parser`, and the patched options:

- `--mpi-list` and `--omp-list` define tested layouts;
- `--min-ncpus`, `--max-ncpus`, and `--min-eff` bound the matrix;
- `options.accept_mpi_omp` or `accept_conf` filters configurations;
- `options.get_workdir(__file__)` supplies the conventional benchmark work directory.

Use `BenchmarkFlow.exclude_from_benchmark` for prerequisite tasks whose timing is not part of the comparison. Do not
silently change cutoffs, meshes, algorithms, convergence criteria, task manager, or executable between scaling
points. Record MPI processes, OpenMP threads, total cores, nodes, affinity/binding, manager/queue adapter, ABINIT
version/build, compiler, math/FFT libraries, and relevant hardware.

`BenchmarkFlow.get_parser` includes successful non-excluded tasks and delegates to the ABINIT timer parser. Check for
missing timing sections and failed tasks before analyzing speedup. Define the reference configuration explicitly:

```text
speedup(N) = T(reference) / T(N)
efficiency(N) = speedup(N) * reference_cores / N
```

Do not compare different scientific outputs as if they were scaling points. Verify representative energies, forces,
iteration counts, or other invariants when a parallel algorithm could change the computation.

Building a benchmark flow does not execute it. `--scheduler` can submit real jobs, and running a local benchmark can
consume substantial CPU time. Obtain the required authorization before execution or submission, bound the resource
matrix, and use a new or explicitly disposable work directory. Never remove an existing benchmark directory merely
because `--remove` is available.

## FFT and timing results

Use `FFTBenchmark.from_file` to analyze existing ABINIT `PROF_*` output and `FFTProf` only when deliberately running
the external FFT profiler. Compare algorithms and thread counts against a clearly identified reference. Preserve
FFT grid and cutoff alignment; do not divide arrays from different problem sizes to compute speedup.

Use `AbinitTimerParser` for ABINIT timing sections. Treat elapsed, CPU, and section percentages as different metrics,
and confirm all expected files parsed. Queue wait time and scheduler overhead should be reported separately from
ABINIT wall time.

## Regression tests and acceptance

Ordinary unit tests should validate parsing, data reduction, configuration generation, and scientific equivalence,
not assert tight wall-clock limits. Shared CI hosts are noisy. Use a dedicated benchmark or a generously bounded
smoke threshold only when a catastrophic complexity regression cannot be detected structurally.

For a claimed optimization, report baseline and candidate samples, median or another stated statistic, variability,
speedup, peak memory when relevant, and unchanged correctness checks. Prefer an input-size series when algorithmic
complexity is the concern. A result within measurement noise is inconclusive, not an improvement.

Tests for `abipy/benchmarks` should construct flows in temporary directories, verify task counts and resource
configurations, and build/pickle without scheduling. Tests for parsers should use maintained small output files. Do
not run the full benchmark matrix or open graphical windows in pytest.

## Completion information

Report the bottleneck demonstrated by the profile, benchmark command and input, environment and hardware, repetitions
and variability, before/after metrics, numerical equivalence checks, memory tradeoffs, and whether any external jobs
were submitted. Keep raw profiles and large benchmark outputs out of the repository unless they are intentional,
small reference artifacts.
