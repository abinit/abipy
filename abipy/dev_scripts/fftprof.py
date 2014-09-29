#!/usr/bin/env python
"""
Python interface to fftprof. Provides objects to benchmark
the FFT libraries used by ABINIT and plot the results with matplotlib.
"""
from __future__ import print_function, division, unicode_literals

import sys

from abipy.htc.fftbench import FFT_Benchmark


def main():
    try:
        prof_files = sys.argv[1:]
    except IndexError:
        raise RuntimeError("Must specify one or more files")

    # Plot the benchmark results saved in the files
    for prof_file in prof_files:
        FFT_Benchmark.from_file(prof_file).plot()

    return 0


if __name__ == "__main__":
    sys.exit(main())
