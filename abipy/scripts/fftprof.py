#!/usr/bin/env python
from __future__ import print_function, division

import sys

from abipy.abitools.fftprof import FFT_Benchmark


def main():
    fileobj = sys.stdin
    bench = FFT_Benchmark.from_file(fileobj)

    bench.plot(exclude_algs=[412], exclude_threads=[3, 5, 7])
    return 0

if __name__ == "__main__":
    sys.exit(main())






