# coding: utf-8
"""Tests for fftprof module."""
import os
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.tools.fftprof import FFTBenchmark


class FftProfTest(AbipyTest):

    def test_fft_benchmark(self):
        """Testing FFt benchmark."""

        # Plot the benchmark results saved in the files
        path = os.path.join(abidata.dirpath, "PROF_fourwf_cplex0_option3_istwfk1")
        bench = FFTBenchmark.from_file(path)

        assert bench.title.strip() == " Benchmark: routine = fourwf, cplex = 0, option= 3, istwfk= 1".strip()
        assert len(bench.tests) == 4
        assert len(bench.tests_with_fftalg(112)) == 2
        assert len(bench.tests_with_fftalg(512)) == 2

        test0 = bench.tests_with_fftalg(112)[0]
        test1 = bench.tests_with_fftalg(112)[1]
        assert len(test0.ecut) == len(test0.ngfft)
        assert str(test0)
        test0.speedup_wrt(test1)

        if self.has_matplotlib():
            #test0.plot_ax()
            assert bench.plot(show=False)
