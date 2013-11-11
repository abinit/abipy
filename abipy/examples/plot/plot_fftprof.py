#!/usr/bin/env python
#
# This example shows how to plot the benchmark results for the FFT libraries produced by fftprof.F90
# To run the graphical interface to fftprof see abipy.gui.demos.demo_fftprof.py
from abipy.abitools.fftprof import FFT_Benchmark
import abipy.data as abidata

# Path to the PROF file produced by fftprof.
prof_filename = abidata.ref_file("PROF_fourwf_cplex0_option3_istwfk1")

# Create the object
bench = FFT_Benchmark.from_file(prof_filename)

# Plot data.
bench.plot()
