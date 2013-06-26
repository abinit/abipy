# This example shows how to plot the benchmark 
# results for the FFT libraries produced by fftprof

from abipy.tests import get_reference_file
from abipy.abitools.fftprof import FFT_Benchmark

# Path to the PROF file produced by fftprof.
prof_filename = get_reference_file("PROF_fourwf_cplex0_option3_istwfk1")

# Create the object
bench = FFT_Benchmark.from_filename(prof_filename)

# Plot data.
bench.plot()
