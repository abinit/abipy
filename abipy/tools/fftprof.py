"""
Python interface to fftprof. Provides objects to benchmark
the FFT libraries used by ABINIT and plot the results with matplotlib.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import tempfile
import numpy as np

from subprocess import Popen, PIPE
from monty.os.path import which
from monty.fnmatch import WildCard
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt

__all__ = [
    "FFTBenchmark",
]

_color_fftalg = {
    112: "r",
    #311: "y",
    312: "b",
    412: "y",
    401: "c",
    #511: "c",
    512: "g",
}

_linestyle_nt = {
    1: "-",
    2: "--",
    3: "-.",
    4: ":",
    5: "-",
    6: "--",
    7: "-.",
    8: ":",
}

_markers_nt = {
    1: "1",
    2: ">",
    3: "2",
    4: "<",
    5: "3",
    6: "*",
    7: "4",
    8: "+",
}


class FFT_Test(object):
    """
    This object stores the wall-time of the FFT
    as a function of the size of the problem.
    """
    def __init__(self, ecut, ngfft, wall_time, info):
        """
        Args:
            ecut: List of cutoff energies in Hartree.
            ngfft: List with FFT divisions.
            wall_time: List of wall_time for the different ecut.
            info: Dictionary with extra information.
        """
        self.ecut = np.asarray(ecut)
        self.necut = len(ecut)
        self.ngfft = np.asarray(ngfft, np.int)
        self.wall_time = np.asarray(wall_time)

        self._info = info
        self.__dict__.update(info)
        self.available = int(self.available)

        self.fftalg = int(self.fftalg)
        self.nthreads = int(self.nthreads)

    def __str__(self):
        """String representation."""
        return "fftalg: %s, nt: %d" % (self.fftalg, self.nthreads)

    def speedup_wrt(self, other):
        """Compare self with other, return a new instance of FFT_Test containing the speedup."""
        info = self._info.copy()
        info["speedup"] = other.wall_time / self.wall_time
        return FFT_Test(self.ecut, self.ngfft, self.wall_time, info)

    def plot_ax(self, ax, ydata=None, fact=1.0):
        """Plot ydata on the axis ax. If data is None, the wall_time is plotted."""
        color = _color_fftalg[self.fftalg]
        linestyle = _linestyle_nt[self.nthreads]
        marker = _markers_nt[self.nthreads]

        xx = self.ecut

        if ydata is None:
            yy = self.wall_time
        else:
            yy = self.__dict__[ydata]

        yy = yy * fact

        line, = ax.plot(xx, yy,
                        color=color, linestyle=linestyle, label=str(self), linewidth=3.0,
                        marker=marker, markersize=10)
        return line


class FFTBenchmark(object):
    """
    Container class storing the results of the FFT benchmark.

    Use the class method ``from_file`` to generate a new instance.
    """
    @classmethod
    def from_file(cls, fileobj):
        return parse_prof_file(fileobj)

    def __init__(self, title, FFT_tests):
        self.title = title

        self._fftalgs = []
        for test in FFT_tests:
            alg = test.fftalg
            if alg not in self._fftalgs:
                self._fftalgs.append(alg)

        self._fftalgs.sort()
        self.tests = [t for t in FFT_tests]

    def iter_fftalgs(self):
        """Iterator over the FFT algorithms."""
        return self._fftalgs.__iter__()

    def tests_with_fftalg(self, fftalg):
        """Return the list of FFT_tests with a given fftalg."""
        lst = []
        for t in self.tests:
            if t.fftalg == fftalg: lst.append(t)
        return lst

    @add_fig_kwargs
    def plot(self, exclude_algs=None, exclude_threads=None, **kwargs):
        """
        Plot the wall-time and the speed-up.
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)

        exc_algs = []
        if exclude_algs is not None:
            exc_algs = [int(alg) for alg in exclude_algs]

        exc_nths = []
        if exclude_threads is not None:
            exc_nths = [int(nth) for nth in exclude_threads]

        for test in self.tests:
            if test.fftalg in exc_algs: continue
            if test.available and test.nthreads not in exc_nths:
                #print("test", test.nthreads)
                test.plot_ax(ax1)

        ax1.legend(loc="upper left")
        ax1.grid(True)
        ax1.set_title(self.title)
        ax1.set_xlabel('Ecut [Hartree]')
        ax1.set_ylabel('WALL time (s)')

        ax2 = fig.add_subplot(2, 1, 2)
        ax2.grid(True)

        ax2.set_xlabel('FFT divisions')
        ax2.set_ylabel('Efficiency')

        # Use FFT divs as labels.
        xticks = ax1.get_xticks()
        #print(xticks)

        t0 = self.tests[0]
        labels = [str(ndiv) for ndiv in t0.ngfft]
        #labels = []
        #for t ndiv in zip(t0.ecut, t0.ngfft):
        #for #labels = [ str(ndiv) for ndiv in t0.ngfft ]
        #print labels

        # Rotate labels.
        ax2.set_xticklabels(labels, fontdict=None, minor=False, rotation=35)

        for fftalg in self.iter_fftalgs():
            if fftalg in exc_algs: continue
            tests = self.tests_with_fftalg(fftalg)
            for t in tests:
                if t.nthreads == 1:
                    ref_test = t
                    break
            else:
                raise ValueError("Ref file not found")

            for t in tests:
                tspeed = t.speedup_wrt(ref_test)
                fact = 1.0 / t.nthreads
                if t.nthreads != 1 and t.available and t.nthreads not in exc_nths:
                    tspeed.plot_ax(ax2, ydata="speedup", fact=fact)

        # Use FFT divs as labels.
        xticks = ax1.get_xticks()
        #print(xticks)

        t0 = self.tests[0]
        labels = []
        for xtick in xticks:
            xecut = float(xtick)
            for ecut, ndiv in zip(t0.ecut, t0.ngfft):
                if abs(ecut - xecut) < 0.1:
                    #print(ecut, xecut, ndiv)
                    labels.append(str(ndiv))
                    break
            else:
                msg = "xecut" + str(xecut) + " not found"
                labels.append("")
                #raise RuntimeError(msg)
                #print("labels:", labels)

        # Set and rotate labels.
        ax2.set_xticklabels(labels, fontdict=None, minor=False, rotation=35)

        ideal = [1.0 for i in range(t0.necut)]
        ax2.plot(t0.ecut, ideal, "b-", linewidth=3.0)

        return fig


def parse_prof_file(fileobj):
    """
    Parse the PROF file generated by fftprof.F90.

    Args: fileobj: String or file-like object.

    Returns: Instance of FFTBenchmark.
    """
    # The file contains
    # 1) An initial header with info on the routine.
    # 2) The list of fftalgs that have been analyzed.
    # Results:
    #   ecut ngfft(1:3) wall_times
    #
    # Example

    # Benchmark: routine = fourwf, cplex = 1, option= 2, istwfk= 1
    #  fftalg = 112, fftcache = 16, ndat = 1, nthreads = 1, available = 1
    #  fftalg = 401, fftcache = 16, ndat = 1, nthreads = 1, available = 1
    #  fftalg = 412, fftcache = 16, ndat = 1, nthreads = 1, available = 1
    #  fftalg = 312, fftcache = 16, ndat = 1, nthreads = 1, available = 1
    #   20.0  91  91  90 0.0307 0.0330 0.0275 0.0283
    #   30.0 101 101 100 0.0395 0.0494 0.0404 0.0333

    if not hasattr(fileobj, "readlines"):
        with open(fileobj, "r") as fh:
            lines = fh.readlines()
    else:
        lines = fileobj.readlines()

    # Parse the header
    title = lines.pop(0)[1:]
    line = lines.pop(0)
    info_of_test = []
    while line[0] == "#":
        line = line.replace("#", "")
        line = line.replace("\n", "")
        tokens = line.split(",")
        tokens = [tk.replace(" ", "") for tk in tokens]
        info = {}
        for tk in tokens:
            kw = tk.split("=")
            info[kw[0]] = kw[1]
        info_of_test.append(info)
        line = lines.pop(0)

    # Parse the columns: ecut ngfft(4:6) wall_time[1] wall_time[2] ...
    ntests = len(info_of_test)
    data = [[] for i in range(ntests)]

    ecut, ngfft = [], []
    while line:
        vals = [float(v) for v in line.split()]
        ecut.append(vals[0])
        ngfft.append(vals[1:4])
        for idx, v in enumerate(vals[4:]):
            data[idx].append(v)
        try:
            line = lines.pop(0)
        except IndexError:
            line = None

    # Instantiate FFTBenchmark.
    fft_tests = []
    for idx, wall_time in enumerate(data):
        info = info_of_test[idx]
        Test = FFT_Test(ecut, ngfft, wall_time, info)
        fft_tests.append(Test)

    return FFTBenchmark(title, fft_tests)


class FFTProfError(Exception):
    """Exceptions raised by FFTprof."""


class FFTProf(object):
    """Wrapper around fftprof Fortran executable."""
    Error = FFTProfError

    def __init__(self, fft_input, executable="fftprof"):
        self.verbose = 1
        self.fft_input = fft_input

        self.executable = which(executable)
        if self.executable is None:
            raise self.Error("Cannot find executable %s in $PATH" % executable)

    @classmethod
    def from_file(cls, filename, executable="fftprof"):
        with open(filename, "r") as fh:
            fft_input = fh.read()

        return cls(fft_input, executable=executable)

    def run(self):
        """Execute fftprof in a subprocess."""
        self.workdir = tempfile.mkdtemp()
        print(self.workdir)

        self.stdin_fname = os.path.join(self.workdir, "fftprof.in")
        self.stdout_fname = os.path.join(self.workdir, "fftprof.out")
        self.stderr_fname = os.path.join(self.workdir, "fftprof.err")

        with open(self.stdin_fname, "w") as fh:
            fh.write(self.fft_input)

        args = [self.executable, "<", self.stdin_fname, ">", self.stdout_fname, "2>", self.stderr_fname]

        cmd_str = " ".join(args)

        p = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.workdir)

        (self.stdout_data, self.stderr_data) = p.communicate()

        self.returncode = p.returncode

        if self.returncode != 0:
            with open(self.stdout_fname, "r") as out, open(self.stderr_fname, "r") as err:
                self.stdout_data = out.read()
                self.stderr_data = err.read()

            if self.verbose:
                print("*** stdout: ***\n", self.stdout_data)
                print("*** stderr: ***\n", self.stderr_data)

            raise self.Error("%s returned %s\n cmd_str: %s" % (self, self.returncode, cmd_str))

        return self.returncode

    def plot(self):
        filepaths = WildCard("PROF_*").filter(os.listdir(self.workdir))
        filepaths = filter(os.path.isfile, [os.path.join(self.workdir, f) for f in filepaths])

        for prof_file in filepaths:
            if self.verbose:
                print("About to plot prof_file: ", prof_file)

            bench = FFTBenchmark.from_file(prof_file)
            bench.plot()
