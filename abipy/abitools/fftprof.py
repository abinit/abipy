from __future__ import print_function, division

import numpy as np

__all__ = [
    "FFT_Benchmark",
]

_color_fftalg = {
    112: "r",
    #311: "y",
    312: "b",
    412: "y",
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
            ecut:
                List of cutoff energies.
            ngfft:
                List of FFT divisions.
            wall_time:
                List of wall_time for the different ecut.
            info:
                Dictionary with extra information.
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
                        marker=marker, markersize=10,
                        )
        return line

######################################################################


class FFT_Benchmark(object):
    """
    Container class storing the results of the FFT benchmark.

    Use the classmethod from_filename to generate an instance.
    """
    @classmethod
    def from_filename(cls, fileobj):
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

    def plot(self, exclude_algs=None, exclude_threads=None, **kwargs):
        """
        Plot the wall-time and the speed-up.

        Args:

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        show            True to show the figure (Default)

        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ===============================================================
        """
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

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
        #labels = list()
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

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(os.path.abspath(savefig))

        return fig

def parse_prof_file(fileobj):
    """
    Parse the PROF file generated by fftprof.F90.

    Args:
        fileobj:
            String or file-like object.

    Returns:
        Instance of FFT_Benchmark.
    """
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
        #
    # Instanciate FFT_Benchmark.
    fft_tests = []
    for (idx, wall_time) in enumerate(data):
        info = info_of_test[idx]
        Test = FFT_Test(ecut, ngfft, wall_time, info)
        fft_tests.append(Test)

    return FFT_Benchmark(title, fft_tests)
