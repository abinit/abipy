from __future__ import print_function, division

import collections
import subprocess

from abipy.tools import which

__all__ = [
    "animator",
]


class Animator(object):

    DEFAULT_OPTIONS = {
        "-delay": "100", # display the next image after pausing <1/100ths of a second>
    }

    def __init__(self):
        self._figures = collections.OrderedDict()

        self.animate_bin = which("animate")

        if self.animate_bin is None:
            raise RuntimeError("Cannot find animate executable in $PATH.\n Please install the ImageMagick suite of tools.")

    def add_figure(self, label, figure):
        """
        Add a figure.

        Args:
            label:
            figure:
        """
        if label in self._figures:
            raise ValueError("label %s is already in %s" % (label, self._figures.keys()))

        self._figures[label] = figure

    def add_figures(self, labels, figure_list):
        """
        Add a list of figures.

        Args:
            labels:
                List of labels.
            figure_list:
                List of Figures
        """
        assert len(labels) == len(figure_list)

        for label, figure in zip(labels, figure_list):
            self.add_figure(label, figure)

    def animate(self, **kwargs):
        figs = list(self._figures.values())

        options = []
        if not kwargs:
            for k,v in self.DEFAULT_OPTIONS.items():
                options.extend([k, v])
        else:
            for k,v in kwargs.items():
                options.extend([k, v])

        command = [self.animate_bin] + options + figs
        print(command)

        retcode = subprocess.call(command)
        return retcode

