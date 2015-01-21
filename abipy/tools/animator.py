# coding: utf-8
from __future__ import print_function, division, unicode_literals

import os
import collections
import subprocess
import tempfile

from monty.os.path import which

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "FilesAnimator",
]


class FilesAnimator(object):
    """
    This object wraps the animate tool of ImageMagick.
    If received a list of files and/or a list of `matplotlib` figures 
    and creates an animation. 

    .. note:
        The module `matplotlib.animation` provides a much more powerful API for the 
        creation of animation. `FilesAnimator` should only be used when we already 
        have a list of files and figures and we want to create a simple animation.
    """
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
        """
        Calls `animate` in a subprocess.

        Returns:
            return code of the subprocess.
        """
        figs = list(self._figures.values())

        options = []
        if not kwargs:
            for k,v in self.DEFAULT_OPTIONS.items():
                options.extend([k, v])
        else:
            for k,v in kwargs.items():
                options.extend([k, v])

        # We have to create files before calling animate.
        # If we have a matplotlib figure we call savefig 
        tmp_dirname = tempfile.mkdtemp()
        files = figs

        files = []
        for i, fig in enumerate(figs):

            if hasattr(fig, "savefig"):
                # matplotlib figure --> save it.
                fname = os.path.join(tmp_dirname, "scratch_" + str(i) + ".png")
                fig.savefig(fname)
            else:
                # Assume filepath.
                fname = fig

            assert os.path.isfile(fname)
            files.append(fname)

        command = [self.animate_bin] + options + files
        msg = "will execute command: %s" % command
        logger.info(msg)
        retcode = subprocess.call(command)

        # Remove the temporary directory.
        #os.rmdir(tmp_dirname)

        return retcode

