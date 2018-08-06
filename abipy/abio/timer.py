"""
This module provides objects for extracting timing data from the ABINIT output files
It also provides tools to analyze and to visualize the parallel efficiency.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.core.mixins import NotebookWriter
from abipy.flowtk import AbinitTimerParser as _Parser


class AbinitTimerParser(_Parser, NotebookWriter):

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_stacked_hist(show=False)
        yield self.plot_efficiency(show=False)
        yield self.plot_efficiency(show=False)
        yield self.plot_pie(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)
        filenames = self.filenames

        nb.cells.extend([
            nbv.new_code_cell("parser = abilab.AbinitTimerParser()\nparser.parse(%s)" % str(filenames)),
            nbv.new_code_cell("display(parser.summarize())"),
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("parser.plot_stacked_hist();"),
            nbv.new_code_cell("parser.plot_efficiency(what='good', nmax=5);"),
            nbv.new_code_cell("parser.plot_efficiency(what='bad', nmax=5);"),
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("parser.plot_pie();"),

            nbv.new_code_cell("""\
for timer in parser.timers():
    print(timer)
    display(timer.get_dataframe())"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)
