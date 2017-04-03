"""
Objects used to extract and plot results from output files in text format.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os

from monty.string import is_string
from abipy.flowtk import EventsParser, NetcdfReader, GroundStateScfCycle, D2DEScfCycle
from abipy.abio.timer import AbinitTimerParser
from abipy.core.mixins import TextFile, AbinitNcFile, NotebookWriter


class AbinitTextFile(TextFile):
    """Class for the ABINIT main output file and the log file."""

    @property
    def events(self):
        """
        List of ABINIT events reported in the file.
        """
        # Parse the file the first time the property is accessed or when mtime is changed.
        stat = os.stat(self.filepath)
        if stat.st_mtime != self._last_mtime or not hasattr(self, "_events"):
            self._events = EventsParser().parse(self.filepath)
        return self._events

    def get_timer(self):
        """
        Timer data.
        """
        timer = AbinitTimerParser()
        timer.parse(self.filepath)
        return timer


class AbinitLogFile(AbinitTextFile, NotebookWriter):
    """Class representing the Abinit log file."""

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abilog = abilab.abiopen('%s')" % self.path),
            nbv.new_code_cell("print(abilog.events)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class AbinitOutputFile(AbinitTextFile, NotebookWriter):
    """Class representing the main output file."""

    def __init__(self, filepath):
        super(AbinitOutputFile, self).__init__(filepath)

        # TODO: Parse header to get important dimensions and variables
        #ndtset
        #offset_dataset
        #dims_dataset
        #vars_dataset
        #pseudos
        #self.seek(0)

    def next_gs_scf_cycle(self):
        """
        Return the next :class:`GroundStateScfCycle` in the file. None if not found.
        """
        return GroundStateScfCycle.from_stream(self)

    def next_d2de_scf_cycle(self):
        """
        Return :class:`GroundStateScfCycle` with information on the GS iterations. None if not found.
        """
        return D2DEScfCycle.from_stream(self)

    def compare_gs_scf_cycles(self, others, show=True):
        """
        Produce and returns a list of `matplotlib` figure comparing the GS self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        # Open file here if we receive a string. Files will be closed before returning
        close_files = []
        for i, other in enumerate(others):
            if is_string(other):
                others[i] = self.__class__.from_file(other)
                close_files.append(i)

        fig, figures = None, []
        while True:
            cycle = self.next_gs_scf_cycle()
            if cycle is None: break

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_gs_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(axlist=fig.axes, show=show and last)
                if last:
                    fig.tight_layout()
                    figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)

        if close_files:
            for i in close_files: others[i].close()

        return figures

    def compare_d2de_scf_cycles(self, others, show=True):
        """
        Produce and returns a `matplotlib` figure comparing the DFPT self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        # Open file here if we receive a string. Files will be closed before returning
        close_files = []
        for i, other in enumerate(others):
            if is_string(other):
                others[i] = self.__class__.from_file(other)
                close_files.append(i)

        fig, figures = None, []
        while True:
            cycle = self.next_d2de_scf_cycle()
            if cycle is None: break

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_d2de_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(axlist=fig.axes, show=show and last)
                if last:
                    fig.tight_layout()
                    figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)

        if close_files:
            for i in close_files: others[i].close()

        return figures

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abo = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abo.events)"),
            nbv.new_code_cell("""
abo.seek(0); icycle = -1
while True:
    gs_cycle = abo.next_gs_scf_cycle()
    if gs_cycle is None: break
    icycle += 1
    gs_cycle.plot(title="SCF cycle no %d" % icycle, tight_layout=True)"""),

        nbv.new_code_cell("""
abo.seek(0); icycle = -1
while True:
    d2de_cycle = abo.next_d2de_scf_cycle()
    if d2de_cycle is None: break
    icycle += 1
    d2de_cycle.plot(title="DFPT cycle no %d" % icycle, tight_layout=True)"""),

       nbv.new_code_cell("""
abo.seek(0)
timer = abo.get_timer()
if timer:
    timer.plot_all()
"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class OutNcFile(AbinitNcFile):
    """
    Class representing the _OUT.nc file containing the dataset results
    produced at the end of the run. The netcdf variables can be accessed
    via instance attribute e.g. `outfile.ecut`. Provides integration with ipython.
    """
    def __init__(self, filepath):
        super(OutNcFile, self).__init__(filepath)
        self.reader = NetcdfReader(filepath)
        self._varscache= {k: None for k in self.reader.rootgrp.variables}

    def __dir__(self):
        """Ipython integration."""
        return sorted(list(self._varscache.keys()))

    def __getattribute__(self, name):
        try:
            return super(OutNcFile, self).__getattribute__(name)
        except AttributeError:
            # Look in self._varscache
            varscache = super(OutNcFile, self).__getattribute__("_varscache")
            if name not in varscache:
                raise AttributeError("Cannot find attribute %s" % name)
            reader = super(OutNcFile, self).__getattribute__("reader")
            if varscache[name] is None:
                varscache[name] = reader.read_value(name)
            return varscache[name]

    def close(self):
        self.reader.close()

    def get_allvars(self):
        """
        Read all netcdf variables present in the file.
        Return dictionary varname --> value
        """
        for k, v in self._varscache.items():
            if v is not None: continue
            self._varscache[k] = self.reader.read_value(k)
        return self._varscache
