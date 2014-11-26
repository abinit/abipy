# coding: utf-8
"""Robots."""
from __future__ import print_function, division, unicode_literals

#import numpy as np
import pandas as pd

from collections import OrderedDict, deque 
from monty.string import is_string
from pymatgen.io.abinitio.flows import AbinitFlow


__all__ = [
    "GsrRobot",
    "SigresRobot",
]


class Robot(object):
    """
    The main function of a `Robot` is facilitating the extraction of the output data produced by
    multiple tasks in a `AbinitFlow`. This is the base class from which all Robot subclasses should derive.
    """
    def __init__(self, *args):
        self._ncfiles, self._do_close = OrderedDict(), OrderedDict()
        self._exceptions = deque(maxlen=100)

        for label, ncfile in args:
            self.add_file(label, ncfile)

    def add_file(self, label, ncfile):
        if is_string(ncfile):
            from abipy.abilab import abiopen
            ncfile = abiopen(ncfile)
            self._do_close[ncfile.filepath] = True

        self._ncfiles[label] = ncfile

    def exceptions(self):
        """List of exceptions."""
        return self._exceptions

    def __iter__(self):
        return iter(self._ncfiles.items())

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.close()

    def __str__(self):
        return "\n".join(["%s --> %s" % (label, ncfile.filepath) for label, ncfile in self])

    @property
    def ncfiles(self):
        return self._ncfiles.values()

    def close(self):
        """Close all the files that have been opened by the Robot"""
        for ncfile in self.ncfiles:
            if self._do_close.pop(ncfile.filepath, False): 
                try:
                    ncfile.close()
                except:
                    pass

    @classmethod
    def _from_flow(cls, flow, open_method, **kwargs):
        if is_string(flow): flow = AbinitFlow.pickle_load(flow)

        items = []
        if hasattr(flow, "iflat_tasks"):
            for task in flow.iflat_tasks():
                if hasattr(task, open_method):
                    ncfile = getattr(task, open_method)()
                    if ncfile is not None: items.append((task.pos_str, ncfile))
        else:
            # Work
            for task in flow:
                if hasattr(task, open_method):
                    ncfile = getattr(task, open_method)()
                    if ncfile is not None: items.append((task.pos_str, ncfile))

        return cls(*items)

    @staticmethod
    def _get_geodict(structure):
        """Return a dictionary with info on structure to be used in pandas dataframe."""
        abc, angles = structure.lattice.abc, structure.lattice.angles
        return dict(
            a=abc[0], b=abc[1], c=abc[2], volume=structure.volume,
            angle0=angles[0], angle1=angles[1], angle2=angles[2],
            formula=structure.formula,
        )

    def _exec_funcs(self, funcs, arg):
        # Execute funcs.
        if not isinstance(funcs, (list, tuple)): funcs = [funcs]
        d = {}
        for func in funcs:
            try:
                key, value = func(arg)
                d[key] = value
            except Exception as exc:
                self._exceptions.append(str(exc))
        return d


class GsrRobot(Robot):
    """This robot analyzes the results contained in multiple GSR files."""
    @classmethod
    def from_flow(cls, flow, **kwargs):
        return cls._from_flow(flow, "open_gsr", **kwargs)

    def get_dataframe(self, **kwargs):
        """
        Return a pandas DataFrame with the most important GS results.

        kwargs:
            attrs:
                List of additional attributes of the gsr file to add to the DataFrame
            funcs:
                Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a GSR_File object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # TODO add more columns
        # Add attributes specified by the users
        attrs = [
            "nsppol", "nspinor", "nspden", "ecut", "pawecutdg",
            "tsmear", "nkpts", "energy", "magnetization", "pressure", "max_force",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, gsr in self:
            row_names.append(label)
            d = {aname: getattr(gsr, aname) for aname in attrs}

            # Add info on structure.
            if kwargs.get("with_geo", True):
                d.update(self._get_geodict(gsr.structure))

            # Execute funcs.
            d.update(self._exec_funcs(kwargs.get("funcs", []), gsr))

            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=rows[0].keys())

    def ebands_plotter(self):
        from abipy import abilab
        plotter = abilab.ElectronBandsPlotter()
        for label, gsr in self:
            plotter.add_ebands(label, gsr.ebands)
        return plotter

    def eos_fit(self, eos_name="murnaghan"):
        """
        Fit E(V)
        For the list of available models, see EOS.MODELS

        TODO: which default? all should return a list of fits
        """
        # Read volumes and energies from the GSR files.
        energies, volumes = [], []
        for label, gsr in self:
            energies.append(gsr.energy)
            volumes.append(gsr.structure.volume)

        # Note that eos.fit expects lengths in Angstrom, and energies in eV.
        from abipy.abilab import EOS 
        return EOS(eos_name=eos_name).fit(volumes, energies, vol_unit="ang^3", ene_unit="eV")


class SigresRobot(Robot):
    """This robot analyzes the results contained in multiple SIGRES files."""

    @classmethod
    def from_flow(cls, flow, **kwargs):
        return cls._from_flow(flow, "open_sigres", **kwargs)

    def merge_dataframes_sk(self, spin, kpoint, **kwargs):
        for i, (label, sigr) in enumerate(self):
            frame = sigr.get_dataframe_sk(spin, kpoint, index=label)
            if i == 0:
                table = frame
            else:
                table = table.append(frame)
        
        return table

    def get_qpgaps_dataframe(self, spin=None, kpoint=None, **kwargs):
        # TODO: Ideally one should select the k-point for which we have the fundamental gap for the given spin
        if spin is None: spin = 0
        if kpoint is None: kpoint = 0

        attrs = [
            "nsppol", "nspinor", "nspden", #"ecut", "pawecutdg",
            #"tsmear", "nkpts",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for i, (label, sigr) in enumerate(self):
            row_names.append(label)
            d = {aname: getattr(sigr, aname) for aname in attrs}
            d.update({"qpgap": sigr.get_qpgap(spin, kpoint)})
            rows.append(d)

            # Add convergence parameters
            d.update(sigr.params)
                                                        
            # Add info on structure.
            if kwargs.get("with_geo", False):
                d.update(self._get_geodict(sigr.structure))

            # Execute funcs.
            d.update(self._exec_funcs(kwargs.get("funcs", []), sigr))

        return pd.DataFrame(rows, index=row_names, columns=rows[0].keys())


# TODO
#class PhononRobot(Robot):
#    #@classmethod
#    #def from_flow(cls, flow, **kwargs):
#    #    return cls._from_flow(flow, "open_sigres", **kwargs)
#
#    def get_ph_dataframe(self, qpoint=None, **kwargs):
#        if kpoint is None: kpoint = 0
#
#        attrs = [
#            "nsppol", "nspinor", "nspden", #"ecut", "pawecutdg",
#            #"tsmear", "nkpts",
#        ] + kwargs.pop("attrs", [])
#
#        rows, row_names = [], []
#        for i, (label, sigr) in enumerate(self):
#            row_names.append(label)
#            d = {aname: getattr(sigr, aname) for aname in attrs}
#            #d.update({"qpgap": sigr.get_qpgap(spin, kpoint)})
#            rows.append(d)
#
#            # Add convergence parameters
#            d.update(sigr.params)
#                                                        
#            # Add info on structure.
#            if kwargs.get("with_geo", False):
#                d.update(self._get_geodict(sigr.structure))
#
#            # Execute funcs.
#            d.update(self._exec_funcs(kwargs.get("funcs", []), sigr))
#
#        return pd.DataFrame(rows, index=row_names, columns=rows[0].keys())
