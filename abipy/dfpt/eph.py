# coding: utf-8
"""This module contains objects for the postprocessing of EPH calculations."""
from __future__ import print_function, division, unicode_literals

import numpy as np

from abipy.core.func1d import Function1D
from abipy.iotools import ETSF_Reader


class EliashbergFunction(object):
    """This object stores the Eliashberg function a2F(w)."""

    def __init__(self, mesh, values):
        values = np.atleast_2d(values)
        self.nsppol = len(values)

        self.a2f_spin, a2f_tot = [], np.zeros(len(mesh))

        for vs in values:
            self.a2f_spin.append(Function1D(mesh, vs))
            a2f_tot += vs

        # Spin dependent and total a2F(w)
        self.a2f_spin = tuple(self.a2f_spin)
        self.a2f = Function1D(mesh, a2f_tot)

    @classmethod
    def from_file(cls, path):
        """Initialize the object from the data stored in a netcdf file."""
        with ETSF_Reader(path) as r:
            mesh = r.get_a2fmesh()
            values = r.get_a2fvalues()
            return cls(mesh, values)

    def get_lambdaw(self, spin=None):
        """Returns lambda(w)."""
        raise NotImplementedError()

    def get_momentum(self, order, spin=None):
        """Computes the momenta of a2F(w)"""
        raise NotImplementedError()

    def get_mcmillan_Tc(self, mustar):
        """
        Computes the critical temperature with the McMillan equation and the input mustar
        """""
        raise NotImplementedError()

    def plot(self, spin=None):
        """
        Plot a2F(w) with matplotlib.
        """""
        raise NotImplementedError()


class EPH_Reader(ETSF_Reader):

    def __init__(self, file):
        super(EPH_Reader, self).__init__(file)

        self.eph_lambda = self.read_value("eph_coupling")
        raise NotImplementedError("")

    def read_qpoints(self):
        """List of q-points where the EPH quanties are computed."""
        return qpoints_factor(self.path)

    def read_lambda(self, spin, qpoint=None): 
        """
        Reads the EPH couping matrix elements 

        spin: Spin index.
        qpoint: Qpoint or int. If None all the matrix elements for this spin are read.
        """
