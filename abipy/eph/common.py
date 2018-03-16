# coding: utf-8
"""
Objects common to the other eph modules.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.functools import lazy_property
from abipy.electrons.ebands import ElectronsReader


class BaseEphReader(ElectronsReader):
    """
    Provides methods common to the netcdf files produced by the EPH code.
    See Abinit docs for the meaning of the variables.
    """

    @lazy_property
    def ddb_ngqpt(self):
        """Q-Mesh for DDB file."""
        return self.read_value("ddb_ngqpt")

    @lazy_property
    def ngqpt(self):
        """Effective Q-mesh used in to compute integrals (ph_linewidts, e-ph self-energy)."""
        return self.read_value("ngqpt")

    @lazy_property
    def ph_ngqpt(self):
        """Q-mesh for Phonon DOS, interpolated A2F ..."""
        return self.read_value("ph_ngqpt")

    @lazy_property
    def eph_ngqpt_fine(self):
        """Q-mesh for interpolated DFPT potentials"""
        return self.read_value("eph_ngqpt_fine")

    @lazy_property
    def common_eph_params(self):
        """
        Read basic parameters (scalars) from the netcdf files produced by the EPH code and cache them
        """
        od = OrderedDict([
            ("ddb_nqbz", np.prod(self.ddb_ngqpt)),
            ("eph_nqbz_fine", np.prod(self.eph_ngqpt_fine)),
            ("ph_nqbz", np.prod(self.ph_ngqpt)),
        ])

        for vname in ["eph_intmeth", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"]:
            value = self.read_value(vname)
            if vname in ("eph_intmeth",):
                value = int(value)
            else:
                value = float(value)
            od[vname] = value

        return od
