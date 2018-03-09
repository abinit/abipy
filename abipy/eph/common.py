# coding: utf-8
"""
Objects common to the other eph modules.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import warnings

from collections import OrderedDict
from abipy.electrons.ebands import ElectronsReader


class BaseEphReader(ElectronsReader):
    """
    Provides methods common to the netcdf files produced by the EPH code.
    """

    def read_base_eph_params(self):
        """
        Read basic parameters from the netcdf files produced by the EPH code.
        See Abinit docs for the meaning of the variables.
        """
        try:
            ddb_ngqpt = self.read_value("ddb_ngqpt")
            eph_ngqpt_fine = self.read_value("eph_ngqpt_fine")
            ddb_nqbz = np.prod(ddb_ngqpt)
            eph_nqbz = np.prod(eph_ngqpt_fine)
        except Exception:
            warnings.warn("ddb_ngqpt, eph_ngqpt_fine not in file.")
            ddb_nqbz = None
            eph_nqbz = None

        od = OrderedDict([
            ("ddb_nqbz", ddb_nqbz),
            ("eph_nqbz", eph_nqbz),
        ])
        for vname in ["eph_intmeth", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"]:
            value = self.read_value(vname, default=None)
            if value is None:
                warnings.warn("%s not in file." % vname)
                continue
            if vname in ("eph_intmeth", ):
                value = int(value)
            else:
                value = float(value)
            od[vname] = value

        return od
