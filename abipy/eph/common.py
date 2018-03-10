# coding: utf-8
"""
Objects common to the other eph modules.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import warnings

from collections import OrderedDict
from monty.functools import lazy_property
from abipy.electrons.ebands import ElectronsReader


class BaseEphReader(ElectronsReader):
    """
    Provides methods common to the netcdf files produced by the EPH code.
    """

    @lazy_property
    def common_eph_params(self):
        """
        Read basic parameters from the netcdf files produced by the EPH code and cache them
        See Abinit docs for the meaning of the variables.
        """
        try:
            ddb_ngqpt = self.read_value("ddb_ngqpt")
            ddb_nqbz = np.prod(ddb_ngqpt)
            eph_ngqpt_fine = self.read_value("eph_ngqpt_fine")
            eph_nqbz_fine = np.prod(eph_ngqpt_fine)
            #ph_ngqpt = self.read_value("ph_ngqpt")
            #ph_nqbz = np.prod(ph_ngqpt)
        except Exception:
            # TODO: Remove
            warnings.warn("ddb_ngqpt, eph_ngqpt_fine not in file.")
            ddb_ngqpt = None
            ddb_nqbz = None
            eph_ngqpt_fine = None
            eph_nqbz_fine = None
            #ph_ngqpt = None
            #ph_nqbz= None

        od = OrderedDict([
            ("ddb_ngqpt", ddb_ngqpt),
            ("ddb_nqbz", ddb_nqbz),
            ("eph_ngqpt_fine", eph_ngqpt_fine),
            ("eph_nqbz_fine", eph_nqbz_fine),
            #("ph_ngqpt", ph_ngqpt),
            #("ph_nqbz", ph_nqbz),
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
