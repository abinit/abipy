# coding: utf-8
"""GSR file."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict, Iterable, defaultdict
from monty.string import is_string, list_strings
from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands
from prettytable import PrettyTable
from .ebands import ElectronsReader

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "HirshfeldCharges",
]


class HirshfeldCharges(Has_Structure, object):

    def __init__(self, charges, structure):
        self.charges = charges
        self._structure = structure

    def structure(self):
        return self._structure

    @classmethod
    def from_cut3d_outfile(cls, filepath, structure):
        charges = []
        with open(filepath, 'rt') as f:
            while True:
                l = f.readline()

                if l is None:
                    raise RuntimeError('The file does not contain Hirshfeld harges')
                elif "Hirshfeld Charge" in l:
                    break

            #skip one line
            f.readline()

            for i in range(len(structure)):
                l = f.readline()
                charges.append(float(l.split()[2]))

        return cls(charges, structure)


