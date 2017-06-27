# coding: utf-8
"""HirshfeldCharges."""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.core.mixins import Has_Structure

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "HirshfeldCharges",
]


class HirshfeldCharges(Has_Structure):

    def __init__(self, charges, structure):
        self.charges = charges
        self._structure = structure

    @property
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
