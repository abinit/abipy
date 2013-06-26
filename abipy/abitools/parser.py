"""Abinit parser."""
from __future__ import print_function, division

import collections
import json
import numpy as np

from abipy.core.structure import Structure

__all__ = [
    "Datasets",
]

_parser_ok = True
try:
    from .ab6_invars import get_ids, Dtsets
    #: list of abinit variable names.
    abinit_vars = get_ids()
except ImportError:
    _parser_ok = False
    #import warnings
    #warnings.warn("Import of ab6_invars failed. abinit Parser is not available")

_abipy2abinit = {
    "rprimd"    : "rprimd_orig",
    "xred"      : "xred_orig",
    "znucl_type": "znucl",
}


class NumPyArangeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist() # or map(int, obj)
        return json.JSONEncoder.default(self, obj)

##########################################################################################


class Dataset(dict):

    def __init__(self, idx, *args, **kwargs):
        super(Dataset, self).__init__(*args, **kwargs)
        self.idx = idx

        #for k,v in _abipy2abinit.iteritems(): self[k] = self[v]
        #if idx>0: # Init CrystalStruct
        #    d = dict()
        #    for key in CrystalStruct._slots:
        #        #
        #        # Select the keyword to read.
        #        try:
        #            abikey = _abipy2abinit[key]
        #        except KeyError:
        #            abikey = key
        #        d[key] = self[abikey]
        #    self.crystal = CrystalStruct(**d)

        # Build pymatgen structure.
        #lattice = Bohr2Ang(self["rprimd"])

        #red_coords = self["xred"]
        #natom = len(red_coords)

        #znucl_type = ncdata.read_value("atomic_numbers")

        ## type_atom[0:natom] --> index Between 1 and number of atom species
        #type_atom = ncdata.read_value("atom_species")

        ## Fortran to C index and float --> int conversion.
        #species = natom * [None]
        #for atom in range(natom):
        #    type_idx = type_atom[atom] - 1
        #    species[atom] = int(znucl_type[type_idx])

        #self.structure = Structure(lattice, species, red_coords)

    #def __str__(self):
    #  return pprint(self)

    def export_crystal(self, filename):
        return self.structure.export(filename)

    def show_bz(self):
        return self.structure.show_bz()

##########################################################################################


class Datasets(collections.Iterable):
    """List of :class:`Dataset`."""

    def __init__(self, filename):
        """
        Args:
            filename: Name of the ABINIT input file.
        """
        # Call the abinit parser.
        dts = Dtsets(filename)

        # Number of datasets (not counting the default one).
        self.ndtset = dts.get_ndtset()

        # Build tuple of Dataset instances.
        self.dtsets = [None] * (self.ndtset+1)

        for idx in range(self.ndtset+1):
            kwargs = {}
            for key in abinit_vars:
                kwargs[key] = dts.get(key, idx)
            self.dtsets[idx] = Dataset(idx, **kwargs)

        self.dtsets = tuple(self.dtsets)

    def __len__(self):
        return self.ndtset + 1

    def __iter__(self):
        return self.dtsets.__iter__()

    def __getitem__(self, slice):
        return self.dtsets[slice]

    #def non_default(self, dt_idx, key)

    def export_json(self, filename):
        with open (filename, 'w') as fh:
            fh.write(json.dumps(self.dtsets, cls=NumPyArangeEncoder))

