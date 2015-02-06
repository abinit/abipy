# coding: utf-8
"""
This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and abipy objects using
Materials Project data. 
See also pymatgen.matproj subpackage.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
https://www.materialsproject.org/dashboard.
"""
from __future__ import division, print_function, unicode_literals, division

from pymatgen.matproj.rest import MPRester, MPRestError
from abipy.core import Structure


class MPConnection(object):
    """
    Connection to the materials project database.
    """
    def __init__(self, mp_key=None):
        if mp_key is None:
            try:
                mp_key = os.environ['MP_KEY']
            except OSError:
                mp_key = raw_input("there is no key for accessing the materials projects database\n"
                                   "please provide one. (if you don't have one visit the materials \n"
                                   "project website to generate one) :")
        if len(mp_key) > 4:
            self.mp_key = mp_key

    def abistructure_from_mp(self, mpid, final=True):
        """
        Get an abipy :class:`Structure` corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string, e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.

        Returns:
            :class:`Structure` object.
        """
        with MPRester(self.mp_key) as mp_database:
            # Get pytmatgen structure and convert it to abipy structure
            structure = mp_database.get_structure_by_material_id(mpid, final=final)
            structure.__class__ = Structure
            return structure
