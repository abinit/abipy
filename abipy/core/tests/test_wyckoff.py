"""Tests for core.density module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import unittest

from collections import OrderedDict
from abipy.core.testing import AbipyTest
from abipy.core.structure import Lattice
from abipy.core.wyckoff import Wyckoff


class TestWyckoff(AbipyTest):
    """Unit tests for Wyckoff."""

    def test_sio2(self):
        """Testing Wyckoff positions with SiO2 (requires sympy)"""
        try:
            import sympy
        except ImportError:
            raise unittest.SkipTest("sympy is not installed")

        # http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?gnum=152
        site2wpos = {}
        #3	b	.2.
        #site2wpos["Si"] = "(x,0,5/6)  (0,x,1/6)   (-x,-x,1/2)"
        #3	a	.2.
        site2wpos["Si"] = "[x,0,1/3], (0,x,2/3), (-x,-x,0)"
        #6	c	1
        site2wpos["O"] = """(x,y,z), (-y,x-y,z+1/3), (-x+y,-x,z+2/3), (y,x,-z),
                            (x-y,-y,-z+2/3), (-x,-x+y,-z+1/3)"""

        wyckoff = Wyckoff(site2wpos)
        repr(wyckoff); str(wyckoff)

        site2params = OrderedDict([
            ("Si", dict(x=0.4763)),
            ("O", dict(x=0.1588, y=0.7439, z=0.4612)),
            #("O", dict(x=0.1588, y=0.7439, z=0.4612, yyy=3)),
        ])

        lattice = Lattice.hexagonal(a=4.971, c=5.473)
        for to_unit_cell in [False]:
        #for to_unit_cell in [False, True]:
            structure = wyckoff.generate_structure(lattice, site2params, to_unit_cell=to_unit_cell)
            print("to_unit_cell %s\n" % to_unit_cell, structure)

            #from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            #spga = SpacegroupAnalyzer(structure)
            #structure = spga.get_refined_structure()
            #print("Refined:", structure)

            wyck_params = wyckoff.find_params(structure)
            #structure.perturb(distance=0.01)
            #wyck_params = wyckoff.find_params(structure)

            print("Error of params")
            wyckoff.error_of_params(wyck_params, structure)

            for site, params in site2params.items():
                print(wyck_params[site])
                assert wyck_params[site] == params
