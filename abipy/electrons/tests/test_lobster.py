"""Tests for electrons.lobster module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata
from abipy.electrons.lobster import LobsterInput, Coxp, ICoxp, LobsterDos
from abipy.core.testing import AbipyTest


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class CoxpTest(AbipyTest):

    def test_coxp(self):
        # Test file based on the GaAs lobster test. Extracted from abinit calculation
        cohp = Coxp.from_file(os.path.join(test_dir, "GaAs_COHPCAR.lobster.gz"))
        self.assertEqual(len(cohp.energies), 401)
        self.assertIn((0,1), cohp.partial)
        self.assertIn(("4s", "4p_x"), cohp.partial[(0,1)])
        self.assertIn(0, cohp.partial[(0,1)][("4s", "4p_x")])
        self.assertEqual(len(cohp.partial[(0,1)][("4s", "4p_x")][0]["single"]), 401)

        self.assertAlmostEqual(cohp.partial[(0,1)][("4s", "4p_x")][0]["single"][200], -0.02075)
        self.assertEqual(len(cohp.site_pairs_total), 2)
        self.assertEqual(len(cohp.site_pairs_partial), 2)

        self.assertAlmostEqual(cohp.functions_pair_lorbitals[(0,1)][("4s", "4p")][0].values[200], -0.06225)
        self.assertAlmostEqual(cohp.functions_pair_morbitals[(0,1)][("4s", "4p_x")][0].values[200], -0.02075)
        self.assertAlmostEqual(cohp.functions_pair[(0,1)][0].values[200], -0.06124)

        coop = Coxp.from_file(os.path.join(test_dir, "GaAs_COOPCAR.lobster.gz"))
        self.assertEqual(len(coop.energies), 401)
        self.assertIn((0, 1), coop.partial)
        self.assertIn(("4s", "4p_x"), coop.partial[(0, 1)])
        self.assertIn(0, coop.partial[(0, 1)][("4s", "4p_x")])
        self.assertEqual(len(coop.partial[(0, 1)][("4s", "4p_x")][0]["single"]), 401)

        self.assertAlmostEqual(coop.partial[(0, 1)][("4s", "4p_x")][0]["single"][200], 0.01466)


class ICoxpTest(AbipyTest):

    def test_icoxp(self):
        icohp = ICoxp.from_file(os.path.join(test_dir, "GaAs_ICOHPLIST.lobster.gz"))
        self.assertIn((0,1), icohp.values)
        self.assertIn(0, icohp.values[(0,1)])

        self.assertAlmostEqual(icohp.values[(0,1)][0]['average'], -4.36062)


class LobsterDosTest(AbipyTest):

    def test_lobsterdos(self):
        ldos = LobsterDos.from_file(os.path.join(test_dir, "GaAs_DOSCAR.lobster.gz"))
        self.assertEqual(len(ldos.energies), 401)
        self.assertIn(1, ldos.pdos)
        self.assertIn("4p_x", ldos.pdos[1])
        self.assertIn(0, ldos.pdos[1]["4p_x"])

        self.assertAlmostEqual(ldos.pdos[1]["4p_x"][0][200], 0.02694)
        self.assertAlmostEqual(ldos.total_dos[0][200], 0.17824)


class LobsterInputTest(AbipyTest):

    def test_lobsterinput(self):
        lin = LobsterInput(atom_pairs=[(0,1), (2,3)], orbitalwise=True)
        lin.set_basis_functions_from_abinit_pseudos([abidata.pseudo("Al.GGA_PBE-JTH.xml")])
        assert lin.basis_functions[0] == "Al 3s  3p"

        lin.write(self.mkdtemp())
