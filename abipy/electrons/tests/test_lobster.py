"""Tests for electrons.lobster module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.abilab import abiopen, LobsterAnalyzer
from abipy.electrons.lobster import LobsterInput

lobster_gaas_dir = os.path.join(abidata.dirpath, "refs", "lobster_gaas")

class CoxpTest(AbipyTest):

    def test_coxp(self):
        """Test files based on the GaAs lobster test. Extracted from abinit calculation."""

        # Test COHPCAR
        with abiopen(os.path.join(lobster_gaas_dir, "GaAs_COHPCAR.lobster.gz")) as cohp:
            repr(cohp); str(cohp)
            assert cohp.to_string(verbose=2)
            assert cohp.nsppol == 1
            assert len(cohp.type_of_index) == 2
            assert cohp.type_of_index[0] == "Ga" and cohp.type_of_index[1] == "As"
            assert cohp.cop_type == "cohp"
            self.assert_almost_equal(cohp.fermie, 2.29843, decimal=4)

            assert len(cohp.energies) == 401
            assert (0, 1) in cohp.partial
            assert ("4s", "4p_x") in cohp.partial[(0, 1)]
            assert 0 in cohp.partial[(0, 1)][("4s", "4p_x")]
            assert len(cohp.partial[(0, 1)][("4s", "4p_x")][0]["single"]) == 401

            self.assertAlmostEqual(cohp.partial[(0, 1)][("4s", "4p_x")][0]["single"][200], -0.02075)
            assert len(cohp.site_pairs_total) == 2
            assert len(cohp.site_pairs_partial) == 2

            self.assertAlmostEqual(cohp.functions_pair_lorbitals[(0, 1)][("4s", "4p")][0].values[200], -0.06225)
            self.assertAlmostEqual(cohp.functions_pair_morbitals[(0, 1)][("4s", "4p_x")][0].values[200], -0.02075)
            self.assertAlmostEqual(cohp.functions_pair[(0, 1)][0].values[200], -0.06124)
            #self.check_average(cohp)

            if self.has_matplotlib():
                assert cohp.plot(title="default values", show=False)
                assert cohp.plot_site_pairs_total(from_site_index=[0, 1], what="single", exchange_xy=True, show=False)
                assert cohp.plot_site_pairs_partial(from_site_index=[0, 1], what="single", exchange_xy=True, show=False)
                assert cohp.plot_average_pairs(with_site_index=[0, 1], what="single", exchange_xy=True, show=False)
                #assert cohp.plot_with_ebands(ebands_kpath, show=False)

            if self.has_nbformat():
                assert cohp.write_notebook(nbpath=self.get_tmpname(text=True))

        # Test COOPCAR
        with abiopen(os.path.join(lobster_gaas_dir, "GaAs_COOPCAR.lobster.gz")) as coop:
            repr(coop); str(coop)
            assert coop.to_string(verbose=2)
            assert coop.nsppol == 1
            assert len(coop.type_of_index) == 2
            assert coop.type_of_index[0] == "Ga" and coop.type_of_index[1] == "As"
            assert coop.cop_type == "coop"
            assert len(coop.energies) == 401
            assert (0, 1) in coop.partial
            assert ("4s", "4p_x") in coop.partial[(0, 1)]
            assert 0 in coop.partial[(0, 1)][("4s", "4p_x")]
            assert len(coop.partial[(0, 1)][("4s", "4p_x")][0]["single"]) ==  401
            self.assertAlmostEqual(coop.partial[(0, 1)][("4s", "4p_x")][0]["single"][200], 0.01466)
            #self.check_average(coop)

            if self.has_matplotlib():
                assert coop.plot(title="default values", show=False)

            if self.has_nbformat():
                assert coop.write_notebook(nbpath=self.get_tmpname(text=True))

    def check_average(self, coxp):
        # averaged should contain the average over all atom pairs.
        # pair data is stored in total[pair][spin][what]
        for what in ("single",): #, "integrated")
            sum_all = np.zeros((coxp.nsppol, len(coxp.energies)))
            for pair, pair_data in coxp.total.items():
                for spin in range(coxp.nsppol):
                    sum_all[spin] += pair_data[spin][what]
            #assert len(coxp.total) == 2
            npairs = len(coxp.total)
            sum_all /= (npairs)
            #sum_all /= (2.0 * npairs)

        for spin in range(coxp.nsppol):
            ref_values = coxp.averaged[spin][what] #* len(coxp.total)
            vals = sum_all[spin]
            # In this range, we loose precision.
            vals = np.where(np.abs(vals) > 1e-3, vals, ref_values)
            #print("ref_values")
            #print(ref_values)
            for x, y in zip(ref_values, vals):
                print(x, y, y/x if abs(x) > 0 else 0)
            self.assert_almost_equal(ref_values, vals, decimal=3)


class ICoxpTest(AbipyTest):

    def test_icoxp(self):
        """Testing ICOHPLIST files."""
        with abiopen(os.path.join(lobster_gaas_dir, "GaAs_ICOHPLIST.lobster.gz")) as icohp:
            repr(icohp); str(icohp)
            assert icohp.to_string(verbose=2)
            assert (0, 1) in icohp.values
            assert 0 in icohp.values[(0, 1)]
            self.assertAlmostEqual(icohp.values[(0, 1)][0]['average'], -4.36062)
            self.assertAlmostEqual(icohp.dataframe.average[0], -4.36062)

            assert icohp.cop_type == "cohp"
            assert len(icohp.type_of_index) == 2
            assert icohp.type_of_index[0] == "Ga" and icohp.type_of_index[1] == "As"

            if self.has_matplotlib():
                assert icohp.plot(title="default values", show=False)

            if self.has_nbformat():
                assert icohp.write_notebook(nbpath=self.get_tmpname(text=True))


class LobsterDosTest(AbipyTest):

    def test_lobsterdos(self):
        """Testing Lobster DOSCAR files."""
        with abiopen(os.path.join(lobster_gaas_dir, "GaAs_DOSCAR.lobster.gz")) as ldos:
            repr(ldos); str(ldos)
            assert ldos.to_string(verbose=2)
            assert len(ldos.energies) == 401
            assert ldos.nsppol == 1
            assert ldos.nsites == 2
            assert 1 in ldos.pdos
            assert "4p_x" in ldos.pdos[1]
            assert 0 in ldos.pdos[1]["4p_x"]
            self.assert_almost_equal(ldos.fermie, 2.29843, decimal=4)
            assert ldos.type_of_index[0] == "Ga" and ldos.type_of_index[1] == "As"

            self.assertAlmostEqual(ldos.pdos[1]["4p_x"][0][200], 0.02694)
            self.assertAlmostEqual(ldos.total_dos[0][200], 0.17824)
            self.check_average(ldos)

            if self.has_matplotlib():
                assert ldos.plot(title="default values", show=False)
                assert ldos.plot_pdos_site(site_index=0, title="default values", show=False)
                assert ldos.plot_pdos_site(site_index=[0, 1], exchange_xy=True, title="default values", show=False)
                #assert cohp.plot_with_ebands(ebands_kpath, show=False)

            if self.has_nbformat():
                assert ldos.write_notebook(nbpath=self.get_tmpname(text=True))

    def check_average(self, icoxp):
        # total_dos should contain the sum over all atoms and orbitals.
        # self.pdos[site_index]["4p_x"][spin]
        sum_all = np.zeros((icoxp.nsppol, len(icoxp.energies)))
        for isite, s in icoxp.pdos.items():
            for orb, d in s.items():
                for spin in range(icoxp.nsppol):
                    sum_all[spin] += d[spin]

        for spin in range(icoxp.nsppol):
            ref_values = icoxp.total_dos[spin]
            vals = sum_all[spin]
            # In this range, we loose precision.
            vals = np.where(np.abs(vals) > 1e-3, vals, ref_values)
            #for x, y in zip(ref_values, vals):
            #    print(x, y, y/x if abs(x) > 0 else 0)
            self.assert_almost_equal(ref_values, vals, decimal=3)


class LobsterAnalyzerTest(AbipyTest):

    def test_lobster_analyzer(self):
        """Testing lobster analyzer."""
        lobana = LobsterAnalyzer.from_dir(lobster_gaas_dir, prefix="GaAs_")
        repr(lobana); str(lobana)
        assert lobana.nsppol == 1
        assert lobana.to_string(verbose=2)

        if self.has_matplotlib():
            assert lobana.plot(title="default values", show=False)
            assert lobana.plot_coxp_with_dos(from_site_index=1, what="coop", title="default values", show=False)
            assert lobana.plot_coxp_with_dos(from_site_index=[1, 2], with_orbitals=True, show=False)

        if self.has_nbformat():
            assert lobana.write_notebook(nbpath=self.get_tmpname(text=True))


class LobsterInputTest(AbipyTest):

    def test_lobsterinput(self):
        """Testing LobsterInput."""
        lin = LobsterInput(atom_pairs=[(0, 1), (2, 3)], orbitalwise=True)
        repr(lin); str(lin)
        assert lin.to_string(verbose=2)
        lin.set_basis_functions_from_abinit_pseudos([abidata.pseudo("Al.GGA_PBE-JTH.xml")])
        # XML parser does not guarantee order.
        assert sorted(lin.basis_functions[0].split()) == sorted("Al 3s 3p".split())
        lin.write(self.mkdtemp())
