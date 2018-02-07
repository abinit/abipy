"""Tests for electrons.bse module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import itertools
import abipy.data as abidata

from abipy import abilab
from abipy.electrons.fatbands import FatBandsFile
from abipy.core.testing import AbipyTest


class TestElectronFatbands(AbipyTest):

    def test_MgB2_fatbands(self):
        """Testing MgB2 fatbands with prtdos 3."""
        fbnc_kpath = FatBandsFile(abidata.ref_file("mgb2_kpath_FATBANDS.nc"))
        repr(fbnc_kpath); str(fbnc_kpath)
        assert fbnc_kpath.to_string(verbose=2)
        assert fbnc_kpath.ebands.kpoints.is_path
        assert not fbnc_kpath.ebands.kpoints.is_ibz
        assert fbnc_kpath.prtdos == 3
        assert fbnc_kpath.prtdosm == 0
        assert fbnc_kpath.mbesslang == 5
        assert fbnc_kpath.pawprtdos == 0
        assert fbnc_kpath.usepaw == 0
        assert fbnc_kpath.nsppol == 1
        assert fbnc_kpath.nspden == 1
        assert fbnc_kpath.nspinor == 1
        assert fbnc_kpath.nkpt == 78
        assert fbnc_kpath.mband == 8
        assert fbnc_kpath.natsph_extra == 0
        assert not fbnc_kpath.ebands.has_metallic_scheme
        assert fbnc_kpath.params["nkpt"] == 78

        if self.has_matplotlib():
            assert fbnc_kpath.plot_fatbands_typeview(tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_lview(tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_siteview(e0=None, view="inequivalent",
                                                     fact=2.0, bdist=[1, 2, 3], show=False)

        if self.has_nbformat():
            fbnc_kpath.write_notebook(nbpath=self.get_tmpname(text=True))

        fbnc_kmesh = FatBandsFile(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc"))
        repr(fbnc_kmesh); str(fbnc_kmesh)
        assert fbnc_kmesh.ebands.kpoints.is_ibz
        assert fbnc_kmesh.ebands.has_metallic_scheme

        if self.has_matplotlib():
            assert fbnc_kmesh.plot_pjdos_typeview(tight_layout=True, show=False)
            assert fbnc_kmesh.plot_pjdos_lview(tight_layout=True, stacked=True, show=False)
            assert fbnc_kmesh.plot_pjdos_lview(tight_layout=True, stacked=False, show=False)
            assert fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, view="type", tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, view="lview", tight_layout=True, show=False)
            assert fbnc_kpath.plot_spilling(e0=0, tight_layout=True, show=False)

            # These calls should return None
            #assert fbnc_kpath.plot_fatbands_spinor() is None
            assert fbnc_kpath.plot_pjdos_lview() is None
            assert fbnc_kpath.plot_pawdos_terms() is None
            #assert fbnc_kpath.plot_pjdos_spinor() is None
            assert fbnc_kpath.plot_fatbands_mview(iatom=0) is None

        fbnc_kpath.close()
        fbnc_kmesh.close()

    def test_nickel_fatbands(self):
        """Testing Nickel fatbands with nsppol = 2 and prtdos 3."""
        lmax = 2
        elims = [-10, 2]

        fbnc_kpath = abilab.abiopen(abidata.ref_file("ni_kpath_FATBANDS.nc"))
        repr(fbnc_kpath); str(fbnc_kpath)
        assert fbnc_kpath.ebands.kpoints.is_path
        assert not fbnc_kpath.ebands.kpoints.is_ibz
        assert fbnc_kpath.prtdos == 3
        assert fbnc_kpath.prtdosm == 1

        assert fbnc_kpath.mbesslang == 5
        assert fbnc_kpath.pawprtdos == 0
        assert fbnc_kpath.usepaw == 0
        assert fbnc_kpath.nsppol == 2
        assert fbnc_kpath.nspden == 2
        assert fbnc_kpath.nspinor == 1
        #assert fbnc_kpath.nkpt == 78
        #assert fbnc_kpath.mband == 8
        assert fbnc_kpath.natsph_extra == 0
        #assert not fbnc_kpath.ebands.has_metallic_scheme

        if self.has_matplotlib():
            assert fbnc_kpath.plot_fatbands_typeview(ylims=elims, lmax=lmax, tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_lview(ylims=elims, lmax=lmax, tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_siteview(e0=None, view="inequivalent",
                                                     fact=2.0, bdist=[1, 2, 3], show=False)

        if self.has_nbformat():
            fbnc_kpath.write_notebook(nbpath=self.get_tmpname(text=True))

        fbnc_kmesh = abilab.abiopen(abidata.ref_file("ni_666k_FATBANDS.nc"))
        repr(fbnc_kmesh); str(fbnc_kmesh)
        assert fbnc_kmesh.ebands.kpoints.is_ibz
        assert fbnc_kmesh.ebands.has_metallic_scheme
        assert fbnc_kmesh.prtdosm == 0

        if self.has_matplotlib():
            for combined_spins, stacked in itertools.product([True, False], [True, False]):
                assert fbnc_kmesh.plot_pjdos_typeview(xlims=elims, combined_spins=combined_spins,
                                                      stacked=stacked, tight_layout=True, show=False)

                assert fbnc_kmesh.plot_pjdos_lview(xlims=elims, combined_spins=combined_spins,
                                                   stacked=stacked, tight_layout=True, show=False)

            assert fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, ylims=elims, lmax=lmax,
                                                       view="type", tight_layout=True, show=False)
            assert fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, ylims=elims, lmax=lmax,
                                                       view="lview", tight_layout=True, show=False)

            # These calls should return None
            #assert fbnc_kpath.plot_fatbands_spinor() is None
            assert fbnc_kpath.plot_pjdos_lview() is None
            assert fbnc_kpath.plot_pawdos_terms() is None
            #assert fbnc_kpath.plot_pjdos_spinor() is None
            # TODO
            assert fbnc_kmesh.plot_fatbands_mview(iatom=0, show=False) is None
            assert fbnc_kpath.plot_fatbands_mview(iatom=0, lmax=2, fact=1.5,
                ylims=[-1.5, +1.0], blist=list(range(4, 10)), show=False) is not None

        fbnc_kpath.close()
        fbnc_kmesh.close()
