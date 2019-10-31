"""Tests for electrons.effmass_analyzer module"""
#import numpy as np
import abipy.data as abidata
#from abipy import abilab

from abipy.core.testing import AbipyTest
from abipy.electrons.effmass_analyzer import EffMassAnalyzer


class EffMassAnalyzerTest(AbipyTest):

    def test_api(self):
        """Testing EffMassAnalyzer."""
        emana = EffMassAnalyzer.from_file(abidata.ref_file("si_nscf_GSR.nc"))
        repr(emana); str(emana)
        assert emana.to_string(verbose=2)

        with self.assertRaises(RuntimeError):
            emana.summarize()

        emana.set_kpoint_band((0, 0, 0), band=3, degtol_ev=0.1)
        emana.summarize()

        segment = emana.segments_spin[0][0]
        repr(segment); str(segment)
        assert segment.to_string(verbose=2)
        df = segment.get_dataframe_with_accuracies()

        #assert len(emana.segments_spin) == 1
        #for segment in emana.segments_spin[0]:
        #    segment.get_effmass_line(acc=2)

        if self.has_matplotlib():
            assert emana.plot_emass(acc=4, show=False)
            assert emana.plot_all_segments(show=False)
            assert emana.segments_spin[0][0] .plot_emass(show=False)

        #if self.has_nbformat():
        #   assert emana.write_notebook(nbpath=self.get_tmpname(text=True))

    #def test_ebands_skw_interpolation(self):
    #    """Testing SKW interpolation."""
    #    si_ebands_kmesh = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))

    #    # Test interpolation.
    #    vertices_names = [((0.0, 0.0, 0.0), "G"), ((0.5, 0.5, 0.0), "M")]
    #    r = si_ebands_kmesh.interpolate(lpratio=10, vertices_names=vertices_names,
    #                                    kmesh=[8, 8, 8], verbose=1)
    #    assert r.ebands_kpath is not None
    #    assert r.ebands_kpath.kpoints.is_path
    #    assert not r.ebands_kpath.kpoints.is_ibz
    #    mpdivs, shifts = r.ebands_kpath.kpoints.mpdivs_shifts
    #    assert mpdivs is None and shifts is None

    #    assert r.ebands_kmesh is not None
    #    assert r.ebands_kmesh.kpoints.is_ibz
    #    assert not r.ebands_kmesh.kpoints.is_path
    #    assert r.ebands_kmesh.kpoints.ksampling is not None
    #    assert r.ebands_kmesh.kpoints.is_mpmesh
    #    mpdivs, shifts = r.ebands_kmesh.kpoints.mpdivs_shifts
    #    self.assert_equal(mpdivs, [8, 8, 8])
    #    self.assert_equal(shifts.flatten(), [0, 0, 0])

    #    # Export it in BXSF format.
    #    r.ebands_kmesh.to_bxsf(self.get_tmpname(text=True))

    #def test_derivatives(self):
    #    """Testing computation of effective masses."""
    #    ebands = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))

    #    # Hack eigens to simulate free-electron bands.
    #    # This should produce all(effective masses == 1)
    #    new_eigens = np.empty(ebands.shape)
    #    branch = 0.5 * units.Ha_to_eV * np.array([(k.norm * units.bohr_to_ang)**2 for k in ebands.kpoints])
    #    for spin in ebands.spins:
    #        for band in range(ebands.mband):
    #            new_eigens[spin, :, band] = branch
    #    ebands._eigens = new_eigens

    #    effm_lines = ebands.effective_masses(spin=0, band=0, acc=2)

    #    # Flatten structure (.flatten does not work in this case)
    #    values = []
    #    for arr in effm_lines:
    #        values.extend(arr)

    #    self.assert_almost_equal(np.array(values), 1.0)

    #    em = ebands.get_effmass_line(spin=0, kpoint=(0, 0, 0), band=0)
    #    repr(em); str(em)
    #    #self.assert_almost_equal(np.array(values), 1.0)
