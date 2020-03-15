"""Tests for electrons.effmass_analyzer module"""
import abipy.data as abidata

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

        emana.select_kpoint_band((0, 0, 0), band=3, spin=0, etol_ev=0.1)
        emana.summarize()

        emana.select_band_edges()
        emana.select_cbm()
        emana.select_vbm(etol_ev=1e-3)
        emana.summarize()

        segment = emana.segments[0]
        repr(segment); str(segment)
        assert segment.to_string(verbose=2)
        df = segment.get_dataframe_with_accuracies(acc_list=(2, 4))

        #assert len(emana.segments) == 1
        #for segment in emana.segments[0]:
        #    segment.get_effmass_line(acc=2)

        if self.has_matplotlib():
            assert emana.plot_emass(acc=4, show=False)
            assert emana.plot_all_segments(show=False)
            assert emana.segments[0].plot_emass(show=False)
