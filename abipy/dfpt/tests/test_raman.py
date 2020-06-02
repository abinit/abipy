"""Tests for Raman module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.raman import Raman


class RamanTest(AbipyTest):

    def test_raman(self):
        """Testsing Raman."""

        r = Raman.from_file(abidata.ref_file("AlAs_nl_dte_anaddb.nc"))

        im = r.get_modes_intensities(temp=300, laser_freq=2.54, non_anal_dir=0)
        self.assertArrayEqual(im.shape, (6, 3, 3))
        self.assertAlmostEqual(im[3, 0, 1], 0.03668238755672564)

        il = r.get_lorentz_intensity(temp=300, laser_freq=2.54, non_anal_dir=None, width=0.001, num=100)
        self.assertAlmostEqual(il[0][1].values[50], 11.558435430746329)

        il = r.get_lorentz_intensity(temp=300, laser_freq=20491, non_anal_dir=None, width=5,
                                     pol_in="x", pol_out="y", units="cm-1", relative=True)
        self.assertAlmostEqual(il.values[501], 0.9991991198326963)

        pi = r.get_powder_intensity(temp=300, laser_freq=2.54, relative=True)
        self.assertArrayEqual(pi.perp.shape, (6,))
        self.assertAlmostEqual(pi.tot[5], 1)

        pil = r.get_powder_lorentz_intensity(temp=300, laser_freq=2.54, non_anal_dir=0, width=0.001, num=100)
        self.assertAlmostEqual(pil.tot.values[50], 109.47538037690873)

        if self.has_matplotlib():
            assert r.plot_intensity(temp=300, laser_freq=20491, non_anal_dir=None, width=5,
                                    value="powder", units="cm-1", relative=True, show=False)
            assert r.plot_intensity(temp=300, laser_freq=2.54, non_anal_dir=0, width=0.0001,
                                    value="xy", units="eV", relative=False, show=False, plot_phfreqs=True)

            assert r.plot_intensity(temp=300, laser_freq=2.54, non_anal_dir=None, width=None,
                                    value="xz", units="eV", relative=False, show=False)

    def test_fail_parse(self):
        with self.assertRaises(ValueError):
            Raman.from_file(abidata.ref_file("alas_anaddb.nc"))
