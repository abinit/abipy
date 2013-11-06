from __future__ import print_function, division

import abipy.data as data

from abipy.core.testing import *
from abipy.waves import WFK_File


class TestWFKFile(AbipyTest):
    """Unit tests for WFKFile."""

    def test_read_wfkfile(self):
        """Read WFK_File and waves from NC example data files."""
        assert data.WFK_NCFILES

        for path in data.WFK_NCFILES:
            wfk = WFK_File(path)
            print(wfk)

            spin, kpoint, band = (0, 0, 0)
            structure = wfk.structure
            wave = wfk.get_wave(spin, kpoint, band)

            self.assertTrue(wave == wave)

            #for ig, (g, ug) in enumerate(wave):
            #    self.assertTrue(np.all(g == wave.gvecs[ig]))
            #    self.assertTrue(np.all(ug == wave.ug[:,ig]))

            #print(wave[0:1])

            # Test the norm
            for space in ["g", "r"]:
                norm2 = wave.norm2(space=space)
                if space == "r": norm2 = norm2 / structure.volume
                self.assert_almost_equal(norm2, 1.0)

            # Test bracket with the same wave.
            for space in ["g", "gsphere", "r"]:
                print(space)
                norm2 = wave.braket(wave, space=space)
                if space == "r": norm2 = norm2 / structure.volume
                self.assert_almost_equal(norm2, 1.0)

            # FFT and FFT^{-1} on the BOX.
            ug_mesh = wave.mesh.fft_r2g(wave.ur)
            same_ur = wave.mesh.fft_g2r(ug_mesh)

            self.assert_almost_equal(wave.ur, same_ur)

            # Back to the sphere
            same_ug = wave.gsphere.fromfftmesh(wave.mesh, ug_mesh)
            self.assert_almost_equal(wave.ug, same_ug)

            if self.which("xcrysden") is not None:
                wave.export_ur2(".xsf", structure)


if __name__ == "__main__":
   import unittest
   unittest.main()
