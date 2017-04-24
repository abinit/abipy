"""Tests for Wfkfile module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.waves import WfkFile


class TestWFKFile(AbipyTest):
    """Unit tests for WfkFile."""

    def test_read_wfkfile(self):
        """Testing WfkFile and waves from NC example data files."""
        filepath = abidata.ref_file("si_nscf_WFK.nc")

        wfk = WfkFile(filepath)
        repr(wfk); str(wfk)
        assert wfk.nsppol == 1 and wfk.nspinor == 1 and wfk.nspden == 1

        spin, kpoint, band = (0, 0, 0)
        with self.assertRaises(ValueError):
            wfk.get_wave(-1, kpoint, band)

        wave = wfk.get_wave(spin, kpoint, band)
        repr(wave); str(wave)
        assert wave.structure is wfk.structure
        assert wave.shape == (wfk.nspinor, wave.npw)
        other_wave = wfk.get_wave(spin, kpoint, band + 1)

        assert wave == wave
        assert wave != other_wave

        for ig, (g, u_g) in enumerate(wave):
            assert np.all(g == wave.gvecs[ig])
            assert np.all(u_g == wave.ug[:,ig])
            if ig == 5: break
        #print(wave[0:1])

        # Test the norm
        for space in ["g", "r"]:
            norm2 = wave.norm2(space=space)
            if space == "r": norm2 = norm2 / wave.structure.volume
            self.assert_almost_equal(norm2, 1.0)

        # Test bracket with the same wave.
        for space in ["g", "gsphere", "r"]:
            norm2 = wave.braket(wave, space=space)
            if space == "r": norm2 = norm2 / wave.structure.volume
            self.assert_almost_equal(norm2, 1.0)

        ur = wave.get_ur_mesh(wave.mesh, copy=False)
        assert ur is wave.ur
        other_mesh = wave.mesh.__class__((8, 8, 8), wave.mesh.vectors)
        ur = wave.get_ur_mesh(other_mesh)
        assert ur.shape == (8, 8, 8)
        self.assert_almost_equal(other_mesh.integrate(ur.conj() * ur) / wave.structure.volume, 1.0)

        visu = wfk.visualize_ur2(spin=0, kpoint=0, band=0, visu_name="vesta")
        assert callable(visu)

        # FFT and FFT^{-1} on the BOX.
        visu = wave.visualize_ur2("xcrysden")
        assert callable(visu)
        ug_mesh = wave.mesh.fft_r2g(wave.ur)
        same_ur = wave.mesh.fft_g2r(ug_mesh)

        self.assert_almost_equal(wave.ur, same_ur)

        # Back to the sphere
        same_ug = wave.gsphere.fromfftmesh(wave.mesh, ug_mesh)
        self.assert_almost_equal(wave.ug, same_ug)

        wave.set_ug(same_ug)
        self.assert_equal(wave.ug, same_ug)
        with self.assertRaises(ValueError):
            wave.set_ug(np.empty(3))

        wave.export_ur2(".xsf")

        if self.has_matplotlib():
            assert wave.plot_line(0, 1, num=100, show=False)

        if self.has_nbformat():
            wfk.write_notebook(nbpath=self.get_tmpname(text=True))

        wfk.close()
