"""Tests for core.skw module"""
from __future__ import print_function, division, unicode_literals

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.skw import SkwInterpolator


class TestSkwInterpolator(AbipyTest):
    """Unit tests for SkwInterpolator."""

    def test_silicon_interpolation(self):
        """Testing interpolation of Si band energies with SKW method."""

        from abipy.abilab import abiopen
        with abiopen(abidata.ref_file("si_scf_GSR.nc")) as gsr:
            # Extract data from GSR
            structure, ebands = gsr.structure, gsr.ebands
            kcoords = [k.frac_coords for k in ebands.kpoints]
            cell = (structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers)

            # Get FM part of symmetry operations.
            abispg = structure.abi_spacegroup
            fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]

            # Build interpolator.
            lpratio, has_timrev = 5.0, True
            skw = SkwInterpolator(lpratio, kcoords, ebands.eigens, ebands.fermie, ebands.nelect, cell,
                                  fm_symrel, has_timrev, filter_params=None, verbose=1)

        repr(skw); print(skw)
        assert skw.occtype == "insulator"
        assert skw.use_cache
        assert skw.lpratio == lpratio
        assert skw.nsppol == 1 and skw.nkpt == 29 and skw.nband == 8
        assert skw.ptg_nsym == 48
        assert skw.nr == 145 and skw.rcut is None and skw.rsigma is None
        self.assert_almost_equal(skw.mae, 7.0e-11)
        assert skw.val_ib == 3 and isinstance(skw.val_ib, int)

        kmesh, is_shift = [8, 8, 8], None

        # Test interpolation routines (low-level API).
        k = skw.get_sampling(kmesh, is_shift)
        assert np.all(k.mesh == kmesh) and np.all(k.shift == 0)
        assert k.nbz == 8 ** 3 and k.nibz == 29 and k.weights.sum() == 1
        assert len(k.bz2ibz) == k.nbz

        # interpolate energies at 3 new k-points
        new_kcoords = [(0, 0, 0), (0.1, 0, 0), (0.12, 0.13, 0.14)]
        new_eigens = skw.interp_kpts(new_kcoords).eigens
        assert new_eigens.shape == (skw.nsppol, len(new_kcoords), skw.nband)

        res1 = skw.interp_kpts(new_kcoords, dk1=True, dk2=False)
        print(res1.dedk)
        assert res1.dedk.shape == (skw.nsppol, len(new_kcoords), skw.nband, 3)
        # Group velocities at Gamma should be zero by symmetry.
        self.assert_almost_equal(res1.dedk[0, 0], 0.0)
        #assert 0
        #res12 = skw.interp_kpts(new_kcoords, dk1=True, dk2=True)
        #print(res12.dedk2)

        # Test interpolation routines (high-level API).
        edos = skw.get_edos(kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None)
        jdos = skw.get_jdos_q0(kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None)
        #nest = skw.get_nesting_at_e0(qpoints, kmesh, e0, width=0.2, is_shift=None)

        # Test pickle
        tmpname = self.get_tmpname(text=True)
        skw.pickle_dump(tmpname)
        new = SkwInterpolator.pickle_load(tmpname)

        # Test plotting API.
        if self.has_matplotlib():
            kmeshes = [[2, 2, 2], [4, 4, 4]]
            widths = [0.2, 0.3]
            is_shift = None
            assert skw.plot_dos_vs_kmeshes(kmeshes, is_shift=is_shift, show=False)
            assert skw.plot_jdosq0_vs_kmeshes(kmeshes, is_shift=is_shift, show=False)

            #assert skw.plot_nesting_vs_widths(widths=widths, kmesh=[4, 4, 4], e0=None, qvertices_names=None,
            #                                  line_density=5, is_shift=is_shift, show=False)

            #assert skw.plot_nesting_vs_kmeshes(width=0.2, kmeshes=kmeshes, e0=None, qvertices_names=None, line_density=20,
            #                                   is_shift=is_shift, show=False):

            #assert skw.plot_nesting_vs_kmeshes(width=0.2, kmeshes=kmesh, e0=None, qvertices_names=None, line_density=20,
            #                                   is_shift=is_shift, show=False):

            #assert skw.plot_group_velocites(kvertices_names=None, line_density=20, ax=None, show=False)
