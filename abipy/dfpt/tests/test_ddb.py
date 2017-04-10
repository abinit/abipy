"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import *
from abipy.dfpt.ddb import DdbFile

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class DdbTest(AbipyTest):

    def test_alas_ddb_1qpt_phonons(self):
        """Testing DDB with one q-point"""
        ddb_fname = os.path.join(test_dir, "AlAs_1qpt_DDB")

        with DdbFile(ddb_fname) as ddb:
            print(ddb)
            # Test qpoints.
            assert np.all(ddb.qpoints[0] == [0.25, 0, 0])
            assert len(ddb.qpoints) == 1

            # Test header
            h = ddb.header
            assert h.version == 100401 and h.ecut == 3
            assert "ecut" in h and h["ecut"] == h.ecut
            assert h.occ == 4 * [2]
            assert h.xred.shape == (h.natom, 3) and h.kpt.shape == (h.nkpt, 3)
            print(h.znucl)
            assert "ecut" in ddb.params

            assert np.all(h.symrel[1].T.ravel() == [0, -1, 1, 0, -1, 0, 1, -1, 0])
            assert np.all(h.symrel[2].T.ravel() == [-1, 0, 0, -1, 0, 1, -1, 1, 0])

            # Test structure
            struct = ddb.structure
            assert struct.formula == "Al1 As1"

            # Test interface with Anaddb.
            print(ddb.qpoints[0])
            assert ddb.qindex(ddb.qpoints[0]) == 0

            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=ddb.qpoints[0], verbose=1)
            assert phbands is not None and hasattr(phbands, "phfreqs")

            # Wrong qpoint
            with self.assertRaises(ValueError):
                ddb.anaget_phmodes_at_qpoint(qpoint=(0,0,0), verbose=1)

            # Wrong ngqpt
            with self.assertRaises(ddb.AnaddbError):
                ddb.anaget_phbst_and_phdos_files(ngqpt=(4,4,4), verbose=1)

            # Cannot compute DOS since we need a mesh.
            with self.assertRaises(ddb.AnaddbError):
                ddb.anaget_phbst_and_phdos_files(verbose=1)

            # Test notebook
            if self.has_nbformat():
                ddb.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_alas_ddb_444_nobecs(self):
        """Testing DDB for AlAs on a 4x4x4x q-mesh without Born effective charges."""
        ddb = DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB"))
        print(ddb)
        print(ddb.header)

        ref_qpoints = np.reshape([
                 0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
                 2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
                 5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
                 2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
                 5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
                -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
                 5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
                -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
        ], (-1, 3))

        assert len(ddb.qpoints) == 8
        for qpt, ref_qpt in zip(ddb.qpoints, ref_qpoints):
            assert qpt == ref_qpt

        for qpoint in ddb.qpoints:
            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=qpoint, verbose=1)
            assert phbands is not None and hasattr(phbands, "phfreqs")

        assert np.all(ddb.guessed_ngqpt == [4, 4, 4])

        # Get bands and Dos
        phbands_file, phdos_file = ddb.anaget_phbst_and_phdos_files(verbose=1)
        phbands, phdos = phbands_file.phbands, phdos_file.phdos

        if self.has_matplotlib():
            phbands.plot_with_phdos(phdos, title="Phonon bands and DOS of %s" % phbands.structure.formula, show=False)

        # Get emacro and becs
        emacro, becs = ddb.anaget_emacro_and_becs(chneut=1, verbose=1)
        assert np.all(becs.values == 0)
        assert np.all(becs.becs == 0)
        print(becs)

        self.assert_almost_equal(phdos.idos.values[-1], 3 * len(ddb.structure), decimal=1)
        phbands_file.close()
        phdos_file.close()

        # Test DOS computation via anaddb.
        c = ddb.anacompare_phdos(nqsmalls=[2, 4, 6], num_cpus=None)
        if self.has_matplotlib():
            c.plotter.plot(show=False)

        # Execute anaddb to compute the interatomic forces.
        ifc = ddb.anaget_ifc()
        assert ifc.structure == ddb.structure
        assert ifc.number_of_atoms == len(ddb.structure)

        if self.has_matplotlib():
            ifc.plot_longitudinal_ifc(show=False)
            ifc.plot_longitudinal_ifc_short_range(show=False)
            ifc.plot_longitudinal_ifc_ewald(show=False)

        ddb.close()
