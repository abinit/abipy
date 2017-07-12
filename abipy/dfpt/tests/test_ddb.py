"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import *
from abipy.dfpt.ddb import DdbFile, DielectricTensorGenerator
from abipy.dfpt.anaddbnc import AnaddbNcFile
from abipy.dfpt.phonons import PhononBands
import abipy.data as abidata

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class DdbTest(AbipyTest):

    def test_alas_ddb_1qpt_phonons(self):
        """Testing DDB with one q-point"""
        with DdbFile(os.path.join(test_dir, "AlAs_1qpt_DDB")) as ddb:
            repr(ddb); print(ddb)
            # Test qpoints.
            assert len(ddb.qpoints) == 1
            assert np.all(ddb.qpoints[0] == [0.25, 0, 0])
            assert ddb.natom == len(ddb.structure)

            # Test header
            h = ddb.header
            assert h.version == 100401 and h.ecut == 3
            assert "ecut" in h and h["ecut"] == h.ecut
            assert "ixc" in ddb.params
            assert ddb.params["ixc"] == 7
            assert h.occ == 4 * [2]
            assert h.xred.shape == (h.natom, 3) and h.kpt.shape == (h.nkpt, 3)
            self.assert_equal(h.znucl, [13, 33])
            assert ddb.version == 100401

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
                ddb.anaget_phmodes_at_qpoint(qpoint=(0, 0, 0), verbose=1)

            # Wrong ngqpt
            with self.assertRaises(ddb.AnaddbError):
                try:
                    ddb.anaget_phbst_and_phdos_files(ngqpt=(4, 4, 4), verbose=1)
                except Exception as exc:
                    # This to test AnaddbError.__str__
                    print(exc)
                    raise

            # Cannot compute DOS since we need a mesh.
            with self.assertRaises(ddb.AnaddbError):
                ddb.anaget_phbst_and_phdos_files(verbose=1)

            # Test notebook
            if self.has_nbformat():
                ddb.write_notebook(nbpath=self.get_tmpname(text=True))

            # Test block parsing.
            blocks = ddb._read_blocks()
            assert len(blocks) == 1
            assert blocks[0]["qpt"] == [0.25, 0, 0]

            lines = blocks[0]["data"]
            assert lines[0].rstrip() == " 2nd derivatives (non-stat.)  - # elements :      36"
            assert lines[2].rstrip() ==  "   1   1   1   1  0.80977066582497D+01 -0.46347282336361D-16"
            assert lines[-1].rstrip() == "   3   2   3   2  0.49482344898401D+01 -0.44885664256253D-17"

            for qpt in ddb.qpoints:
                assert ddb.get_block_for_qpoint(qpt)
                assert ddb.get_block_for_qpoint(qpt.frac_coords)

            assert ddb.replace_block_for_qpoint(ddb.qpoints[0], blocks[0]["data"])

            # Write new DDB file.
            tmp_file = nbpath=self.get_tmpname(text=True)
            ddb.write(tmp_file)
            with DdbFile(tmp_file) as new_ddb:
                assert ddb.qpoints == new_ddb.qpoints
                # Call anaddb to check if we can read new DDB
                phbands = new_ddb.anaget_phmodes_at_qpoint(qpoint=new_ddb.qpoints[0], verbose=1)
                assert phbands is not None and hasattr(phbands, "phfreqs")

    def test_alas_ddb_444_nobecs(self):
        """Testing DDB for AlAs on a 4x4x4x q-mesh without Born effective charges."""
        ddb = DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB"))
        repr(ddb); str(ddb)
        print(ddb.header)

        # TODO
        #assert ddb.has_phonon_terms()
        #assert not ddb.has_bec_terms()
        #assert not ddb.has_emacro_terms()

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
            assert phbands.plot_with_phdos(phdos, show=False,
                title="Phonon bands and DOS of %s" % phbands.structure.formula)

        # Get emacro and becs
        emacro, becs = ddb.anaget_emacro_and_becs(chneut=1, verbose=1)
        assert np.all(becs.values == 0)
        assert np.all(becs.values == 0)
        repr(becs); str(becs)
        assert becs.to_string(verbose=1)

        self.assert_almost_equal(phdos.idos.values[-1], 3 * len(ddb.structure), decimal=1)
        phbands_file.close()
        phdos_file.close()

        # Test DOS computation via anaddb.
        c = ddb.anacompare_phdos(nqsmalls=[2, 4, 6], num_cpus=1)
        assert c.phdoses and c.plotter is not None
        if self.has_matplotlib():
            assert c.plotter.combiplot(show=False)

        # Execute anaddb to compute the interatomic forces.
        ifc = ddb.anaget_ifc()
        str(ifc); repr(ifc)
        assert ifc.structure == ddb.structure
        assert ifc.number_of_atoms == len(ddb.structure)

        if self.has_matplotlib():
            assert ifc.plot_longitudinal_ifc(show=False)
            assert ifc.plot_longitudinal_ifc_short_range(show=False)
            assert ifc.plot_longitudinal_ifc_ewald(show=False)

        ddb.close()


class DielectricTensorGeneratorTest(AbipyTest):

    def test_base(self):
        """Testing DielectricTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")
        phbstnc_fname = abidata.ref_file("AlAs_nl_dte_PHBST.nc")

        d = DielectricTensorGenerator.from_files(phbstnc_fname, anaddbnc_fname)
        repr(d); str(d)

        self.assertAlmostEqual(d.tensor_at_frequency(0.001, units='Ha')[0,0], 11.917178775812721)

        d = DielectricTensorGenerator.from_objects(PhononBands.from_file(phbstnc_fname),
                                                   AnaddbNcFile.from_file(anaddbnc_fname))

        self.assertAlmostEqual(d.tensor_at_frequency(0.001, units='Ha')[0,0], 11.917178775812721)

        if self.has_matplotlib():
            assert d.plot_vs_w(0.0001, 0.01, 10, units="Ha", show=False)