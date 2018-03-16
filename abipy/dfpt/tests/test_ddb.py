"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.dfpt.ddb import DdbFile, DielectricTensorGenerator
from abipy.dfpt.anaddbnc import AnaddbNcFile
from abipy.dfpt.phonons import PhononBands


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
            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=ddb.qpoints[0], lo_to_splitting=False, verbose=1)

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
        assert str(ddb.header)
        assert ddb.to_string(verbose=2)
        assert ddb.header["nkpt"] == 256
        assert ddb.header.nsym == 24 and ddb.header.ntypat == 2
        self.assert_equal(ddb.header.znucl, [13, 33])
        self.assert_equal(ddb.header.acell, [1, 1, 1])
        self.assert_equal(ddb.header.ngfft, [10, 10, 10])
        self.assert_equal(ddb.header.spinat, 0.0)
        #assert ddb.header.occ.shape = (ddb.header.nsppol, ddb.header.nkpt, ddb.header.nsppol)

        assert not ddb.has_qpoint([0.345, 0.456, 0.567])
        assert ddb.has_qpoint([0, 0, 0])
        for qpoint in ddb.qpoints:
            assert ddb.has_qpoint(qpoint)
            assert ddb.has_qpoint(qpoint.frac_coords)
            assert qpoint in ddb.computed_dynmat
            assert len(ddb.computed_dynmat[qpoint].index[0]) == 4

        assert ddb.has_bec_terms(select="at_least_one")
        assert not ddb.has_bec_terms(select="all")
        assert not ddb.has_emacro_terms()
        assert not ddb.has_lo_to_data()

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

        assert ddb.view_phononwebsite(verbose=1, dryrun=True) == 0

        if self.has_matplotlib():
            assert phbands.plot_with_phdos(phdos, show=False,
                title="Phonon bands and DOS of %s" % phbands.structure.formula)
            assert phbands_file.plot_phbands(show=False)

        # Get emacro and becs
        emacro, becs = ddb.anaget_emacro_and_becs(chneut=1, verbose=1)
        assert np.all(becs.values == 0)
        #assert np.all(emacro.values == 0)
        repr(becs); str(becs)
        assert becs.to_string(verbose=2)

        self.assert_almost_equal(phdos.idos.values[-1], 3 * len(ddb.structure), decimal=1)
        phbands_file.close()
        phdos_file.close()

        # Test DOS computation via anaddb.
        c = ddb.anacompare_phdos(nqsmalls=[2, 3, 4], dipdip=0, num_cpus=1, verbose=2)
        assert c.phdoses and c.plotter is not None
        if self.has_matplotlib():
            assert c.plotter.combiplot(show=False)

        # Use threads and gaussian DOS.
        c = ddb.anacompare_phdos(nqsmalls=[2, 3, 4], dos_method="gaussian", dipdip=0, asr=0,
                num_cpus=2, verbose=2)
        assert c.phdoses and c.plotter is not None

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

    def test_zno_gamma_ddb_with_becs(self):
        """Testing DDB for ZnO: Gamma only, with Born effective charges and E_macro."""
        with DdbFile(os.path.join(test_dir, "ZnO_gamma_becs_DDB")) as ddb:
            repr(ddb); str(ddb)
            assert str(ddb.header)
            assert ddb.to_string(verbose=2)
            assert ddb.header["nkpt"] == 486
            assert ddb.header.nband == 22 and ddb.header.occopt == 1
            self.assert_equal(ddb.header.typat, [1, 1, 2, 2])
            assert len(ddb.header.wtk) == ddb.header.nkpt
            #assert ddb.header.occ.shape = (ddb.header.nsppol, ddb.header.nkpt, ddb.header.nsppol)

            assert not ddb.has_qpoint([0.345, 0.456, 0.567])
            assert ddb.has_qpoint([0, 0, 0])
            assert len(ddb.qpoints) == 1
            for qpoint in ddb.qpoints:
                assert ddb.has_qpoint(qpoint)
                assert ddb.has_qpoint(qpoint.frac_coords)
                assert qpoint in ddb.computed_dynmat
                assert len(ddb.computed_dynmat[qpoint].index[0]) == 4

            # Test Lru_cache as well
            assert ddb.has_bec_terms(select="at_least_one")
            assert ddb.has_bec_terms(select="at_least_one")
            assert not ddb.has_bec_terms(select="all")
            assert not ddb.has_bec_terms(select="all")
            assert ddb.has_emacro_terms()
            assert ddb.has_lo_to_data()

            # Get emacro and becs
            emacro, becs = ddb.anaget_emacro_and_becs(chneut=1, verbose=1)
            ref_becs_values = [
                [[  2.15646571e+00,   0.00000000e+00,   3.26402110e-25],
                 [  0.00000000e+00,   2.15646571e+00,  -5.46500204e-24],
                 [ -5.66391495e-25,  -6.54012564e-25,   2.19362823e+00]],
                [[  2.15646571e+00,   0.00000000e+00,   1.19680774e-24],
                 [  0.00000000e+00,   2.15646571e+00,   8.10327888e-24],
                 [ -1.69917448e-24,  -1.30802513e-24,   2.19362823e+00]],
                [[ -2.15646571e+00,   6.66133815e-16,  -1.84961196e-24],
                 [  8.88178420e-16,  -2.15646571e+00,   2.82672519e-24],
                 [ -3.39834897e-24,  -3.27006282e-25,  -2.19362823e+00]],
                [[ -2.15646571e+00,  -6.66133815e-16,   3.26402110e-25],
                 [ -8.88178420e-16,  -2.15646571e+00,  -5.46500204e-24],
                 [  5.66391495e-24,   2.28904397e-24,  -2.19362823e+00]]
                ]

            self.assert_almost_equal(becs.values, ref_becs_values)
            #self.assert_almost_equal(emacro.values, ref_emacro_values)
            repr(becs); str(becs)
            assert becs.to_string(verbose=2)

            # get the dielectric tensor generator from anaddb
            dtg = ddb.anaget_dielectric_tensor_generator(verbose=2)
            assert dtg is not None and hasattr(dtg, "phfreqs")

    def test_mgo_becs_emacro(self):
        """
        Testing DDB for MgO with with Born effective charges and E_macro.
        Large breaking of the ASR.
        """
        with abilab.abiopen(abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")) as ddb:
            assert ddb.structure.formula == "Mg1 O1"
            assert len(ddb.qpoints) == 72
            assert ddb.has_emacro_terms()
            assert ddb.has_bec_terms()

            if self.has_matplotlib():
                plotter = ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(0, 1), dipdip=1,
                    nqsmall=2, ndivsm=5, dos_method="tetra", ngqpt=None, verbose=2)
                str(plotter)
                assert plotter.combiplot(show=False)

                # Test nqsmall == 0
                plotter = ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(0, 1), dipdip=1,
                    nqsmall=0, ndivsm=5, dos_method="tetra", ngqpt=None, verbose=2)
                assert plotter.gridplot(show=False)

                plotter = ddb.anacompare_dipdip(chneut_list=(0, 1), asr=1,
                    nqsmall=0, ndivsm=5, dos_method="gaussian", ngqpt=None, verbose=2)
                assert plotter.gridplot(show=False)

    def test_mgb2_ddbs_ngkpt_tsmear(self):
        """Testing multiple DDB files and gridplot_with_hue."""
        paths = [
            #"mgb2_444k_0.01tsmear_DDB",
            #"mgb2_444k_0.02tsmear_DDB",
            #"mgb2_444k_0.04tsmear_DDB",
            "mgb2_888k_0.01tsmear_DDB",
            #"mgb2_888k_0.02tsmear_DDB",
            "mgb2_888k_0.04tsmear_DDB",
            "mgb2_121212k_0.01tsmear_DDB",
            #"mgb2_121212k_0.02tsmear_DDB",
            "mgb2_121212k_0.04tsmear_DDB",
        ]
        paths = [os.path.join(abidata.dirpath, "refs", "mgb2_phonons_nkpt_tsmear", f) for f in paths]

        robot = abilab.DdbRobot.from_files(paths)
        robot.remap_labels(lambda ddb: "nkpt: %s, tsmear: %.3f" % (ddb.header["nkpt"], ddb.header["tsmear"]))

        # Invoke anaddb to get bands and doses
        r = robot.anaget_phonon_plotters(nqsmall=2)

        data = robot.get_dataframe_at_qpoint(qpoint=(0, 0, 0), units="meV", with_geo=False)
        assert "tsmear" in data
        self.assert_equal(data["ixc"].values, 1)

        if self.has_matplotlib():
            assert r.phbands_plotter.gridplot_with_hue("nkpt", with_dos=True, show=False)
            assert r.phbands_plotter.gridplot_with_hue("nkpt", with_dos=False, show=False)

        robot.close()

    def test_ddb_from_mprester(self):
        """Test creation methods for DdbFile and DdbRobot from MP REST API."""
        #ddb = abilab.DdbFile.from_mpid("mp-1138")
        ddb = abilab.DdbFile.from_mpid("mp-149")
        assert ddb.structure.formula == "Si2"
        self.assert_equal(ddb.guessed_ngqpt, [9, 9, 9])
        assert ddb.header["version"] == 100401
        assert ddb.header["ixc"] == -116133

        mpid_list = ["mp-149", "mp-1138"]
        robot = abilab.DdbRobot.from_mpid_list(mpid_list)
        assert len(robot) == len(mpid_list)
        assert robot.abifiles[1].structure.formula == "Li1 F1"
        assert robot.abifiles[1].header["ixc"] == -116133


class DielectricTensorGeneratorTest(AbipyTest):

    def test_base(self):
        """Testing DielectricTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")
        phbstnc_fname = abidata.ref_file("AlAs_nl_dte_PHBST.nc")

        d = DielectricTensorGenerator.from_files(phbstnc_fname, anaddbnc_fname)
        repr(d); str(d)

        self.assertAlmostEqual(d.tensor_at_frequency(0.001, units='Ha')[0,0], 11.917178540635028)

        d = DielectricTensorGenerator.from_objects(PhononBands.from_file(phbstnc_fname),
                                                   AnaddbNcFile.from_file(anaddbnc_fname))

        self.assertAlmostEqual(d.tensor_at_frequency(0.001, units='Ha')[0,0], 11.917178540635028)

        if self.has_matplotlib():
            assert d.plot_vs_w(w_min=0.0001, w_max=0.01, num=10, units="Ha", show=False)
            assert d.plot_vs_w(w_min=0, w_max=None, num=10, units="cm-1", show=False)
            for comp in ["diag", "all", "diag_av"]:
                assert d.plot_vs_w(num=10, component=comp, units="cm-1", show=False)


class DdbRobotTest(AbipyTest):

    def test_ddb_robot(self):
        """Testing DDB robots."""
        assert not abilab.DdbRobot.class_handles_filename("foo_DDB.nc")
        assert abilab.DdbRobot.class_handles_filename("foo_DDB")

        path = abidata.ref_file("refs/znse_phonons/ZnSe_hex_qpt_DDB")
        robot = abilab.DdbRobot.from_files(path)
        robot.add_file("same_ddb", path)
        repr(robot); str(robot)
        assert robot.to_string(verbose=2)
        assert len(robot) == 2
        assert robot.EXT == "DDB"

        data = robot.get_dataframe_at_qpoint(qpoint=[0, 0, 0], asr=2, chneut=1,
                dipdip=0, with_geo=True, abspath=True)
        assert "mode1" in data and "alpha" in data

        r = robot.anaget_phonon_plotters(nqsmall=2, ndivsm=2, dipdip=0, verbose=2)
        if self.has_matplotlib():
            assert r.phbands_plotter.gridplot(show=False)
            assert r.phdos_plotter.gridplot(show=False)

        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()


class PhononComputationTest(AbipyTest):

    def test_phonon_computation(self):
        """Testing if pjdoses compute by anaddb integrate to 3*natom"""
        path = os.path.join(abidata.dirpath, "refs", "mgb2_phonons_nkpt_tsmear", "mgb2_121212k_0.04tsmear_DDB")
        ddb = abilab.abiopen(path)

        #for dos_method in ("gaussian",):
        #for dos_method in ("tetra",):
        for dos_method in ("tetra", "gaussian"):
            # Get phonon bands and Dos with anaddb.
            phbands_file, phdos_file = ddb.anaget_phbst_and_phdos_files(nqsmall=4, ndivsm=2,
                dipdip=0, chneut=0, dos_method=dos_method, lo_to_splitting=False, verbose=1)

            phbands, phdos = phbands_file.phbands, phdos_file.phdos
            natom3 = len(phbands.structure) * 3

            # Test that amu is present with correct values.
            assert phbands.amu is not None
            self.assert_almost_equal(phbands.amu[12.0], 0.24305e+02)
            self.assert_almost_equal(phbands.amu[5.0], 0.10811e+02)
            self.assert_almost_equal(phbands.amu_symbol["Mg"], phbands.amu[12.0])
            self.assert_almost_equal(phbands.amu_symbol["B"],  phbands.amu[5.0])

            # Total PHDOS should integrate to 3 * natom
            # Note that anaddb does not renormalize the DOS so we have to increate the tolerance.
            #E       Arrays are not almost equal to 2 decimals
            #E        ACTUAL: 8.9825274146312282
            #E        DESIRED: 9
            self.assert_almost_equal(phdos.integral_value, natom3, decimal=1)

            # Test convertion to eigenvectors. Verify that they are orthonormal
            cidentity = np.eye(natom3, dtype=np.complex)
            eig = phbands.dyn_mat_eigenvect
            for iq in range(phbands.nqpt):
                #print("About to test iq", iq, np.dot(eig[iq].T.conjugate(), eig[iq]))
                #assert np.allclose(np.dot(eig[iq], eig[iq].T), cidentity , atol=1e-5, rtol=1e-3)
                self.assert_almost_equal(np.dot(eig[iq].conjugate().T, eig[iq]), cidentity) #, decimal=1)

            # Summing projected DOSes over types should give the total DOS.
            pj_sum = sum(pjdos.integral_value for pjdos in phdos_file.pjdos_symbol.values())
            self.assert_almost_equal(phdos.integral_value, pj_sum)

            # Summing projected DOSes over types and directions should give the total DOS.
            # phdos_rc_type[ntypat, 3, nomega]
            values = phdos_file.reader.read_value("pjdos_rc_type").sum(axis=(0, 1))
            tot_dos = abilab.Function1D(phdos.mesh, values)
            self.assert_almost_equal(phdos.integral_value, tot_dos.integral_value)

            phbands_file.close()
            phdos_file.close()

        ddb.close()
