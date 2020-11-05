"""Tests for dfpt converters"""
import os
import tempfile
import numpy as np
import abipy.core.abinit_units as abu

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.dfpt.ddb import DdbFile

from abipy.dfpt.converters import abinit_to_phonopy, phonopy_to_abinit, tdep_to_abinit
from abipy.dfpt.converters import born_to_lotosplitting

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.phonopy import get_phonopy_structure

try:
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS, parse_BORN
except ImportError:
    Phonopy = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


def find_anaddbnc_in_dir(dirpath):
    for n in ("anaddb.nc", "run_anaddb.nc"):
        p = os.path.join(dirpath, n)
        if os.path.isfile(p): return p

    raise RuntimeError(f"Cannot file anaddb.nc file in {dirpath}")


class ConverterTest(AbipyTest):

    def test_ddb_phonopy_ddb_phfreqs(self):
        """
        Checks the phonon frequencies when converting starting from a DDB then to phonopy
        and back to DDB.
        """
        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()
        ngqpt = [4, 4, 4]
        with DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB")) as ddb:
            qpoints = ddb.qpoints
            structure = ddb.structure
            phonon = ddb.anaget_phonopy_ifc(ngqpt=ngqpt, asr=0, dipdip=0, chneut=0,
                                            output_dir_path=tmp_dir, set_masses=True)

            orig_phbands = ddb.anaget_phmodes_at_qpoints(qpoints=qpoints, asr=0, dipdip=0, chneut=0,
                                                         lo_to_splitting=False)

        # remove nac_params to be sure that this does not enter into play.
        # There should be no nac_params anyway here.
        phonon.nac_params = None
        phfreqs_phonopy = np.array([phonon.get_frequencies(q.frac_coords) for q in qpoints])

        out_ddb_path = os.path.join(tmp_dir, "out_DDB")

        scm = np.eye(3) * ngqpt
        ddb_conv = phonopy_to_abinit(unit_cell=structure, supercell_matrix=scm, qpt_list=qpoints.frac_coords,
                                     out_ddb_path=out_ddb_path, force_constants=os.path.join(tmp_dir, "FORCE_CONSTANTS"),
                                     born=None, primitive_matrix=np.eye(3), symprec=1e-5, tolsym=None)

        conv_phbands = ddb_conv.anaget_phmodes_at_qpoints(qpoints=qpoints, asr=0, dipdip=0, chneut=0,
                                                          lo_to_splitting=False)

        self.assertArrayAlmostEqual(orig_phbands.phfreqs, conv_phbands.phfreqs, decimal=5)
        self.assertArrayAlmostEqual(orig_phbands.phfreqs * abu.eV_to_THz, phfreqs_phonopy, decimal=3)

    def test_ddb_phonopy_ddb_becs(self):
        """
        Checks the phonon frequencies, BECs and dielectric tensor when converting starting from
        a DDB then to phonopy and back to DDB.
        """

        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()
        orig_run_ana = os.path.join(tmp_dir, "anaddb_orig")
        conv_run_ana = os.path.join(tmp_dir, "anaddb_conv")
        ngqpt = [1, 1, 1]
        with DdbFile(os.path.join(test_dir, "ZnO_gamma_becs_DDB")) as ddb:
            qpoints = ddb.qpoints
            structure = ddb.structure
            phonon = ddb.anaget_phonopy_ifc(ngqpt=ngqpt, asr=0, dipdip=0, chneut=0,
                                            output_dir_path=tmp_dir, set_masses=True)

            orig_phbands = ddb.anaget_phmodes_at_qpoints(qpoints=qpoints, asr=0, dipdip=0, chneut=0,
                                                         lo_to_splitting=False, workdir=orig_run_ana)
            ananc_orig = abilab.abiopen(find_anaddbnc_in_dir(orig_run_ana))

        nac = phonon.nac_params
        # remove nac_params to be sure that this does not enter into play.
        phonon.nac_params = None
        phfreqs_phonopy = np.array([phonon.get_frequencies(q.frac_coords) for q in qpoints])

        out_ddb_path = os.path.join(tmp_dir, "out_DDB")

        scm = np.eye(3) * ngqpt
        ddb_conv = phonopy_to_abinit(unit_cell=structure, supercell_matrix=scm, qpt_list=qpoints.frac_coords,
                                     out_ddb_path=out_ddb_path, force_constants=os.path.join(tmp_dir, "FORCE_CONSTANTS"),
                                     born=os.path.join(tmp_dir, "BORN"), primitive_matrix=np.eye(3), symprec=1e-5, tolsym=None)

        conv_phbands = ddb_conv.anaget_phmodes_at_qpoints(qpoints=qpoints, asr=0, dipdip=0, chneut=0,
                                                          lo_to_splitting=False, workdir=conv_run_ana)
        ananc_conv = abilab.abiopen(find_anaddbnc_in_dir(conv_run_ana))

        self.assertArrayAlmostEqual(orig_phbands.phfreqs, conv_phbands.phfreqs, decimal=5)
        self.assertArrayAlmostEqual(orig_phbands.phfreqs * abu.eV_to_THz, phfreqs_phonopy, decimal=3)
        self.assertArrayAlmostEqual(ananc_orig.becs.values, nac["born"], decimal=5)
        self.assertArrayAlmostEqual(ananc_conv.becs.values, nac["born"], decimal=5)
        self.assertArrayAlmostEqual(ananc_orig.epsinf, nac["dielectric"], decimal=5)
        self.assertArrayAlmostEqual(ananc_conv.epsinf, nac["dielectric"], decimal=5)

    def test_phonopy_ddb_phonopy(self):

        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()

        ngqpt = [2, 2, 2]
        scm = np.eye(3) * ngqpt

        tmp_dir = tempfile.mkdtemp()
        out_ddb_path = os.path.join(tmp_dir, "out_DDB")

        unit_cell = abilab.Structure.from_file(os.path.join(test_dir, "MgO_phonopy_POSCAR"))
        fc_path = os.path.join(test_dir, "MgO_phonopy_FORCE_CONSTANTS")
        born_path = os.path.join(test_dir, "MgO_phonopy_BORN")

        ddb = phonopy_to_abinit(unit_cell=unit_cell, supercell_matrix=scm, ngqpt=ngqpt,
                                out_ddb_path=out_ddb_path, force_constants=fc_path,
                                born=born_path, primitive_matrix=np.eye(3), symprec=1e-5,
                                tolsym=None)

        qpoints = ddb.qpoints
        phonon_conv = ddb.anaget_phonopy_ifc(ngqpt=ngqpt, asr=0, dipdip=0, chneut=0,
                                             output_dir_path=tmp_dir, set_masses=True)

        conv_nac = phonon_conv.nac_params
        phonon_conv.nac_params = None

        orig_fc = parse_FORCE_CONSTANTS(fc_path)
        orig_nac = parse_BORN(phonon_conv.primitive, filename=born_path)
        # create the Phonopy object here to take advantage of the information in phonopy_conv
        phonon_orig = Phonopy(unitcell=phonon_conv.unitcell, supercell_matrix=scm, primitive_matrix=np.eye(3),
                              nac_params=None)
        phonon_orig.set_force_constants(orig_fc)

        phfreqs_phonopy_orig = np.array([phonon_orig.get_frequencies(q.frac_coords) for q in qpoints])
        phfreqs_phonopy_conv = np.array([phonon_conv.get_frequencies(q.frac_coords) for q in qpoints])

        run_ana = os.path.join(tmp_dir, "run_anaddb")
        phbands = ddb.anaget_phmodes_at_qpoints(qpoints=qpoints, asr=0, dipdip=0, chneut=0,
                                                lo_to_splitting=False, workdir=run_ana)
        ananc = abilab.abiopen(find_anaddbnc_in_dir(run_ana))

        self.assertArrayAlmostEqual(phbands.phfreqs * abu.eV_to_THz, phfreqs_phonopy_orig, decimal=3)
        self.assertArrayAlmostEqual(phbands.phfreqs * abu.eV_to_THz, phfreqs_phonopy_conv, decimal=3)
        self.assertArrayAlmostEqual(ananc.becs.values, orig_nac["born"], decimal=5)
        self.assertArrayAlmostEqual(ananc.becs.values, conv_nac["born"], decimal=5)
        self.assertArrayAlmostEqual(ananc.epsinf, orig_nac["dielectric"], decimal=5)
        self.assertArrayAlmostEqual(ananc.epsinf, conv_nac["dielectric"], decimal=5)

    def test_tdep_lotosplitting(self):
        tmp_dir = tempfile.mkdtemp()
        primitive = get_phonopy_structure(abilab.Structure.from_file(os.path.join(test_dir, "MgO_phonopy_POSCAR")))
        born_path = parse_BORN(primitive, filename=os.path.join(test_dir, "MgO_phonopy_BORN"))
        loto_path = os.path.join(tmp_dir, "infile.lotosplitting")
        born_to_lotosplitting(born_path, loto_path)
        self.assertTrue(os.path.isfile(loto_path))

    def test_tdep_ddb(self):
        tmp_dir = tempfile.mkdtemp()
        out_ddb_path = os.path.join(tmp_dir, "out_DDB")

        # need to pass through Poscar, since Structure.from_file fails due to tdep file names
        unit_cell = Poscar.from_file(os.path.join(test_dir, "Si_tdep.ucposcar")).structure
        supercell = Poscar.from_file(os.path.join(test_dir, "Si_tdep.ssposcar")).structure
        fc_path = os.path.join(test_dir, "Si_tdep.forceconstant")

        ddb = tdep_to_abinit(unit_cell=unit_cell, fc_path=fc_path, supercell_matrix=np.eye(3)*5,
                             supercell=supercell, out_ddb_path=out_ddb_path, ngqpt=[5, 5, 5],
                             primitive_matrix=np.eye(3), lotosplitting_path=None)

        phbands = ddb.anaget_phmodes_at_qpoints(qpoints=[[0, 0, 0]], asr=0, dipdip=0, chneut=0,
                                                lo_to_splitting=False)
        self.assertAlmostEqual(phbands.phfreqs[0, 3], 0.062997, places=3)
