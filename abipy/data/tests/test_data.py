"""Tests for data module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata

from abipy.core.testing import AbipyTest


class DataModuleTest(AbipyTest):
    """Test data module."""

    def test_data_api(self):
        """Testing abipy.data API."""
        ncpseudo_h = abidata.pseudo("01h.pspgth")
        assert ncpseudo_h.isnc
        paw_table = abidata.pseudos("1h.paw", "28ni.paw")
        assert paw_table.allpaw
        assert os.path.isfile(abidata.cif_file("al.cif"))
        assert os.path.isfile(abidata.pyscript("plot_bz.py"))

        d = abidata.get_mp_structures_dict()
        assert isinstance(d, dict) and d is abidata.get_mp_structures_dict()

        structure = abidata.structure_from_cif("gan2.cif")
        assert hasattr(structure, "to_abivars")

        structure = abidata.structure_from_mpid("mp-4820")
        assert hasattr(structure, "to_abivars")
        with self.assertRaises(KeyError):
            abidata.structure_from_mpid("foobar")


class FilesGeneratorTest(AbipyTest):

    def test_abinit_files_generator(self):
        """Testing AbinitFilesGenerator."""
        class MyGenerator(abidata.AbinitFilesGenerator):
            """This class generates the output files used in the unit tests and in the examples."""
            # Subclasses must define the following class attributes:
            # List of pseudos (basenames) in abipy/data/pseudos
            pseudos = ["14si.pspnc"]

            # Mapping old_name --> new_name for the output files that must be preserved.
            files_to_save = {
                "out_DS1_DEN.nc": "si_DEN.nc",
                "out_DS2_GSR.nc": "si_nscf_GSR.nc",
                "out_DS2_WFK.nc": "si_nscf_WFK.nc",
                "out_DS1_GSR.nc": "si_scf_GSR.nc",
                "out_DS1_WFK.nc": "si_scf_WFK.nc",
            }

        abgen = MyGenerator()
        print(abgen)
        assert abgen.executable == "abinit"
        assert len(abgen.make_filesfile_str())

    def test_anaddb_files_generator(self):
        """Testing AnaddbFilesGenerator."""
        class MyGenerator(abidata.AnaddbFilesGenerator):
            """This class generates the output files used in the unit tests and in the examples."""

            # Mapping old_name --> new_name for the output files that must be preserved.
            files_to_save = {
                "out_PHBST.nc": "trf2_5.out_PHBST.nc",
                "out_PHDOS.nc": "trf2_5.out_PHDOS.nc",
            }

            in_ddb = "trf2_3.ddb.out"

        anagen = MyGenerator()
        print(anagen)
        assert anagen.executable == "anaddb"
        assert len(anagen.make_filesfile_str())
