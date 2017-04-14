"""Test for GSR module"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata
import abipy.core

from pprint import pprint
from abipy.core.testing import AbipyTest
from abipy.electrons.gsr import GsrReader, GsrFile


class GSRReaderTestCase(AbipyTest):

    def test_read_Si2(self):
        """Test the reading of GSR file."""
        path = abidata.ref_file("si_scf_GSR.nc")

        ref_dims = {
            "number_of_spins": 1
        }

        ref_int_values = {
            "space_group": 227,
        }

        ref_float_values = {
            "etotal": -8.8652767680604807,
        #    "primitive_vectors": np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125,
        #                                     5.125, 5.125, 0], (3,3)),
        }

        with GsrReader(path) as r:
            assert r.ngroups == 1
            varnames = r.read_varnames()
            assert varnames and "ecut" in varnames

            # Test dimensions.
            for dimname, int_ref in ref_dims.items():
                value = r.read_dimvalue(dimname)
                self.assert_equal(value, int_ref)

            # Test int variables
            for varname, int_ref in ref_int_values.items():
                value = r.read_value(varname)
                self.assert_equal(value, int_ref)

            # Test float variables
            for varname, float_ref in ref_float_values.items():
                value = r.read_value(varname)
                self.assert_almost_equal(value, float_ref)

            # Reading non-existent variables or dims should raise a subclass of NetcdReder.
            with self.assertRaises(GsrReader.Error): r.read_value("foobar")
            with self.assertRaises(GsrReader.Error): r.read_dimvalue("foobar")

            r.print_tree()
            for group in r.walk_tree():
                print("group: " + str(group))

            # Initialize pymatgen structure from GSR.
            structure = r.read_structure()
            assert isinstance(structure, abipy.core.Structure)


class GSRFileTestCase(AbipyTest):

    def test_gsr_silicon(self):
        """spin unpolarized GSR file"""

        with GsrFile(abidata.ref_file("si_scf_GSR.nc")) as gsr:
            assert gsr.basename == "si_scf_GSR.nc"
            assert gsr.relpath == os.path.relpath(abidata.ref_file("si_scf_GSR.nc"))
            assert gsr.filetype
            assert gsr.filestat()
            assert len(gsr.ncdump())
            repr(gsr)
            str(gsr)
            str(gsr.ebands)
            assert gsr.filepath == abidata.ref_file("si_scf_GSR.nc")
            assert gsr.nsppol == 1
            assert gsr.mband == 8 and gsr.nband == 8 and gsr.nelect == 8 and len(gsr.kpoints) == 29
            self.assert_almost_equal(gsr.energy.to("Ha"), -8.86527676798556)
            self.assert_almost_equal(gsr.energy_per_atom * len(gsr.structure), gsr.energy)

            # Test energy_terms
            eterms = gsr.energy_terms
            repr(eterms)
            str(eterms)
            assert eterms.to_string(with_doc=True)
            self.assert_almost_equal(eterms.e_xc.to("Ha"), -3.51815936301812)
            self.assert_almost_equal(eterms.e_nonlocalpsp.to("Ha"), 1.91660690901782)
            self.assert_almost_equal(eterms.e_kinetic.to("Ha"), 2.96421325671218)
            self.assert_almost_equal(eterms.e_fermie.to("Ha"), 0.205739364929368)

            # Forces and stress
            self.assert_almost_equal(gsr.cart_forces.to("Ha bohr^-1").flat,
               [-1.14726679671674e-28, -3.76037290483622e-29, 5.65937773808884e-29,
                 1.14726679671674e-28, 3.76037290483622e-29, -5.65937773808884e-29])

            self.assert_almost_equal(gsr.max_force, 0)
            assert gsr.force_stats()
            assert gsr.residm > 0
            assert str(gsr.xc) == "LDA_XC_TETER93"

            #self.assert_almost_equal(gsr.cart_stress_tensor.flat,
            # Cartesian components of stress tensor (hartree/bohr^3)
            #  sigma(1 1)=  1.77139311E-04  sigma(3 2)=  0.00000000E+00
            #  sigma(2 2)=  1.77139311E-04  sigma(3 1)=  0.00000000E+00
            #  sigma(3 3)=  1.77139311E-04  sigma(2 1)=  2.67294316E-15
            self.assert_almost_equal(gsr.pressure, -5.21162150)

            # Test as_dict
            d = gsr.as_dict()
            assert d and "structure" in d
            #import json
            #with open("hello.json", "w") as fp:
            #    json.dump(gsr.as_dict(), fp)
            #assert 0

            # Test pymatgen computed_entries
            for inc_structure in (True, False):
                e = gsr.get_computed_entry(inc_structure=inc_structure)
                str(e)
                d = e.as_dict()
                if inc_structure: assert "structure" in d
                assert d["energy"] == gsr.energy
                assert gsr.energy == e.energy

            if self.has_nbformat():
                gsr.write_notebook(nbpath=self.get_tmpname(text=True))
