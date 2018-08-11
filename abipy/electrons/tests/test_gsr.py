"""Test for GSR module"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata
import abipy.core
import abipy.core.abinit_units as abu

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
            repr(gsr); str(gsr)
            assert gsr.to_string(verbose=2)
            assert gsr.abinit_version == "8.0.6"
            str(gsr.ebands)
            assert gsr.filepath == abidata.ref_file("si_scf_GSR.nc")
            assert gsr.nsppol == 1
            assert gsr.mband == 8 and gsr.nband == 8 and gsr.nelect == 8 and len(gsr.kpoints) == 29
            assert gsr.mband == gsr.hdr.mband
            assert "nelect" in gsr.hdr and gsr.nelect == gsr.hdr.nelect
            self.assert_almost_equal(gsr.energy.to("Ha"), -8.86527676798556)
            self.assert_almost_equal(gsr.energy_per_atom * len(gsr.structure), gsr.energy)

            assert gsr.params["nband"] == 8
            assert gsr.params["nkpt"] == 29

            # Test energy_terms
            eterms = gsr.energy_terms
            repr(eterms); str(eterms)
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
            for i in range(3):
                self.assert_almost_equal(gsr.cart_stress_tensor[0, 0], 1.77139311E-04 * abu.HaBohr3_GPa)
            self.assert_almost_equal(gsr.pressure, -5.211617575719521)

            # Test pymatgen computed_entries
            for inc_structure in (True, False):
                e = gsr.get_computed_entry(inc_structure=inc_structure)
                str(e)
                d = e.as_dict()
                if inc_structure: assert "structure" in d
                assert d["energy"] == gsr.energy
                assert gsr.energy == e.energy

            if self.has_matplotlib():
                assert gsr.plot_ebands(show=False)
                assert gsr.plot_ebands_with_edos(edos=gsr.get_edos(), show=False)

            if self.has_nbformat():
                gsr.write_notebook(nbpath=self.get_tmpname(text=True))


class GsrRobotTest(AbipyTest):

    def test_gsr_robot(self):
        """Testing GSR robot"""
        from abipy import abilab
        gsr_ibz = abidata.ref_file("si_scf_GSR.nc")
        robot = abilab.GsrRobot()
        robot.add_file("gsr0", gsr_ibz)
        assert len(robot.abifiles) == 1
        assert "gsr0" in robot.keys()
        assert "gsr0" in robot.labels
        assert robot.EXT == "GSR"
        repr(robot); str(robot)
        assert robot.to_string(verbose=2)

	# Cannot have same label
        with self.assertRaises(ValueError):
            robot.add_file("gsr0", gsr_ibz)

        assert len(robot) == 1 and not robot.exceptions
        robot.add_file("gsr1", abilab.abiopen(gsr_ibz))
        assert len(robot) == 2
        robot.show_files()
        assert not robot.has_different_structures()
        with self.assertRaises(AttributeError):
            robot.is_sortable("foobar", raise_exc=True)
        assert not robot.is_sortable("foobar")
        # Test different syntax.
        assert robot.is_sortable("nkpt")         # gsr.nkpt
        assert robot.is_sortable("ebands.nkpt")  # gsr.ebands.nkpt
        assert robot.is_sortable("ecut")         # in gsr.params

        dfs = robot.get_structure_dataframes()
        assert dfs.lattice is not None
        assert dfs.coords is not None
        assert len(dfs.structures) == len(robot)

        assert len(robot.get_ebands_plotter(kselect="path")) == 0
        filter_abifile = lambda gsr: gsr.ebands.kpoints.is_ibz
        assert len(robot.get_ebands_plotter(filter_abifile=filter_abifile)) == 2

        ebands_plotter = robot.get_ebands_plotter()
        edos_plotter = robot.get_edos_plotter()

        if self.has_matplotlib():
            assert ebands_plotter.gridplot(show=False)
            assert robot.combiplot_ebands(show=False)
            assert robot.gridplot_ebands(show=False)
            assert robot.boxplot_ebands(show=False)
            assert robot.combiboxplot_ebands(show=False)

            assert edos_plotter.gridplot(show=False)
            assert robot.combiplot_edos(show=False)
            assert robot.gridplot_edos(show=False)

            assert robot.plot_gsr_convergence(show=False)
            assert robot.plot_gsr_convergence(sortby="nkpt", hue="tsmear", show=False)
            y_vars = ["energy", "structure.lattice.a", "structure.volume"]
            assert robot.plot_convergence_items(y_vars, sortby="nkpt", hue="tsmear", show=False)

            assert robot.plot_egaps(show=False)
            assert robot.plot_egaps(sortby="nkpt", hue="tsmear")
            assert robot.gridplot_with_hue("structure.formula", show=False)

	# Get pandas dataframe.
        df = robot.get_dataframe()
        assert "energy" in df
        self.assert_equal(df["ecut"].values, 6.0)
        self.assert_almost_equal(df["energy"].values, -241.2364683)

        df_params = robot.get_params_dataframe()
        assert "nband" in df_params

        assert "alpha" in robot.get_lattice_dataframe()
        assert hasattr(robot.get_coords_dataframe(), "keys")

        eterms_df = robot.get_energyterms_dataframe(iref=0)
        assert "energy" in eterms_df
        assert "ecut" in eterms_df
        assert "nkpt" in eterms_df

        if self.has_matplotlib():
            assert robot.plot_xy_with_hue(df, x="nkpt", y="pressure", hue="a", show=False)

        # Note: This is not a real EOS since we have a single volume.
        # But testing is better than not testing.
        r = robot.get_eos_fits_dataframe()
        assert hasattr(r, "fits") and hasattr(r, "dataframe")

        if self.has_matplotlib():
            assert robot.gridplot_eos(show=False)

        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()

        # Test other class methods
        filepath = abidata.ref_file("si_scf_GSR.nc")
        robot = abilab.GsrRobot.from_dirs(os.path.dirname(filepath), abspath=True)
        assert len(robot) == 2
        assert robot[filepath].filepath == filepath

        # Test from_dir_glob
        pattern = "%s/*/si_ebands/" % abidata.dirpath
        same_robot = abilab.GsrRobot.from_dir_glob(pattern, abspath=True)
        assert len(same_robot) == 2
        assert set(robot.labels) == set(same_robot.labels)

        robot.close()
        same_robot.close()
