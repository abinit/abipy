"""Tests for msqdos module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab
from abipy.core.symmetries import indsym_from_symrel


class MsqdTest(AbipyTest):

    def test_from_ddb(self):
        """Testing MsqDos from DDB file."""
        self.skip_if_abinit_not_ge("8.10.2")

        filepath = os.path.join(abidata.dirpath, "refs", "mp-7000_DDB.bz2")
        with abilab.abiopen(filepath) as ddb:
            phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(nqsmall=2, ndivsm=1, mpi_procs=2)
            msqd_dos = phdos_file.msqd_dos
            # Read indsym from file and convert from F to C
            indsym = phdos_file.reader.read_value("indsym")
            indsym[:, :, 3] -= 1
            phbst_file.close()
            phdos_file.close()

        repr(msqd_dos); str(msqd_dos)
        print(msqd_dos.to_string(verbose=2))
        for fmt in ("cartesian", "cif", "ustar", "beta"): #, "B"):
            df = msqd_dos.get_dataframe(temp=100, view="all", select_symbols="Si", fmt=fmt)
            abilab.print_dataframe(df, title="Format: %s" % fmt)

        # Equivalent atoms should have same determinant.
        df = msqd_dos.get_dataframe(temp=300, view="all", fmt="cartesian")
        for _, group in df.groupby(by="element"):
            self.assert_almost_equal(group["determinant"].values, group["determinant"].values[0])

        # Compare Abinit indsym with AbiPy version.
        abispg = msqd_dos.structure.abi_spacegroup
        abipy_indsym = indsym_from_symrel(abispg.symrel, abispg.tnons, msqd_dos.structure, tolsym=1e-8)
        assert np.all(abipy_indsym == indsym)
        assert np.all(abipy_indsym == msqd_dos.structure.indsym)

        cif_string = msqd_dos.get_cif_string(temp=300)
        #print(cif_string)
        self.assertMultiLineEqual(cif_string, """\
# generated using pymatgen
data_SiO2
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   4.94990566
_cell_length_b   4.94990566
_cell_length_c   5.44089660
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   120.00000000
_symmetry_Int_Tables_number   1
_chemical_formula_structural   SiO2
_chemical_formula_sum   'Si3 O6'
_cell_volume   115.45026866
_cell_formula_units_Z   3
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Si  Si1  1  0.528855  0.000000  0.833333  1
  Si  Si2  1  0.471145  0.471145  0.500000  1
  Si  Si3  1  0.000000  0.528855  0.166667  1
  O  O4  1  0.413167  0.147706  0.620242  1
  O  O5  1  0.852294  0.265462  0.953576  1
  O  O6  1  0.734538  0.586833  0.286909  1
  O  O7  1  0.265462  0.852294  0.046424  1
  O  O8  1  0.147706  0.413167  0.379758  1
  O  O9  1  0.586833  0.734538  0.713091  1
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_23
_atom_site_aniso_U_13
_atom_site_aniso_U_12
Si1    0.00850    0.00695    0.00611   -0.00019   -0.00009    0.00348
Si2    0.00850    0.00850    0.00611    0.00009   -0.00009    0.00502
Si3    0.00695    0.00850    0.00611    0.00009    0.00019    0.00348
O4    0.01916    0.01120    0.01353    0.00249   -0.00411    0.00762
O5    0.01120    0.01512    0.01353   -0.00660   -0.00249    0.00358
O6    0.01512    0.01916    0.01353    0.00411    0.00660    0.01153
O7    0.01512    0.01120    0.01353    0.00249    0.00660    0.00358
O8    0.01120    0.01916    0.01353    0.00411   -0.00249    0.00762
O9    0.01916    0.01512    0.01353   -0.00660   -0.00411    0.01153""")

        filepath = msqd_dos.write_cif_file(filepath=None, temp=300)
        assert filepath.endswith(".cif")
        same_structure = abilab.abiopen(filepath)
        assert same_structure.formula == msqd_dos.structure.formula
        assert len(same_structure) == len(msqd_dos.structure)
        # NB: lattice.matrix and cart_coords are not necessarily the
        # same when we read the structure from CIF because the lattice
        # is initialized from_angles_and_lenghts
        #self.assert_almost_equal(same_structure.lattice.matrix, msqd_dos.structure.lattice.matrix)
        for s1, s2 in zip(same_structure, msqd_dos.structure):
            assert s1.specie.symbol == s2.specie.symbol
            self.assert_almost_equal(s1.frac_coords, s2.frac_coords, decimal=5)
            #self.assert_almost_equal(s1.coords, s2.coords, decimal=5)

        maxerr = msqd_dos.check_site_symmetries(temp=300, verbose=1)
        assert maxerr < 1e-10

        if self.has_matplotlib():
            assert msqd_dos.plot(show=False)
            assert msqd_dos.plot(view="all", show=False)
            assert msqd_dos.plot_tensor(show=False)
            assert msqd_dos.plot_tensor(view="all", show=False)
            assert msqd_dos.plot_uiso(show=False)
            assert msqd_dos.plot_uiso(view="all", show=False)
            assert msqd_dos.plot_uiso(view="all", what="vel", show=False)
