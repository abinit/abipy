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

        #filepath = "~/git_repos/abipy/lifetimes/alphaSiO2/run.abo_PHDOS.nc"
        filepath = os.path.join(abidata.dirpath, "refs", "mp-7000_DDB.bz2")
        with abilab.abiopen(filepath) as ddb:
            phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(nqsmall=2, ndivsm=1, mpi_procs=2)
            msqd_dos = phdos_file.msqd_dos
            # Read indsym from file and convert from F to C
            indsym = phdos_file.reader.read_value("indsym")
            indsym[:, :, 3] -= 1

            phbst_file.close()
            phdos_file.close()

        print(msqd_dos)
        #for fmt in ("cartesian", "cif", "ustar", "beta", "B"):
        for fmt in ("cartesian", "cif"):
            #print("Format:", fmt)
            df = msqd_dos.get_dataframe(temp=100, view="all", select_symbols="Si", fmt=fmt)
            abilab.print_dataframe(df)

        # Compare Abinit indsym with AbiPy version.
        abispg = msqd_dos.structure.abi_spacegroup
        abipy_indsym = indsym_from_symrel(abispg.symrel, abispg.tnons, msqd_dos.structure, tolsym=1e-8)
        assert np.all(abipy_indsym == indsym)
        assert np.all(abipy_indsym == msqd_dos.structure.indsym)

        # Equivalent atoms should have same determinant.
        #self.assert_almost_equal(df["determinant"].values, df["determinant"].values[0])

        #cif_string = msqd_dos.get_cif_string(self, temp=300)
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
            assert msqd_dos.plot_tensor(show=False)
            assert msqd_dos.plot_uiso(show=False)
            #assert msqd_dos.plot(view="all")
            #assert msqd_dos.plot_tensor(view="all")
            #assert msqd_dos.plot_uiso(view="all")
            #assert msqd_dos.plot_uiso(view="all", what="vel")
