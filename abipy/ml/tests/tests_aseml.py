"""Tests for aseml module"""
import numpy as np
import abipy.data as abidata
import abipy.ml.aseml as aseml

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure


class AbimlTest(AbipyTest):

    def test_mlrelaxer(self):
        """Testing MlRelaxer."""
        ml_relaxer = aseml.MlRelaxer(atoms=Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms(),
                                     relax_mode="cell",
                                     fmax=0.01,
                                     pressure=0.0,
                                     steps=100,
                                     optimizer="BFGS",
                                     nn_name="emt",
                                     verbose=1,
                                     workdir=None,
                                     prefix=None
                                     )

        assert ml_relaxer.to_string(verbose=1)
        print(ml_relaxer)
        ml_relaxer.run()
        ml_relaxer.rmtree()

    def test_mlmd(self):
        """Testing MlMd."""
        md = aseml.MlMd(atoms=Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms(),
                        temperature=300,
                        pressure=0,
                        timestep=1,
                        steps=30,
                        loginterval=10,
                        ensemble="nvt",
                        nn_name="emt",
                        verbose=1,
                        workdir=None,
                        prefix=None
                 )

        assert md.to_string(verbose=1)
        print(md)
        md.run()
        md.rmtree()

    def test_mlneb(self):
        """Testing MlNeb."""
        initial_atoms = Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms()
        final_atoms = initial_atoms.copy()
        final_atoms.positions += 0.1

        neb = aseml.MlNeb(initial_atoms,
                          final_atoms,
                          nimages=3,
                          neb_method="aseneb",
                          climb=False,
                          optimizer="BFGS",
                          relax_mode="ions",
                          fmax=0.01,
                          pressure=0,
                          nn_name="emt",
                          verbose=1,
                          workdir=None,
                          prefix=None,
                          )


        assert neb.to_string(verbose=1)
        print(neb)
        neb.run()
        neb.rmtree()

    def test_multi_neb(self):
        """Testing MultiMlNeb."""
        initial_atoms = Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms()
        second_atoms = initial_atoms.copy()
        second_atoms.positions += 0.1
        final_atoms = initial_atoms.copy()
        final_atoms.positions += 0.2

        atoms_list = [initial_atoms, second_atoms, final_atoms]

        multi_neb = aseml.MultiMlNeb(atoms_list=atoms_list,
                                     nimages=3,
                                     neb_method="improvedtangent",
                                     climb=False,
                                     optimizer="FIRE",
                                     relax_mode="cell",
                                     fmax=0.01,
                                     pressure=0,
                                     nn_name="emt",
                                     verbose=1,
                                     workdir=None,
                                     prefix=None
                                     )

        assert multi_neb.to_string(verbose=1)
        print(multi_neb)
        #multi_neb.run()
        #multi_neb.rmtree()

    # TODO
    #def test_mlordered(self):
    #    """Testing MlOrdered."""

    #def test_mlvalidate(self):
    #    """Testing MlValidateWithAbinitio."""
    #    filepaths = abidata.ref_file("sic_relax_HIST.nc")

    #    validator = aseml.MlValidateWithAbinitio(filepaths,
    #                                             nn_names="emt",
    #                                             traj_range=range(0, 5, 2),
    #                                             verbose=1,
    #                                             workdir=None,
    #                                             prefix=None
    #                                             )

    #    assert validator.to_string(verbose=1)
    #    print(validator)
    #    validator.run(nprocs=1, print_dataframes=True)
    #    validator.rmtree()

    def test_mleos(self):
        """Testing MlEos."""
        atoms = Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms()

        ml_eos = aseml.MlEos(atoms=atoms,
                             vol_scales=np.arange(0.95, 1.06, 0.01),
                             relax_mode="cell",
                             fmax=0.01,
                             pressure=0,
                             steps=100,
                             optimizer="FIRE",
                             nn_name="emt",
                             verbose=1,
                             workdir=None,
                             prefix=None)

        assert ml_eos.to_string(verbose=1)
        print(ml_eos)
        ml_eos.run()
        ml_eos.rmtree()

    def test_mlcompare_nns(self):
        """Testing MlEos."""
        atoms = Structure.as_structure(abidata.cif_file("al.cif")).to_ase_atoms()

        ml_compare = aseml.MlCompareNNs(atoms=atoms,
                                        nn_names=["emt", "emt"],
                                        num_tests=1,
                                        rattle=0.01,
                                        stdev_rvol=0.01,
                                        verbose=1,
                                        workdir=None,
                                        prefix=None
                                        )

        assert ml_compare.to_string(verbose=1)
        print(ml_compare)
        ml_compare.run()
        ml_compare.rmtree()
