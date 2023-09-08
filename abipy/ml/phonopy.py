"""
Based on a similar implementation available at: https://github.com/modl-uclouvain/randomcarbon/blob/main/randomcarbon/run/phonon.py
"""
from __future__ import annotations

import numpy as np
import logging

#from monty.json import MSONable, MontyDecoder
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms
from pymatgen.io.phonopy import get_phonopy_structure
from monty.dev import requires
from abipy.core.structure import Structure
from abipy.dfpt.ddb import DdbFile
from abipy.ml.aseml import RX_MODE, CalcBuilder, AseResults, MlBase, get_atoms, relax_atoms, dataframe_from_results_list
try:
    import phonopy
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except ImportError:
    phonopy = None
    Phonopy = None


@requires(phonopy, "phonopy should be installed to calculate phonons")
def get_phonopy(structure: Structure,
                supercell_matrix,
                calculator: Calculator,
                distance=0.03,
                primitive_matrix=None) -> Phonopy:

    structure = Structure.as_structure(structure)
    unitcell = get_phonopy_structure(structure)
    #if supercell_matrix is None:
    #    supercell_matrix = np.eye(3)

    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, primitive_matrix=primitive_matrix)

    phonon.generate_displacements(distance=distance)

    forces = []
    nsc = len(phonon.supercells_with_displacements)
    print(f"Calculating forces for {nsc} supercell configurations ...")
    for i, sc in enumerate(phonon.supercells_with_displacements):
        a = Atoms(symbols=sc.symbols,
                  positions=sc.positions,
                  masses=sc.masses,
                  cell=sc.cell, 
                  pbc=True,
                  constraint=None,
                  calculator=calculator)
        #print(f"{i+1} of {nsc}")
        forces.append(a.get_forces())

    phonon.set_forces(forces)
    phonon.produce_force_constants()

    return phonon

    #from phonopy.units import Hartree, Bohr
    #from phonopy.structure.symmetry import symmetrize_borns_and_epsilon 
    #borns_, epsilon_ = symmetrize_borns_and_epsilon(
    #borns,
    #epsilon,
    #unitcell,
    #primitive_matrix=[[0, 0.5, 0.5],
    #                  [0.5, 0, 0.5],
    #                  [0.5, 0.5, 0]],
    #supercell_matrix=np.diag([2, 2, 2]),
    #symprec=1e-5)
    #nac_params = {'born': borns_,
    #              'factor': Hartree * Bohr,
    #              'dielectric': epsilon_}


class MlPhononsWithDDB(MlBase):
    """
    Compute phonons with phonopy and ML potential starting from a DDB file.
    and compare the results.
    """
    def __init__(self, ddb_filepath: str, supercell, distance, qmesh, asr, nqpath,
                 relax_mode, fmax, pressure, steps, optimizer, nn_name,
                 verbose, workdir, prefix=None):
        """
        Args:
            ddb_filepath: DDB filename.
            supercell: tuple with supercell dimension.
            distance:
            qmesh: q-mesh for phonon-DOS
            asr: Enforce acoustic sum-rule.
            nqpath: Number of q-point along the q-path.
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
            verbose: whether to print stdout.
            optimizer: name of the ASE optimizer to use for relaxation.
            nn_name: String defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)

        self.ddb_filepath = ddb_filepath
        with DdbFile(self.ddb_filepath) as ddb:
            self.initial_atoms = get_atoms(ddb.structure)
        self.supercell = supercell
        self.distance = float(distance)
        self.qmesh = qmesh
        self.asr = asr
        self.nqpath = nqpath
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.steps = steps
        self.optimizer = optimizer
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     ddb_path   = {self.ddb_filepath}
     supercell  = {self.supercell}
     distance   = {self.distance}
     qmesh      = {self.qmesh}
     asr        = {self.asr}
     nqpath     = {self.nqpath}
     relax_mode = {self.relax_mode}
     fmax       = {self.fmax}
     steps      = {self.steps}
     optimizer  = {self.optimizer}
     pressure   = {self.pressure}
     nn_name    = {self.nn_name}
     workdir    = {self.workdir}
     verbose    = {self.verbose}

=== ATOMS ===

{self.initial_atoms}
"""
        return s

    def run(self) -> None:
        """Run MlPhononsWithDDB."""
        workdir = self.workdir

        calculator = CalcBuilder(self.nn_name).get_calculator()
        atoms = self.initial_atoms.copy()
        atoms.calc = calculator

        if self.relax_mode != RX_MODE.no:
            print(f"Relaxing input DDB atoms with relax mode: {self.relax_mode}.")
            relax_kws = dict(optimizer=self.optimizer,
                             relax_mode=self.relax_mode,
                             fmax=self.fmax,
                             pressure=self.pressure,
                             steps=self.steps,
                             traj_path=self.get_path("relax.traj", "ASE relax trajectory"),
                             verbose=self.verbose,
                            )

            relax = relax_atoms(atoms, **relax_kws)
            r0, r1 = AseResults.from_traj_inds(relax.traj, 0, -1)
            df = dataframe_from_results_list(["DDB_initial", "DDB_relaxed"], [r0, r1])
            print(df, end=2*"\n")

        # Call phonopy to compute phonons with finite difference and ML potential.
        # Include non-analytical term if dipoles are available in the DDB file.
        phonon = get_phonopy(atoms, self.supercell, calculator, distance=self.distance, primitive_matrix=None)

        with DdbFile(self.ddb_filepath) as ddb:
            if ddb.has_lo_to_data:
                print("Activating dipolar term in phonopy calculation using BECS and Zeff taken from DDB.")
                out = ddb.anaget_epsinf_and_becs(chneut=1)
                becs = out.becs.values
                epsinf = out.epsinf
                # according to the phonopy website 14.399652 is not the coefficient for abinit
                # probably it relies on the other conventions in the output.
                phonon.nac_params = {"born": becs, "dielectric": epsinf, "factor": 14.399652}

            # Call anaddb to compute ab-initio phonons from the DDB.
            #ddb.anaget_phmodes_at_qpoints(
            #   qpoints=None, asr=2, chneut=1, dipdip=1, dipquad=1, quadquad=1,
            #   ifcflag=0, ngqpt=None, workdir=None, mpi_procs=1, manager=None, verbose=0,
            #   lo_to_splitting=False,
            #   spell_check=True, directions=None, anaddb_kwargs=None, return_input=False)


            #ddb.anaget_phbst_and_phdos_files(self, nqsmall=10, qppa=None, ndivsm=20, line_density=None, asr=2, chneut=1,
            #                                 dipdip=1, dipquad=1, quadquad=1,
            #                                 dos_method="tetra", lo_to_splitting="automatic", ngqpt=None, qptbounds=None,
            #                                 anaddb_kwargs=None, with_phonopy_obj=False, verbose=0, spell_check=True,
            #                                 mpi_procs=1, workdir=None, manager=None, return_input=False):



            # Band structure part.
            #qpath = [[[0, 0, 0], [0.5, 0, 0.5], [0.5, 0.25, 0.75], [0.375, 0.375, 0.75],
            #        [0, 0, 0], [0.5, 0.5, 0.5], [0.625, 0.25, 0.625],
            #        [0.5, 0.25, 0.75], [0.5, 0.5, 0.5], [0.375, 0.375, 0.75]],
            #        [[0.625, 0.25, 0.625], [0.5, 0, 0.5]]]
            #labels = ["$\\Gamma$", "X", "W", "K",
            #        "$\\Gamma$", "L", "U",
            #        "W", "L", "K",
            #        "U", "X"]

            labels = [q.name for q in ddb.structure.hsym_kpoints]
            qpath = [q.frac_coords.tolist() for q in ddb.structure.hsym_kpoints]
            print(f"{qpath=}\n{labels=}")

            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections, get_band_qpoints
            rec_lattice = ddb.structure.reciprocal_lattice.matrix.T
            qpoints, connections = get_band_qpoints_and_path_connections(qpath, npoints=51, rec_lattice=rec_lattice)
            qpoints, connections = get_band_qpoints(qpath, npoints=51, rec_lattice=rec_lattice)
            print(f"{qpoints=}\n{connections=}")
            #phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)

            #plt.figure(figsize=(20,20))
            #plt = phonon.plot_band_structure()
            #plt.show() 
            #plt.savefig(f"{folder_path}/{mp_id}_phonon_band_structure.png")
            #plt.close()
            #bands_dict = phonon.get_band_structure_dict()
            #qpoints = bands_dict['qpoints']
            #frequencies = bands_dict['frequencies']

            #phonon.save(settings={'force_constants': True})
