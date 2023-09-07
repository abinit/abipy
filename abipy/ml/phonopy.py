"""
Based on a similar implementation available at: https://github.com/modl-uclouvain/randomcarbon/blob/main/randomcarbon/run/phonon.py
"""
import numpy as np
import logging

from typing import Callable, Union, List
from pymatgen.io.ase import AseAtomsAdaptor
from monty.json import MSONable, MontyDecoder
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms
#from pymatgen.core import Structure, Site, Element, Composition
from pymatgen.io.phonopy import get_phonopy_structure
#from m3gnet.models import M3GNet, Potential, Relaxer, M3GNetCalculator
from phonopy import Phonopy
#from typing import Union, List
from monty.dev import requires
#from pymatgen.core.structure import Structure
#from pymatgen.io.phonopy import get_phonopy_structure
#from ase.atoms import Atoms
from abipy.dfpt.ddb import DdbFile
from abipy.ml.aseml import RX_MODE, CalcBuilder, AseResults
try:
    import phonopy
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except ImportError:
    phonopy = None
    Phonopy = None


#def get_phonons(structure: Structure, calculator: Union[Calculator, Factory], constraints: list = None,
@requires(phonopy, "phonopy should be installed to calculate phonons")
def get_phonopy(structure: Structure,
                calculator: Calculator,
                distance=0.03,
                constraints: list = None,
                supercell_matrix: List[List[int]] = None,
                primitive_matrix: List[List[float]] = None) -> Phonopy:

    unitcell = get_phonopy_structure(structure)
    if supercell_matrix is None:
        supercell_matrix = np.eye(3)

    phonon = Phonopy(unitcell,
                     supercell_matrix=supercell_matrix,
                     primitive_matrix=primitive_matrix)

    phonon.generate_displacements(distance=distance)
    supercells = phonon.supercells_with_displacements

    supercells_atoms = []
    for sc in supercells:
        a = Atoms(symbols=sc.symbols,
                  positions=sc.positions,
                  masses=sc.masses,
                  cell=sc.cell, pbc=True,
                  constraint=None,
                  calculator=calculator)
        if constraints:
            tmp_constraints = []
            for i, c in enumerate(constraints):
                if isinstance(c, Factory):
                    tmp_constraints.append(c.generate(atoms=a))
                else:
                    tmp_constraints.append(c)
            a.set_constraint(tmp_constraints)
        supercells_atoms.append(a)

    forces = []
    for i, sca in enumerate(supercells_atoms):
        #logger.debug(f"calculating forces for supercell {i+1} of {len(supercells_atoms)}")
        forces.append(sca.get_forces())

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

    # Band structure part.
    #path = [[[0, 0, 0], [0.5, 0, 0.5], [0.5, 0.25, 0.75], [0.375, 0.375, 0.75],
    #        [0, 0, 0], [0.5, 0.5, 0.5], [0.625, 0.25, 0.625],
    #        [0.5, 0.25, 0.75], [0.5, 0.5, 0.5], [0.375, 0.375, 0.75]],
    #        [[0.625, 0.25, 0.625], [0.5, 0, 0.5]]]
    #labels = ["$\\Gamma$", "X", "W", "K",
    #        "$\\Gamma$", "L", "U",
    #        "W", "L", "K",
    #        "U", "X"]

    #qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
    #phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
    #plt.figure(figsize=(20,20))
    #phonon.plot_band_structure()
    #plt.savefig(f"{folder_path}/{mp_id}_phonon_band_structure.png")
    #plt.close()
    #bands_dict = phonon.get_band_structure_dict()
    #qpoints = bands_dict['qpoints']
    #frequencies = bands_dict['frequencies']



class MlPhononsWithDDB(_MlBase):
    """
    Compute phonons with phonopy and a set of ML potentials starting from a DDB file.
    """

    def __init__(self, ddb_filepath: str, supercell, qmesh, asr, nqpath,
                 relax_mode, fmax, pressure, steps, optimizer, nn_name,
                 verbose, workdir, prefix=None):
        """
        Args:
            ddb_filepath: DDB filename.
            supercell: tuple with supercell dimension.
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
        self.ddb = DdbFile(ddb_file)
        self.initial_atoms = get_atoms(ddb.structure)
        self.supercell = supercell
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

     ddb_path   = {self.ddb.filepath}
     supercell  = {self.supercell}
     qmesh      = {self.qmesh}
     asr        = {self.asr}
     nqpath     = {self.nqpath}
     relax_mode = {self.relax_mode}
     fmax       = {self.fmax}
     steps      = {self.steps}
     optimizer  = {self.optimizer}
     pressure   = {self.pressure}
     nn_names   = {self.nn_names}
     workdir    = {self.workdir}
     verbose    = {self.verbose}

=== ATOMS ===

{self.initial_atoms}
"""
        return s

    def run(self) -> None:
        """Run MlPhononsWithDDB."""
        workdir = self.workdir

        atoms = self.initial_atoms.copy()
        calculator = CalcBuilder(self.nn_name).get_calculator()
        atoms.calc = calculator

        if self.relax != RX_MODE.no:
            print(f"Relaxing DDB atoms with relax mode: {self.relax_mode}.")
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

        phonon = get_phonopy(structure, calculator, distance, supercell, primitive_matrix=None)
        #phonon.save(settings={'force_constants': True})

        # Call phonopy to compute phonons with finite difference and ML potential.
        # Include non-analytical term if dipoles are available in the DDB file.
        if self.dbb.has_lo_to_data:
            print("Activating dipolar term in phonopy calculation using BECS and Zeff taken from DDB.")

        # Call anaddb to compute ab-initio phonon frequencies from the DDB.
        #self.ddb.anaget_phmodes_at_qpoints(
        #   qpoints=None, asr=2, chneut=1, dipdip=1, dipquad=1, quadquad=1,
        #   ifcflag=0, ngqpt=None, workdir=None, mpi_procs=1, manager=None, verbose=0,
        #   lo_to_splitting=False,
        #   spell_check=True, directions=None, anaddb_kwargs=None, return_input=False)
