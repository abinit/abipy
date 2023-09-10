"""
Based on a similar implementation available at: https://github.com/modl-uclouvain/randomcarbon/blob/main/randomcarbon/run/phonon.py
"""
from __future__ import annotations

import numpy as np
import abipy.core.abinit_units as abu

#from monty.json import MSONable, MontyDecoder
from monty.dev import requires
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms
from pymatgen.io.phonopy import get_phonopy_structure
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

    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, primitive_matrix=primitive_matrix)

    phonon.generate_displacements(distance=distance)

    forces = []
    nsc = len(phonon.supercells_with_displacements)
    #print(f"Calculating forces for {nsc} supercell configurations ...")
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


class MlPhononsWithDDB(MlBase):
    """
    Compute phonons with phonopy and a ML potential starting from a DDB file
    and compare the results.
    """
    def __init__(self, ddb_filepath, distance, qmesh, asr, dipdip, nqpath,
                 relax_mode, fmax, pressure, steps, optimizer, nn_name,
                 verbose, workdir, prefix=None, supercell=None):
        """
        Args:
            ddb_filepath: DDB filepath.
            supercell: with supercell dimensions. None to use the supercell from the DDB file.
            distance:
            qmesh: q-mesh for phonon-DOS
            asr: Enforce acoustic sum-rule.
            dipdip: Treatment of dipole-dipole term.
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
        self.distance = float(distance)
        self.qmesh = qmesh
        self.asr = asr
        self.dipdip = dipdip
        self.nqpath = nqpath
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.steps = steps
        self.optimizer = optimizer
        self.nn_name = nn_name
        self.verbose = verbose
        
        with DdbFile(self.ddb_filepath) as ddb:
            self.initial_atoms = get_atoms(ddb.structure)
            self.supercell = supercell
            if self.supercell is None:
                # Take supercell from the q-mesh used in the DDB
                self.supercell = np.eye(3) * ddb.guessed_ngqpt

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     ddb_path   = {self.ddb_filepath}
     supercell  = {self.supercell}
     distance   = {self.distance}
     qmesh      = {self.qmesh}
     asr        = {self.asr}
     dipdip     = {self.dipdip}
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
        self._run_nn_name(self.nn_name)
        self._finalize()

    def _run_nn_name(self, nn_name) -> None:
        """Run MlPhononsWithDDB."""
        workdir = self.workdir

        calculator = CalcBuilder(nn_name).get_calculator()
        atoms = self.initial_atoms.copy()
        atoms.calc = calculator
        natom = len(atoms)

        if self.relax_mode == RX_MODE.cell:
            raise RuntimeError("Öne should not relax the cell when comparing ML phonons with a DDB file!")

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
        print("Calling get_phonopy with nn_name:", nn_name)
        phonon = get_phonopy(atoms, self.supercell, calculator, distance=self.distance, primitive_matrix=None)
        show = True

        with DdbFile(self.ddb_filepath) as ddb:
            if ddb.has_lo_to_data and self.dipdip != 0:
                print("Activating dipolar term in phonopy using BECS and eps_inf taken from DDB.")
                out = ddb.anaget_epsinf_and_becs()
                # according to the phonopy website 14.399652 is not the coefficient for abinit
                # probably it relies on the other conventions in the output.
                self.nac_params = {"born": out.becs.values, "dielectric": out.epsinf, "factor": 14.399652}
                phonon.nac_params = self.nac_params

            # Call anaddb to compute ab-initio phonons from the DDB.
            print("Starting anaddb ph-bands computation...")
            with ddb.anaget_phbst_and_phdos_files(
                    nqsmall=0, qppa=None, ndivsm=20, line_density=None, asr=self.asr, chneut=1,
                    dipdip=self.dipdip, dipquad=1, quadquad=1, verbose=self.verbose) as g:
                phbst_file, phdos_file = g[0], g[1]
                self.abi_phbands = phbst_file.phbands
                py_qpoints = [[qpt.frac_coords for qpt in self.abi_phbands.qpoints]]

            print("Starting phonopy ph-bands computation...")
            phonon.run_band_structure(py_qpoints, with_eigenvectors=False)

            plt = phonon.plot_band_structure()
            if show: plt.show() 
            plt.savefig(workdir / f"phonopy_{nn_name}_phbands.png")
            plt.close()

            bands_dict = phonon.get_band_structure_dict()
            nqpt = 0
            abi_qpoints, py_phfreqs = [], []
            for q_list, w_list in zip(bands_dict['qpoints'], bands_dict['frequencies']): # bands_dict['eigenvectors'])
                nqpt += len(q_list)
                abi_qpoints.extend(q_list)
                py_phfreqs.extend(w_list)

            abi_qpoints = np.reshape(abi_qpoints, (nqpt, 3))
            py_phfreqs = np.reshape(py_phfreqs, (nqpt, 3*natom)) / abu.eV_to_THz

            from abipy.dfpt.phonons import PhononBands, PhononBandsPlotter
            py_phbands = PhononBands(self.abi_phbands.structure, self.abi_phbands.qpoints, py_phfreqs, 
                                     # FIXME
                                     self.abi_phbands.phdispl_cart, 
                                     non_anal_ph=None, 
                                     amu=self.abi_phbands.amu,
                                     epsinf=self.abi_phbands.epsinf,
                                     zcart=self.abi_phbands.zcart,
                                     )
            
            ph_plotter = PhononBandsPlotter(key_phbands=[
                (f"phonopy with {nn_name}", py_phbands),
                ("ABINIT DDB", self.abi_phbands),
            ])
            ph_plotter.combiplot(show=show, savefig=str(workdir / f"combiplot_{nn_name}.png"))

            mabs_diff_mev = np.abs(py_phbands.phfreqs - self.abi_phbands.phfreqs).mean() * 1000
            print("Mean absolute difference in meV:", mabs_diff_mev)

