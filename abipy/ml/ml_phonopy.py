"""

"""
from __future__ import annotations

import numpy as np
import json
import abipy.core.abinit_units as abu

#from monty.json import MSONable, MontyDecoder
from monty.dev import requires
from monty.string import marquee, list_strings
from monty.termcolor import cprint
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms
from pymatgen.io.phonopy import get_phonopy_structure
from abipy.core.structure import Structure
from abipy.dfpt.ddb import DdbFile
from abipy.tools.context_managers import Timer
from abipy.ml.aseml import RX_MODE, CalcBuilder, AseResults, MlBase, get_atoms, relax_atoms, dataframe_from_results_list
try:
    import phonopy
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except ImportError:
    phonopy = None
    Phonopy = None



def cprint_traceback() -> None:
    """Returns a string with the traceback."""
    import traceback
    cprint(traceback.format_exc(), color="red")


@requires(phonopy, "phonopy should be installed to calculate phonons")
def get_phonopy(structure: Structure,
                supercell_matrix,
                calculator: Calculator,
                distance=0.01,
                primitive_matrix=None,
                remove_drift=False,
                ) -> Phonopy:
    """
    Args:
        structure:
        supercell_matrix
        calculator:
        distance:
        primitive_matrix
        remove_drift:

    Based on a similar implementation available at: https://github.com/modl-uclouvain/randomcarbon/blob/main/randomcarbon/run/phonon.py
    """
    unitcell = Structure.as_structure(structure).get_phonopy_atoms()

    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, primitive_matrix=primitive_matrix)
    phonon.generate_displacements(distance=distance)

    forces_list = []
    nsc = len(phonon.supercells_with_displacements)
    #print(f"Calculating forces for {nsc} supercell configurations ...")

    for i, sc in enumerate(phonon.supercells_with_displacements):
        #print(f"{i+1} of {nsc}")
        a = Atoms(symbols=sc.symbols,
                  positions=sc.positions,
                  masses=sc.masses,
                  cell=sc.cell, 
                  pbc=True,
                  calculator=calculator)

        forces = a.get_forces()
        if remove_drift:
            drift_force = forces.sum(axis=0)
            for force in forces:
                force -= drift_force / forces.shape[0]
        forces_list.append(forces)

    phonon.produce_force_constants(forces_list)
    return phonon


class MlPhonopyWithDDB(MlBase):
    """
    Compute phonons with phonopy and a ML potential starting from a DDB file
    and compare the results.
    """
    def __init__(self, distance, asr, dipdip, line_density, qppa,
                 relax_mode, fmax, pressure, steps, optimizer, nn_names,
                 verbose, workdir, prefix=None, 
                 ddb_filepath=None, supercell=None):
        """
        Args:
            distance:
            asr: Enforce acoustic sum-rule. Abinit variable.
            dipdip: Treatment of dipole-dipole term. Abinit variable.
            line_density: Defines the a density of k-points per reciprocal atom to plot the phonon dispersion.
            qppa: Defines the homogeneous q-mesh used for the DOS in units of q-points per atom.
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
            pressure: Target pressure.
            steps: max number of steps for relaxation.
            optimizer: name of the ASE optimizer to use for relaxation.
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
            workdir: Working directory, None to generate temporary directory automatically.
            prefix: Prefix for workdir
            ddb_filepath: DDB filepath.
            supercell: with supercell dimensions. None to use the supercell from the DDB file.
        """
        super().__init__(workdir, prefix)

        self.distance = float(distance)
        self.asr = asr
        self.dipdip = dipdip
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.steps = steps
        self.optimizer = optimizer
        self.nn_names = list_strings(nn_names)
        self.verbose = verbose
        
        self.supercell = supercell
        self.ddb_filepath = ddb_filepath

        with DdbFile(self.ddb_filepath) as ddb:
            self.initial_atoms = get_atoms(ddb.structure)
            if self.supercell is None:
                # Take supercell from the q-mesh used in the DDB
                self.supercell = np.eye(3) * ddb.guessed_ngqpt

        if self.supercell is None:
            raise ValueError("Supercell must be specified in input!")

        self.line_density = line_density
        self.qppa = qppa
        self.abi_nac_params = {}

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        supercell = self.supercell.tolist()
        s = f"""\

{self.__class__.__name__} parameters:

     ddb_path      = {self.ddb_filepath}
     supercell     = {supercell}
     distance      = {self.distance}
     qppa          = {self.qppa}
     asr           = {self.asr}
     dipdip        = {self.dipdip}
     line_density  = {self.line_density}
     relax_mode    = {self.relax_mode}
     fmax          = {self.fmax}
     steps         = {self.steps}
     optimizer     = {self.optimizer}
     pressure      = {self.pressure}
     nn_names      = {self.nn_names}
     workdir       = {self.workdir}
     verbose       = {self.verbose}

=== ATOMS ===

{self.initial_atoms}
"""
        return s

    def run(self) -> None:
        """Run MlPhonopyWithDDB."""

        # Run anaddb computation and save results in self.
        with DdbFile(self.ddb_filepath) as ddb, Timer(header="Starting anaddb ph-bands computation...") as timer:
            if ddb.has_lo_to_data and self.dipdip != 0:
                # according to the phonopy website 14.399652 is not the coefficient for abinit
                # probably it relies on the other conventions in the output.
                out = ddb.anaget_epsinf_and_becs()
                self.abi_nac_params = {"born": out.becs.values, "dielectric": out.epsinf, "factor": 14.399652}

            # ab-initio phonons from the DDB.
            with ddb.anaget_phbst_and_phdos_files(
                        qppa=self.qppa, line_density=self.line_density, 
                        asr=self.asr, dipdip=self.dipdip, verbose=self.verbose) as g:

                phbst_file, phdos_file = g[0], g[1]
                self.abi_phbands = phbst_file.phbands
                # The q-points passed to phonopy.
                self.py_qpoints = [[qpt.frac_coords for qpt in self.abi_phbands.qpoints]]

        data = {}
        data["ddb"] = self.abi_phbands.get_phfreqs_stats_dict()
        d_nn = data["nn_names"] = {}
        for nn_name in self.nn_names:
            try:
                d = self._run_nn_name(nn_name)
                d_nn[nn_name] = d
                d_nn[nn_name]["exception"] = None
            except Exception as exc:
                cprint_traceback()
                d_nn[nn_name] = {}
                d_nn[nn_name]["exception"] = str(exc)

        self.write_json("data.json", data, info="JSON file with final results.")
        self._finalize()

    def _run_nn_name(self, nn_name) -> None:
        """
        Run calculation for a single NN potential.
        """
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
                             traj_path=self.get_path(f"{nn_name}_relax.traj", "ASE relax trajectory"),
                             verbose=self.verbose,
                            )

            relax = relax_atoms(atoms, **relax_kws)
            relax.summarize(["DDB_initial", "DDB_relaxed"])

        # Call phonopy to compute phonons with finite difference and ML potential.
        # Include non-analytical term if dipoles are available in the DDB file.
        with Timer(header=f"Calling get_phonopy with {nn_name=}", footer=""):
            phonon = get_phonopy(atoms, self.supercell, calculator, distance=self.distance, primitive_matrix=None)
            if self.abi_nac_params:
                print("Including dipolar term in phonopy using BECS and eps_inf taken from DDB.")
                phonon.nac_params = self.abi_nac_params

        force_constants = phonon.get_force_constants()
        #from phonopy.file_IO import write_FORCE_CONSTANTS
        #write_FORCE_CONSTANTS(force_constants, filename=str(self.workdir / f"{nn_name}_FORCE_CONSTANTS"))

        """
        from abipy.dfpt.converters import phonopy_to_abinit
        with DdbFile(self.ddb_filepath) as ddb:
            zion = ddb.header["zion"]
            py_ddb = phonopy_to_abinit(unit_cell=Structure.as_structure(atoms.copy()), 
                                       supercell_matrix=self.supercell, 
                                       out_ddb_path=str(self.workdir / f"{nn_name}_DDB"),
                                       ngqpt=ddb.guessed_ngqpt,
                                       force_constants=force_constants, 
                                       born=None,
                                       primitive_matrix=np.eye(3), 
                                       verbose=self.verbose,
                                       )
            # Call anaddb to compute ab-initio phonons from the DDB.
            with py_ddb.anaget_phbst_and_phdos_files(
                        qppa=self.qppa, line_density=self.line_density, 
                        asr=self.asr, dipdip=self.dipdip, verbose=self.verbose) as g:

                phbst_file, phdos_file = g[0], g[1]
                py_phbands = phbst_file.phbands

            py_ddb.close()
        """

        show = False
        from abipy.dfpt.phonons import PhononBands, PhononBandsPlotter

        with Timer(header="Starting phonopy ph-bands computation...", footer=""):
            phonon.run_band_structure(self.py_qpoints, with_eigenvectors=True)

        plt = phonon.plot_band_structure()
        plt.savefig(workdir / f"phonopy_{nn_name}_phbands.png")
        if show: plt.show() 
        plt.close()

        bands_dict = phonon.get_band_structure_dict()
        nqpt = 0
        py_phfreqs, py_displ_cart = [], []
        for q_list, w_list, eig_list in zip(bands_dict['qpoints'], bands_dict['frequencies'], bands_dict['eigenvectors']):
            nqpt += len(q_list)
            py_phfreqs.extend(w_list)
            #print(eig_list)
            py_displ_cart.extend(eig_list)

        py_phfreqs = np.reshape(py_phfreqs, (nqpt, 3*natom)) / abu.eV_to_THz
        py_displ_cart = np.reshape(py_displ_cart, (nqpt, 3*natom, 3*natom))

        # Build abipy phonon bands from phonopy results.
        py_phbands = PhononBands(self.abi_phbands.structure, self.abi_phbands.qpoints, py_phfreqs, 
                                 # FIXME: Use phononopy displacement
                                 self.abi_phbands.phdispl_cart, 
                                 non_anal_ph=None, 
                                 amu=self.abi_phbands.amu,
                                 epsinf=self.abi_phbands.epsinf,
                                 zcart=self.abi_phbands.zcart,
                                 )

        # Compute diff stats.
        mabs_wdiff_ev = np.abs(py_phbands.phfreqs - self.abi_phbands.phfreqs).mean()
        
        ph_plotter = PhononBandsPlotter(key_phbands=[
            (f"phonopy with {nn_name}", py_phbands),
            ("ABINIT DDB", self.abi_phbands),
        ])
        mae_str = f"MAE {1000 * mabs_wdiff_ev:.3f} meV"
        print(mae_str)
        latex_formula = self.abi_phbands.structure.latex_formula
        ph_plotter.combiplot(show=show, title=f"{latex_formula}: {mae_str}",
                             savefig=str(workdir / f"combiplot_{nn_name}.png"))

        data = dict(
            mabs_wdiff_ev=mabs_wdiff_ev,
            **py_phbands.get_phfreqs_stats_dict()
        )

        return data