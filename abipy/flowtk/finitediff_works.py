# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

import json
import itertools
#import pickle
import dataclasses
import numpy as np

#from monty.json import MSONable
#from monty.string import list_strings #, marquee
from pymatgen.analysis.eos import EOS
from abipy.core.structure import Structure
from abipy.tools.numtools import build_mesh
from abipy.tools.derivatives import central_fdiff_weights # finite_diff
from abipy.abio.inputs import AbinitInput
from abipy.tools.serialization import HasPickleIO
#from abipy.electrons.gsr import GsrRobot
#from abipy.tools.serialization import mjson_write, pmg_serialize
from .works import Work


class FiniteDiffForcesData(HasPickleIO):

    def __init__(self,
                 iatom,
                 step,
                 mesh_type,
                 deltas,
                 ix0,
                 energies_ev,
                 structures,
                 cart_forces_list,
                 cart_stress_tensor_list,
                 ):

        self.iatom = iatom
        self.step = step
        self.mesh_type = mesh_type
        self.deltas = deltas
        self.ix0 = ix0
        self.energies_ev = energies_ev
        self.structures = structures
        self.cart_forces_list = cart_forces_list
        self.cart_stress_tensor_list = cart_stress_tensor_list

    def get_results_nn_name(self, nn_name: str, with_delta: bool) -> dict:
        """
        Args:
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
        """
        from abipy.ml.aseml import CalcBuilder
        calc = CalcBuilder(nn_name).get_calculator()
        energies_ev, forces_list, stress_list = [], [], []
        for structure in self.structures:
            atoms = structure.to_ase_atoms(calc=calc)
            energies_ev.append(float(atoms.get_potential_energy()))
            forces_list.append(atoms.get_forces())
            stress_list.append(atoms.get_stress(voigt=False))

        return dict(energies_ev=np.array(energies_ev),
                    forces_list=np.array(forces_list),
                    stress_list=np.array(stress_list),
                    )

    #def compare_with_nn_names(self, nn_names):
    #    for nn_name in list_strings(nn_names):
    #        d = self.get_results_nn_name(nn_name)



class FiniteDiffForcesWork(Work):
    """
    Work for the computation of forces with finite differences.
    """

    @classmethod
    def from_scf_input(cls, scf_input, iatom, direction,
                       frac_coords=True, num_points=5, step=0.01, mesh_type="centered", manager=None):
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            iatom: Index of the atom to displace.
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                Cartesian coordinates.
            num_points:
            step:
            mesh_type: Generate a linear mesh of step `step` that is centered on x0 if
                mesh_type == "centered" or a mesh that starts/ends at x0 if mesh_type is `>`/`<`.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)

        work.iatom = int(iatom)
        work.step = float(step)
        work.mesh_type = mesh_type
        work.deltas, work.ix0 = build_mesh(0.0, num_points, step, mesh_type)

        structure = scf_input.structure
        norm = structure.lattice.norm(direction, frac_coords=frac_coords)
        versor = np.array(direction) / norm

        for delta in work.deltas:
            vector = delta * versor
            new_structure = structure.copy()
            new_structure.translate_sites([work.iatom], delta * versor, frac_coords=frac_coords, to_unit_cell=False)
            new_input = scf_input.new_with_structure(new_structure)
            work.register_scf_task(new_input)

        return work

    def get_data(self) -> FiniteDiffForcesData:
        """
	    Read data from the GSR files, and produce a JSON file
        in the outdata directory of the work.
        """
        energies_ev, cart_forces_list, cart_stress_tensor_list = [], [], []
        for task in self:
            with task.open_gsr() as gsr:
                energies_ev.append(float(gsr.energy))
                cart_forces_list.append(np.array(gsr.cart_forces))
                cart_stress_tensor_list.append(np.array(gsr.cart_stress_tensor))

        return FiniteDiffForcesData(
            iatom=self.iatom,
            step=self.step,
            mesh_type=self.mesh_type,
            deltas=self.deltas,
            ix0=self.ix0,
            energies_ev=energies_ev,
            structures=[task.input.structure for task in self],
            cart_forces_list=cart_forces_list,
            cart_stress_tensor_list=cart_stress_tensor_list,
        )

    def on_all_ok(self):
        """
        This method is called when all tasks reach S_OK.
        """
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        #mjson_write(data, self.outdir.path_in("forces.json"))
        return super().on_all_ok()


#class FiniteDiffStressWork(Work):
#    """
#    Work for the computation of the stress tensor with finite differences
#    """
#
#    @classmethod
#    def from_scf_input(cls, scf_input, delta=1e-4, ecutsm=0.5, manager=None):
#        """
#        Build a EosWork from an AbinitInput representing a GS SCF calculation.
#
#        Args:
#            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
#	        delta:
#            ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
#                value will be used else the vale in `scf_input`.
#            manager: TaskManager instance. Use default if None.
#
#        Return: EosWork instance.
#        """
#	    scf_input = scf_input.deepcopy()
#
#        if "ecutsm" not in scf_input:
#            scf_input.set_vars(ecutsm=ecutsm)
#            print("Input does not define ecutsm input variable.\n",
#                  "A default value of %s will be added to all the EOS inputs" % ecutsm)
#
#        new_work = cls(manager=manager)
#        self.unstrained_task = new_work.register_scf_task(scf_input)
#
#        voigt_inds = [(0, 0)]
#
#        self.delta = float(delta)
#        self.tasks_voigt = {}
#        for voigt in voigt_inds:
#            self.tasks_voigt[voigt] = defaultdict(list)
#            for isign in (-1, +1):
#            # Apply strain to the lattice.
#                    strain = np.zeros((3, 3))
#                    strain[voigt] = float(isign) * delta
#                    new_structure = scf_input.structure.deepcopy()
#                    new_structure.apply_strain(strain)
#                    new_input = scf_input.new_with_structure(new_structure)
#
#                    # Perform GS calculations with strained cell.
#                    task = new_work.register_scf_task(new_input)
#                    self.tasks_voigt[voigt].append(task)
#
#            return new_work
#
#    def get_data(self):
#        """
#        This method is called when all tasks reach S_OK.
#	     It reads the energies and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#        with self.unstrained_task.open_gsr() as gsr0:
#            e0, v0 = gsr0.energy, gsr0.structure.volume
#            cart_stress_tensor = gsr0.cart_stress_tensor
#
#        fd_tensor = np.empty((3, 3))
#        for voigt in voigt_inds:
#            tasks = self.tasks_voigt[voigt]
#            energies_ev, volumes = np.empty(len(tasks), np.empty(len(tasks))
#            for i, task in enumerate(task):
#                    with task.open_gsr() as gsr:
#                       energies_ev[i] = float(gsr.energy))
#
#            d = (e_plus - e_minus) / (2 * self.delta * v0)
#                fd_tensor[voigt] = d
#                fd_tensor[voigt[1], voigt[0]] = d
#
#        data = {
#                 "dfpt_cart_stress_tensor": cart_stress_tensor,
#                 "finite_diff_cart_stress_tensor": fd_cart_stress_tensor,
#            }
#
#            with open(self.outdir.path_in("stress.json"), "wt") as fh:
#                json.dump(data, fh, indent=4, sort_keys=True)
#
#            return data
#
#    def on_all_ok(self):
#        """
#        This method is called when all tasks reach S_OK. It reads the energies
#        and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#        self.compute_and_write_stress()
#        return dict(returncode=0, message="EOS computed and file written")


class FdDynMagneticChargeWork(Work):
    r"""
    Work for the computation of the dynamical magnetic charges with finite differences.

    The dynamical magnetic charges are defined as:

        Z_jv^m=Ω_0 (∂M_v)/(∂u_j ) = (∂F_j)/(∂H_v ) = Ω_0 (∂^2 E)/(∂H_β ∂u_i ).#(6)
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       berryopt: int,
                       num_points: int,
                       delta_h: float = 0.01,
                       relax_opts: dict | None = None,
                       manager=None) -> FdDynMagneticChargeWork:
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF calculation used as template to generate the other inputs.
            berryopt: Abinit input variable.
            num_points:
            delta_h: Finite difference step for the magnetic field in a.u.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.berryopt = berryopt
        work.num_points = num_points
        work.delta_h = delta_h
        work.h_values, work.ib0 = build_mesh(0.0, num_points, delta_h, "centered")
        work.h_cart_dirs = np.array([
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])

        work.scf_input_template = scf_input.deepcopy()
        relax_input = scf_input.make_relax_input()
        work.initial_relax_task = work.register_relax_task(relax_input)

        return work

    def on_ok(self, sender):
        """
        This method is called when one task reaches status `S_OK`.
        """
        if sender == self.initial_relax_task:
            # Get relaxed structure
            self.relaxed_structure = sender.get_final_structure()
            scf_input = self.scf_input_template.new_with_structure(self.relaxed_structure)
            scf_input["berryopt"] = self.berryopt

            # Build new GS tasks with zeemanfield.
            shape = len(self.h_cart_dirs), len(self.h_values)
            self.tasks_hdir_h = np.empty(shape, dtype=object)
            task_h0 = None
            for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
                for ih, h_val in enumerate(self.h_values):
                    is_h0 = abs(h_val) < 1e-16
                    new_inp = scf_input.new_with_vars(zeemanfield=h_val * h_cart_dir)
                    if is_h0:
                        new_inp["berryopt"] = 0
                        # Avoid computing H=0 3 times.
                        if task_h0 is None:
                            task_h0 = self.register_scf_task(new_inp)
                        self.tasks_hdir_h[hdir, ih] = task_h0
                    else:
                        self.tasks_hdir_h[hdir, ih] = self.register_scf_task(new_inp)

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def get_data(self):
        """
	    Read data from the GSR files, and produce a JSON file in the outdata directory of the work.
        """
        natom = len(self[0].input.structure)

        data = {
            "input_structure": self[0].input.structure,
            "relaxed_structure": self.relaxed_structure,
            "delta_h": self.delta_h,
            "h_values": self.h_values,
            "h_cart_dirs": self.h_cart_dirs,
        }

        # Read forces from the GSR files.
        ndirs = len(self.h_cart_dirs)
        cart_forces_dh = np.empty((ndirs, self.num_points, natom, 3))
        for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
            for ih, h_val in enumerate(self.h_values):
                with work.tasks_hdir_h[hdir, ih].open_gsr() as gsr:
                    cart_forces_dh[hdir, ih] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
        data["cart_forces_dh"] = cart_forces_dh

        # Finite difference: ∂F_i/∂H_β
        # Use all stencils compatible with input num_points so that we can monitor the convergence.
        data["zm_adh_npts"] = {}
        for acc, weights in central_fdiff_weights[1].items():
            if self.num_points < len(weights): continue
            np = acc // 2
            zm_adh = np.zeros((natom, 3, 3))
            for hdir, iat_idir, iat in itertools.product(range(3), range(3), range(len(natom))):
                fvals_h = cart_forces_dh[hdir, :, iat, iat_dir]
                zm_adh[iat, iat_dir, hdir] = np.sum(fvals_h[self.ib0-np:self.ib0+np] * weights) / self.delta_h
            data["zm_adh_npts"][len(weights)] = zm_adh

        return DynMagCharges(**data)
        #return data
        #mjson_write(data, self.outdir.path_in("dynmag_charges.json"))


@dataclasses.dataclass(kw_only=True)
class DynMagCharges:

    input_structure: Structure
    relaxed_structure: Structure
    delta_h: float
    h_values: np.ndarray
    h_cart_dirs: np.ndarray
    cart_forces: np.ndarray
    zm_adh_npts: np.array

