# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

import json

from pymatgen.analysis.eos import EOS
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.electrons.gsr import GsrRobot
from .works import Work

#class FiniteDiffStressWork(Work):
#    """
#    Work for the computation of the stress tensor with finite difference
#    """
#
#    @classmethod
#    def from_scf_input(cls, scf_input, delta=1e-4, ecutsm=0.5, manager=None):
#        """
#        Build a EosWork from an AbinitInput representing a SCF-GS calculation.
#
#        Args:
#            scf_input: AbinitInput for SCF-GS used as template to generate the other inputs.
#	    delta:
#            ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
#                value will be used else the vale in `scf_input`.
#            manager: TaskManager instance. Use default if None.
#
#        Return: EosWork instance.
#        """
#	scf_input = scf_input.deepcopy()
#
#        if "ecutsm" not in scf_input:
#            scf_input.set_vars(ecutsm=ecutsm)
#            print("Input does not define ecutsm input variable.\n",
#                  "A default value of %s will be added to all the EOS inputs" % ecutsm)
#
#        new_work = cls(manager=manager)
#        self.unstrained_task = new_work.register_scf_task(scf_input)
#
#	voigt_inds = [(0, 0)]
#
#	self.delta = float(delta)
#	self.tasks_voigt = {}
#	for voigt in voigt_inds:
#	    self.tasks_voigt[voigt] = defaultdict(list)
#	    for isign in (-1, +1):
#		# Apply strain to the lattice.
#                strain = np.zeros((3, 3))
#                strain[voigt] = float(isign) * delta
#                new_structure = scf_input.structure.deepcopy()
#                new_structure.apply_strain(strain)
#                new_input = scf_input.new_with_structure(new_structure)
#
#                # Perform GS calculations with strained cell.
#                task = new_work.register_scf_task(new_input)
#                self.tasks_voigt[voigt].append(task)
#
#        return new_work
#
#    def compute_and_write_stress(self):
#        """
#        This method is called when all tasks reach S_OK.
#	It reads the energies and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#	with self.unstrained_task.open_gsr() as gsr0:
#	    e0, v0 = gsr0.energy, gsr0.structure.volume
#	    cart_stress_tensor = gsr0.cart_stress_tensor
#
#	fd_tensor = np.empty((3, 3))
#	for voigt in voigt_inds:
#	    tasks = self.tasks_voigt[voigt]
#	    energies_ev, volumes = np.empty(len(tasks), np.empty(len(tasks))
#	    for i, task in enumerate(task):
#                with task.open_gsr() as gsr:
#                   energies_ev[i] = float(gsr.energy))
#
#	    d = (e_plus - e_minus) / (2 * self.delta * v0)
#            fd_tensor[voigt] = d
#            fd_tensor[voigt[1], voigt[0]] = d
#
#	data = {
#             "dfpt_cart_stress_tensor": cart_stress_tensor,
#             "finite_diff_cart_stress_tensor": fd_cart_stress_tensor,
#        }
#
#        with open(self.outdir.path_in("stress.json"), "wt") as fh:
#            json.dump(data, fh, indent=4, sort_keys=True)
#
#        return data
#
#    def on_all_ok(self):
#        """
#        This method is called when all tasks reach S_OK. It reads the energies
#        and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#        self.compute_and_write_stress()
#        return dict(returncode=0, message="EOS computed and file written")