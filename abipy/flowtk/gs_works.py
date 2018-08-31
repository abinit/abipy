# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import print_function, division, unicode_literals, absolute_import

import json
import numpy as np

from .works import Work
from abipy.core.structure import Structure


class EosWork(Work):
    """
    Work for the computation of the Equation of State
    The EOS is obtained by computing E(V) for several volumes around the input V0,
    The initial volumes are obtained by rescaling the input lattice vectors so that
    length proportions and angles are preserved.
    This guess is exact for cubic materials while other Bravais lattices require
    a constant-volume optimization of the cell geometry.

    If lattice_type=="cubic" and atomic positions are fixed by symmetry.
    use can use move_atoms=False to perform standard GS-SCF calculations.
    In all the other cases, E(V) is obtained by relaxing the atomic positions at fixed volume.

    The E(V) points are fitted at the end of the work and the results are saved in the
    `eos_data.json` file produced in the `outdata` directory.
    The file contains the energies, the volumes and the values of V0, B0, B1 obtained
    with different EOS models.
    """

    @classmethod
    def from_scf_input(cls, scf_input, npoints=4, deltap_vol=0.25, ecutsm=0.5, move_atoms=True,
                       manager=None):
        """
        Build a EosWork from an AbinitInput representing a SCF-GS calculation.

        Args:
            scf_input: AbinitInput for SCF-GS used as template to generate the other inputs.
            npoints: Number of volumes generated on the right (left) of the equilibrium volume
                The total number of points is therefore 2 * n + 1.
            deltap_vol: Step of the linear mesh given in relative percentage of the equilibrium volume
                The step is thus: v0 * deltap_vol / 100.
            ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
                value will be used else the vale in `scf_input`.
            move_atoms: If True, a structural relaxation of ions is performed for each volume
                This is needed if the atomic positions are non fixed by symmetry.
            manager: TaskManager instance. Use default if None.

        Return: EosWork instance.
        """
        new_work = cls(manager=manager)

        structure = scf_input.structure
        lattice_type = structure.spget_lattice_type()
        assert lattice_type is not None

        dvol = structure.volume * deltap_vol / 100
        v0 = structure.volume - dvol * npoints
        new_work.input_volumes = [v0 + ipt * dvol for ipt in range(2 * npoints + 1)]

        if "ecutsm" not in scf_input:
            print("Input does not define ecutsm input variable.\n",
                  "A default value of %s will be added to all the EOS inputs" % ecutsm)

        for vol in new_work.input_volumes:
            # Build structure with new volume and generate new input.
            new_lattice = structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            new_input = scf_input.new_with_structure(new_structure)

            # Add ecutsm if not already present.
            new_input.set_vars_ifnotin(ecutsm=ecutsm)

            if lattice_type == "cubic" and not move_atoms:
                # Perform GS calculations without moving atoms, cells do not need to be relaxed.
                new_input.pop_vars(["ionmov", "optcell", "ntime"])
                new_work.register_scf_task(new_input)
            else:
                # Constant-volume optimization of cell geometry + atoms.
                # (modify acell and rprim under constraint - normalize the vectors of rprim to generate the acell)
                # In principle one could take into account the symmetry of the lattice...
                new_input.set_vars_ifnotin(ionmov=2, ntime=50, optcell=3, dilatmx=1.05)
                new_work.register_relax_task(new_input)

        return new_work

    def getnwrite_eosdata(self, write_json=True):
        """
        This method is called when all tasks reach S_OK. It reads the energies
        and the volumes from the GSR file, computes the EOS and produce a
        JSON file `eos_data.json` in outdata.
        """
        energies_ev, volumes = [], []
        for task in self:
            with task.open_gsr() as gsr:
                volumes.append(float(gsr.structure.volume))
                energies_ev.append(float(gsr.energy))

        from pymatgen.analysis.eos import EOS
        eos_data = {"input_volumes_ang3": self.input_volumes,
                    "volumes_ang3": volumes,
                    "energies_ev": energies_ev}

        for model in EOS.MODELS:
            if model in ("deltafactor", "numerical_eos"): continue
            try:
                fit = EOS(model).fit(volumes, energies_ev)
                eos_data[model] = {k: float(v) for k, v in fit.results.items()}
            except Exception as exc:
                eos_data[model] = {"exception": str(exc)}

        if write_json:
            with open(self.outdir.path_in("eos_data.json"), "wt") as fh:
                json.dump(eos_data, fh, indent=4, sort_keys=True)

        return eos_data

    def on_all_ok(self):
        """
        This method is called when all tasks reach S_OK. It reads the energies
        and the volumes from the GSR file, computes the EOS and produce a
        JSON file `eos_data.json` in outdata.
        """
        self.getnwrite_eosdata()
        return dict(returncode=0, message="EOS computed and file written")


#class FiniteDiffStressWork(Work):
#    """
#    Work for the computation of the stress tensor
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
