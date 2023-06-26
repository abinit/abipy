# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

import json

from pymatgen.analysis.eos import EOS
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.electrons.gsr import GsrRobot
from .works import Work

__all__ = [
    "GsKmeshConvWork"
    "GsKmeshTsmearConvWork",
    "EosWork",
]


class GsKmeshConvWork(Work):
    """
    This work performs convergence studies of GS properties
    with respect to the k-mesh

    It produces ...

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsKmeshConvWork
    """

    @classmethod
    def from_scf_input(cls, scf_input: AbinitInput, nksmall_list: list) -> GsKmeshConvWork:
        """
        Build the work from a `scf_input` for a GS SCF run and a list
        with the smallest number of divisions for the k-mesh.
        """
        work = cls()
        for nksmall in nksmall_list:
            new_inp = scf_input.deepcopy()
            new_inp.set_autokmesh(nksmall)
            work.register_scf_task(new_inp)

        return work

    def on_all_ok(self):
        """
        This method is called when all tasks in the GsKmeshTsmearConvWork have reached S_OK.
        """
        with GsrRobot.from_work(self) as gsr_robot:
            df = gsr_robot.get_dataframe(with_geo=False)
            # Write excel file.
            basename = self.__class__.__name__
            df.to_excel(self.outdir.path_in(f"{basename}.xlsx"))

            with gsr_robot.get_pyscript(self.outdir.path_in("gsr_robot.py")) as script:
                script.add_text("""
#item = "energy_per_atom"
#robot.plot_convergence(item, sortby="nkpt", abs_conv=1e-3)

items = ["energy_per_atom", "pressure", "max_force"]
robot.plot_convergence_items(items, sortby="nkpt")

abs_conv = {
"energy_per_atom": 1e-3,
"pressure": 1e-2,
"max_force": 1e-4,
}
items = abs_conv.keys()
robot.plot_convergence_items(items, sortby="nkpt", abs_conv=abs_conv)
""")

        return super().on_all_ok()



class GsKmeshTsmearConvWork(Work):
    """
    This work performs convergence studies of GS properties
    with respect to the k-mesh and the electronic smearing.

    It produces ...

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsKmeshTsmearConvWork
    """

    @classmethod
    def from_scf_input(cls, scf_input: AbinitInput, nksmall_list: list, tsmear_list: list) -> GsKmeshTsmearConvWork:
        """
        Build the work from a `scf_input` for a GS SCF run including `occopt`
        and a list with the smallest number of divisions for the k-mesh.
        """
        occopt = scf_input.get("occopt", default=None)
        if occopt is None or occopt <= 0:
            raise ValueError(f"scf_input should define occopt but found: {occopt}")

        work = cls()
        for tsmear in tsmear_list:
            for nksmall in nksmall_list:
                new_inp = scf_input.new_with_vars(tsmear=tsmear)
                new_inp.set_autokmesh(nksmall)
                work.register_scf_task(new_inp)

        return work

    def on_all_ok(self):
        """
        This method is called when all tasks in the GsKmeshTsmearConvWork have reached S_OK.
        """
        with GsrRobot.from_work(self) as gsr_robot:
            df = gsr_robot.get_dataframe(with_geo=False)
            # Write excel file.
            basename = self.__class__.__name__
            df.to_excel(self.outdir.path_in(f"{basename}.xlsx"))

            with gsr_robot.get_pyscript(self.outdir.path_in("gsr_robot.py")) as script:
                script.add_text("""
#item = "energy_per_atom"
#robot.plot_convergence(item, sortby="nkpt", abs_conv=1e-3)

items = ["energy_per_atom", "pressure", "max_force"]
robot.plot_convergence_items(items, sortby="nkpt")

abs_conv = {
"energy_per_atom": 1e-3,
"pressure": 1e-2,
"max_force": 1e-4,
}
items = abs_conv.keys()
robot.plot_convergence_items(items, sortby="nkpt", hue="tsmear", abs_conv=abs_conv)
""")

        return super().on_all_ok()


class EosWork(Work):
    """
    Work to compute the Equation of State.
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

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EosWork
    """

    @classmethod
    def from_scf_input(cls, scf_input: AbinitInput,
                       npoints=4, deltap_vol=0.25, ecutsm=0.5, move_atoms=True,
                       manager=None) -> EosWork:
        """
        Build an EosWork from an AbinitInput for GS-SCF.

        Args:
            scf_input: AbinitInput for GS-SCF used as template to generate the other inputs.
            npoints: Number of volumes generated on the right (left) of the equilibrium volume
                The total number of points is therefore 2 * n + 1.
            deltap_vol: Step of the linear mesh given in relative percentage of the equilibrium volume
                The step is thus: v0 * deltap_vol / 100.
            ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
                value will be used else the vale in `scf_input`.
            move_atoms: If True, a structural relaxation of ions is performed for each volume
                This is needed if the atomic positions are non fixed by symmetry.
            manager: TaskManager instance. Use default if None.
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

    @classmethod
    def from_inputs(cls, inputs: list, manager=None) -> EosWork:
        """
        Advanced interface to build an EosWork from an list of AbinitInputs.
        """
        new_work = cls(manager=manager)
        new_work.input_volumes = [inp.structure.lattice.volume for inp in inputs]
        for scf_inp in inputs:
            new_work.register_scf_task(scf_inp)

        return new_work

    def get_and_write_eosdata(self, write_json=True) -> dict:
        """
        Compute the EOS and produce a JSON file `eos_data.json` in outdata.
        """
        energies_ev, volumes = [], []
        for task in self:
            with task.open_gsr() as gsr:
                energies_ev.append(float(gsr.energy))
                volumes.append(float(gsr.structure.volume))

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
        This method is called when all tasks have reached S_OK.
        It reads the energies and the volumes from the GSR file, computes the EOS
        and produce a JSON file `eos_data.json` in outdata.
        """
        self.get_and_write_eosdata()
        return super().on_all_ok()
