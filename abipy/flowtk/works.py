# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from pymatgen.io.abinit.works import Work


class NscfDdksWork(Work):
    """
    This work requires a DEN file and computes the KS energies with a non self-consistent task
    with a dense k-mesh and empty states.
    This task is then followed by the computation of the DDK matrix elements with nstep = 1
    (the first order change of the wavefunctions is not converged but we only need the matrix elements)
    Mainly used to prepare optic calculations or other post-processing steps requiring the DDKs.
    """

    @classmethod
    def from_scf_task(cls, scf_task, ddk_ngkpt, ddk_shiftk, ddk_nband, manager=None):
        """
        Args:
            scf_task: GS task. Must produce the DEN file required for the NSCF run.
            ddk_ngkpt: k-mesh used for the NSCF run and the non self-consistent DDK tasks.
            ddk_shiftk: k-mesh shifts
            ddk_nband: Number of bands (occupied + empty) used in the NSCF task and the DDKs tasks.
            manager: TaskManager instance. Use default if None.

        Return: NscfDdksWork instance
        """
        new = cls(manager=manager)

        # NSCF task with nband states and points in the IBZ (note kptopt = 1)
        nscf_inp0 = scf_task.input.deepcopy()
        nscf_inp0.set_vars(nband=ddk_nband, prtwf=1)
        nscf_inp0.set_kmesh(ddk_ngkpt, ddk_shiftk, kptopt=1)
        nscf_task0 = new.register_nscf_task(nscf_inp0, deps={scf_task: "DEN"})

        # NSCF run with nband states and points in the IBZ defined by time-reversal only (as required by DDK)
        # This is gonna be quick because Abinit will symmetrize states from the previous WFK file.
        # Time-reversal symmetry can be used in optic.
        #nscf_inp1 = nscf_inp0.deepcopy()
        #nscf_inp0.set_kmesh(ddk_ngkpt, ddk_shiftk, kptopt=2)
        #nscf_task1 = new.register_nscf_task(nscf_inp1)

        # This is the task producing the KS energies for optic
        new.task_with_ks_energies = nscf_task0

        # Build task for one-shot DDKs (note kptopt 2)
        ddk_inputs = nscf_inp0.make_ddk_inputs(kptopt=2)
        new.ddk_tasks = []
        for ddk_inp in ddk_inputs:
            # FIXME: prtwfk should be set to 0 but need to replace DDK.nc
            ddk_inp.set_vars(nstep=1, nline=0, prtwf=1)
            #new.register_ddk_task(ddk_inp, deps={nscf_task0: "WFK"})
            # FIXME: Here I have a conflict with DDK.nc and DDK
            t = new.register_task(ddk_inp, deps={nscf_task0: "WFK"})
            new.ddk_tasks.append(t)

        return new


class EosWork(Work):
    """
    Work for the computation of the Equation of State
    The EOS is obtained by computing E(V) for several volumes around the input V0,
    The initial volumes are obtained by rescaling the input lattice vectors so that
    length proportions and angles are preserved.
    This guess is exact for cubic/rhomboedral materials while other Bravais lattices require
    a constant-volume optimization of the cell geometry.

    If lattice_type in ("cubic", "rhomboedral") and atomic positions are fixed by symmetry.
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

        from abipy.core.structure import Structure
        for vol in new_work.input_volumes:
            # Build structure with new volume and generate new input.
            new_lattice = structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            new_input = scf_input.new_with_structure(new_structure)

            # Add ecutsm if not already present.
            new_input.set_vars_ifnotin(ecutsm=ecutsm)

            if lattice_type in ("cubic", "rhomboedral") and not move_atoms:
                # Perform GS calculations without moving atoms, cells do not need to be relaxed.
                new_input.pop_vars(["ionmov", "optcell", "ntime"])
                new_work.register_scf_task(new_input)
            else:
                # Constant-volume optimization of cell geometry + atoms.
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
            import json
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
