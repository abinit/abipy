# coding: utf-8
"""Work subclasses for the computation of luminiscent properties."""

from .works import Work, PhononWork


class LumiWork(Work):
    """
    This Work implements Fig 1 of https://arxiv.org/abs/2010.00423.

    Client code is responsible for the preparation of the supercell and
    of the GS SCF input files with the fixed electronic occupations associated to the two configurations.
    By default, the work computes four total energies corresponding to the Ag, Ag*, Ae*, Ae configurations.
    Optionally, one can activate the computation of four electronic band structures
    and phonons for Ag and Ae. See docstring of from_scf_inputs for further info.
    """

    @classmethod
    def from_scf_inputs(cls, gs_scf_inp, exc_scf_inp, relax_kwargs, ndivsm=0, nb_extra=10,
                        tolwfr=1e-20, ngqpt=None, manager=None):
        """
        Args:
            gs_scf_inp: |AbinitInput| representing a GS SCF run for the ground-state.
            exc_scf_inp: |AbinitInput| representing a GS SCF run for the excited-state.
            relax_kwargs: Dictonary with input variables to be added to gs_scf_inp and exc_scf_inp.
                when generating input files for structural relaxations.
            ndivsm: Activates the computation of band structure if different from zero.
                if > 0, it's the number of divisions for the smallest segment of the path (Abinit variable).
                if < 0, it's interpreted as the pymatgen `line_density` parameter in which the number of points
                in the segment is proportional to its length. Typical value: -20.
                This option is the recommended one if the k-path contains two high symmetry k-points that are very close
                as ndivsm > 0 may produce a very large number of wavevectors.
            nb_extra: Number of extra bands added to the input nband when computing band structures (ndivsm != 0).
            tolwfr: Tolerance on the residuals used for the NSCF band structure calculations.
            ngqpt: If not None, activates computation of phonons with DFPT using this q-mesh.
                Usually Gamma-only that is ngqpt = [1, 1, 1].
            manager: |TaskManager| of the work. If None, the manager is initialized from the config file.
        """
        new = cls(manager=manager)

        # Templates for GS SCF calculations.
        new.gs_scf_inp = gs_scf_inp
        new.exc_scf_inp = exc_scf_inp

        # Save paramenters for the generation of input files at runtime.
        new.relax_kwargs = relax_kwargs
        new.ndivsm = int(ndivsm)
        new.tolwfr = tolwfr
        new.nb_extra = int(nb_extra)
        new.ngqpt = ngqpt

        # Relaxation for the Ag configuration.
        new.gs_relax_task = new.register_relax_task(gs_scf_inp.new_with_vars(relax_kwargs))

        if new.ngqpt is not None:
            # Make sure we have the WFK file if we are gonna compute phonons.
            new.gs_relax_task.input["prtwf"] = 1

        # Internal counter used in on_all_ok to drive the differ steps of the calculations.
        new.iteration_step = 0

        # JSON-compatible dictionary storing the most important results of the Work.
        # Results will be stored in the `lumi.json` file in the outdata directory of the Work
        # so that one can easily implement additional post-processing tools.
        new.json_data = {}

        return new

    def on_all_ok(self):
        """
        This method is called when all the tasks in the work have reached S_OK.

        This is the section in which we implement most of the workflow logic at runtime.
        since we need to generate input files with relaxed structures.
        """

        if self.iteration_step == 0:
            print("in iteration step 0")
            self.iteration_step += 1

            # Get Ag relaxed structure and energy.
            with self.gs_relax_task.open_gsr() as gsr:
                ag_relaxed_structure = gsr.structure
                self.json_data["Ag_energy_ev"] = float(gsr.energy)
                self.json_data["Ag_gsr_filepath"] = gsr.filepath

            # Build GS SCF input for the Ag* configuration:
            # use same structure as Ag but with excited occupation factors.
            exc_scf_inp = self.exc_scf_inp.new_with_structure(ag_relaxed_structure)
            self.agstar_scf_task = self.register_scf_task(exc_scf_inp)

            if self.ndivsm != 0:
                # Compute band structure for Ag* configuration.
                self.agstar_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                             tolwfr=self.tolwfr, nb_extra=self.nb_extra)

            #if self.ngqpt is not None:
            #    Create new work to compute phonons for Ag configuration.
            #    TODO: Make sure we can connect DFPT with RelaxTask.
            #    self.ae_scf_task["prtwf"] = 1
            #    ph_work = PhononWork.from_scf_task(scf_task, self.ngqpt, is_ngqpt=True, tolerance=None, with_becs=True,
            #                                       ddk_tolerance=None, manager=None)
            #    self.flow.register_work(ph_work)
            #    Save location of final DDB file in json_data.
            #    self.json_data["Ag_ddb_filepath"] = ph_work.outdir.path_in("out_DDB")

            # Relax geometry with excited configuration starting from Ag*.
            self.agstar_relax_task = self.register_relax_task(exc_scf_inp.new_with_vars(self.relax_kwargs))

            return self.postpone_on_all_ok()

        elif self.iteration_step == 1:
            print("in iteration step 1")
            self.iteration_step += 1

            # Get Ag* total energy.
            with self.agstar_scf_task.open_gsr() as gsr:
                self.json_data["Agstar_energy_ev"] = float(gsr.energy)
                self.json_data["Agstar_gsr_filepath"] = gsr.filepath

            # Get Ae* relaxed structure and total energy.
            with self.agstar_relax_task.open_gsr() as gsr:
                aestar_relaxed_structure = gsr.structure
                self.json_data["Aestar_energy_ev"] = float(gsr.energy)
                self.json_data["Aestar_gsr_filepath"] = gsr.filepath

            # Build GS SCF input for the Ae configuration:
            # same structure as Ae* but with gs occupation factors.
            ae_scf_inp = self.gs_scf_inp.new_with_structure(aestar_relaxed_structure)
            self.ae_scf_task = self.register_scf_task(ae_scf_inp)

            if self.ndivsm != 0:
                # Compute band structure for Ae configuration.
                self.ae_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                         tolwfr=self.tolwfr, nb_extra=self.nb_extra)

            if self.ngqpt is not None:
                # Create new work to compute phonons for the Ae configuration. Make sure prtwf == 1
                self.ae_scf_task.input["prtwf"] = 1
                ph_work = PhononWork.from_scf_task(self.ae_scf_task, self.ngqpt, is_ngqpt=True,
                                                   tolerance=None, with_becs=True, ddk_tolerance=None)
                self.flow.register_work(ph_work)
                # Save location of final DDB file in json_data.
                self.json_data["Ae_ddb_filepath"] = ph_work.outdir.path_in("out_DDB")

            return self.postpone_on_all_ok()

        elif self.iteration_step == 2:
            print("in iteration step 2")
            self.iteration_step += 1

            # Get Ag* total energy.
            with self.ae_scf_task.open_gsr() as gsr:
                ae_structure = gsr.structure
                self.json_data["Ae_energy_ev"] = float(gsr.energy)
                self.json_data["Ae_gsr_filepath"] = gsr.filepath

            # Write json file in the outdire of the work
            self.write_json_in_outdir("lumi.json", self.json_data)

            return super().on_all_ok()
