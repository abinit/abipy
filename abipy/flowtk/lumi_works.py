# coding: utf-8
"""Work subclasses for the computation of luminiscent properties."""
from __future__ import annotations

from .works import Work
from abipy.abilab import abiopen
from abipy.lumi.deltaSCF import DeltaSCF

class LumiWork(Work):
    """
    This Work implements Fig 1 of https://arxiv.org/abs/2010.00423.

    Client code is responsible for the preparation of the supercell and
    of the GS SCF input files with the fixed electronic occupations associated to the two configurations.
    By default, the work computes the two relaxed structures and the four total energies
    corresponding to the Ag, Ag*, Ae*, Ae configurations. Optionally, one can activate the computation of
    four electronic band structures. See docstring of from_scf_inputs for further info.
    """

    @classmethod
    def from_scf_inputs(cls, gs_scf_inp, ex_scf_inp, relax_kwargs_gs, relax_kwargs_ex, ndivsm=0, nb_extra=10, 
                        tolwfr=1e-12, four_points=True, meta=None, manager=None) -> LumiWork:
        """
        Args:
            gs_scf_inp: |AbinitInput| representing a GS SCF run for the ground-state.
            ex_scf_inp: |AbinitInput| representing a GS SCF run for the excited-state.
            relax_kwargs_gs: Dictonary with input variables to be added to gs_scf_inp
                when generating input files for ground state structural relaxations.
            relax_kwargs_ex: Dictonary with input variables to be added to ex_scf_inp
                when generating input files for excited state structural relaxations.
            ndivsm: Activates the computation of band structure if different from zero.
                if > 0, it's the number of divisions for the smallest segment of the path (Abinit variable).
                if < 0, it's interpreted as the pymatgen `line_density` parameter in which the number of points
                in the segment is proportional to its length. Typical value: -20.
                This option is the recommended one if the k-path contains two high symmetry k-points that are very close
                as ndivsm > 0 may produce a very large number of wavevectors.
            nb_extra: Number of extra bands added to the input nband when computing band structures (ndivsm != 0).
            tolwfr: Tolerance of the residuals used for the NSCF band structure calculations.
            four_points : if True, compute the two relaxations and the four points energies.
                If false, only the two relaxations.
            meta : dict corresponding to the metadata of a lumiwork (supercell size, dopant type,...)
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
        """
        new = cls(manager=manager)

        # Templates for GS SCF calculations.
        new.gs_scf_inp = gs_scf_inp
        new.ex_scf_inp = ex_scf_inp

        # Save paramenters for the generation of input files at runtime.
        new.relax_kwargs_gs = relax_kwargs_gs
        new.relax_kwargs_ex = relax_kwargs_ex

        new.four_points = four_points
        new.ndivsm = int(ndivsm)
        new.tolwfr = tolwfr
        new.nb_extra = int(nb_extra)
        new.meta=meta

        # Relaxation for the Ag configuration.`
        new.gs_relax_task = new.register_relax_task(gs_scf_inp.new_with_vars(relax_kwargs_gs))

        # Internal counter used in on_all_ok to drive the differ steps of the calculations.
        new.iteration_step = 0

        # JSON-compatible dictionary storing the most important results of the Work. (relaxations and 4 pts results
        # are separated)
        # Results will be stored in the `lumi_4pts.json` file in the outdata directory of the Work
        # so that one can easily implement additional post-processing tools.
        new.json_data = {}

        return new

    def on_all_ok(self):
        """
        This method is called when all the works in the flow have reached S_OK.

        This is the section in which we implement most of the workflow logic at runtime.
        since we need to generate input files with relaxed structures.
        """
            # Get Ag relaxed structure
        #with self.gs_relax_task.open_gsr() as gsr:
           # ag_relaxed_structure = gsr.structure
        with abiopen(self.gs_relax_task.output_file.path) as relax_gs_abo:
           ag_relaxed_structure = relax_gs_abo.final_structure

        if self.iteration_step == 0:
            print("in iteration step 0")
            self.iteration_step += 1

            # Relax geometry with excited configuration starting from Ag*.
            relax_ex_inp=self.ex_scf_inp.new_with_vars(self.relax_kwargs_ex)
            relax_ex_inp_2 = relax_ex_inp.new_with_structure(ag_relaxed_structure)
            self.ex_relax_task = self.register_relax_task(relax_ex_inp_2)#,deps={self.gs_relax_task: "DEN"})

            # if only the two relaxation, go to results writing step directly
            if self.four_points==False:
                self.iteration_step = 2

            return self.postpone_on_all_ok()

        elif self.iteration_step == 1:
            print("in iteration step 1")
            self.iteration_step += 1

            ##### SCF task for each configuration and optionnaly their electronic band structures ####

            # Build GS SCF input for the Ag configuration:
            # use same structure as Ag with ground occupation factors.
            ag_scf_inp = self.gs_scf_inp.new_with_structure(ag_relaxed_structure)
            self.ag_scf_task = self.register_scf_task((ag_scf_inp),deps={self.gs_relax_task: "DEN"})

            # Build GS SCF input for the Ag* configuration:
            # use same structure as Ag but with excited occupation factors.
            agstar_scf_inp = self.ex_scf_inp.new_with_structure(ag_relaxed_structure)
            self.agstar_scf_task = self.register_scf_task((agstar_scf_inp),deps={self.gs_relax_task: "DEN"})

            # Get Aestar relaxed structure.
            #with self.ex_relax_task.open_gsr() as gsr:
            #    aestar_relaxed_structure = gsr.structure
            with abiopen(self.ex_relax_task.output_file.path) as relax_ex_abo:
                aestar_relaxed_structure = relax_ex_abo.final_structure

            # Build ex SCF input for the Aestar configuration:
            # use same structure as Aestar with excited occupation factors.
            aestar_scf_inp = self.ex_scf_inp.new_with_structure(aestar_relaxed_structure)
            self.aestar_scf_task = self.register_scf_task((aestar_scf_inp),deps={self.ex_relax_task: "DEN"})

            # Build GS SCF task for the Ae configuration:
            # use same structure as Aestar but with ground occupation factors.
            ae_scf_inp = self.gs_scf_inp.new_with_structure(aestar_relaxed_structure)
            self.ae_scf_task = self.register_scf_task((ae_scf_inp),deps={self.ex_relax_task: "DEN"})


            if self.ndivsm != 0:
                # Compute band structure for Ag configuration.
                self.ag_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                             tolwfr=self.tolwfr, nb_extra=self.nb_extra)

                # Compute band structure for Agstar configuration.
                self.agstar_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                             tolwfr=self.tolwfr, nb_extra=self.nb_extra)

                # Compute band structure for aestar configuration.
                self.aestar_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                             tolwfr=self.tolwfr, nb_extra=self.nb_extra)

                # Compute band structure for Ae configuration.
                self.ae_scf_task.add_ebands_task_to_work(self, ndivsm=self.ndivsm,
                                                             tolwfr=self.tolwfr, nb_extra=self.nb_extra)

            return self.postpone_on_all_ok()

        elif self.iteration_step == 2:

            print("in iteration step 2")
            self.iteration_step += 1
            ##### Writing the results in json files #####

            self.json_data["meta"] = self.meta

            #with self.gs_relax_task.open_gsr() as gsr:
            self.json_data["gs_relax_filepath"]=self.gs_relax_task.gsr_path

            #with self.ex_relax_task.open_gsr() as gsr:
            self.json_data["ex_relax_filepath"]=self.ex_relax_task.gsr_path


            if self.four_points == True:
                # Get Ag total energy.
                #with self.ag_scf_task.open_gsr() as gsr:
                self.json_data["Ag_gsr_filepath"] = self.ag_scf_task.gsr_path

                # Get Agstar total energy.
                #with self.ex_relax_task.open_gsr() as gsr:
                self.json_data["Agstar_gsr_filepath"] = self.agstar_scf_task.gsr_path

                # Get Aestar total energy.
                #with self.aestar_scf_task.open_gsr() as gsr:
                self.json_data["Aestar_gsr_filepath"] = self.aestar_scf_task.gsr_path

                # Get Aestar total energy.
                #with self.ae_scf_task.open_gsr() as gsr:
                self.json_data["Ae_gsr_filepath"] = self.ae_scf_task.gsr_path

            # Write json file in the outdir of the work
            self.write_json_in_outdir("lumi.json", self.json_data)

            # Build deltascf results 
            delta_scf = DeltaSCF.from_four_points_file([self.ag_scf_task.gsr_path,
                                                        self.agstar_scf_task.gsr_path,
                                                        self.aestar_scf_task.gsr_path,
                                                        self.ae_scf_task.gsr_path])

            # Create dict with all post-processed results
            d = delta_scf.get_dict_results()

            # save d in json format.
            self.write_json_in_outdir("Delta_SCF.json", d)

            return super().on_all_ok()

class LumiWork_relaxations(Work):
    """
    This Work implements the ground and excited state relaxations only.

    The relaxations run simultaneously. No task creation at run-time. 
    """

    @classmethod
    def from_scf_inputs(cls, gs_scf_inp, ex_scf_inp, relax_kwargs_gs, relax_kwargs_ex, meta=None, manager=None):
        """
        Args:
            gs_scf_inp: |AbinitInput| representing a GS SCF run for the ground-state.
            ex_scf_inp: |AbinitInput| representing a GS SCF run for the excited-state.
            relax_kwargs_gs: Dictonary with input variables to be added to gs_scf_inp
                when generating input files for ground state structural relaxations.
            relax_kwargs_ex: Dictonary with input variables to be added to ex_scf_inp
                when generating input files for excited state structural relaxations.
            meta : dict corresponding to the metadata of a lumiwork (supercell size, dopant type,...)
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
        """
        new = cls(manager=manager)

        # Templates for GS SCF calculations.
        new.gs_scf_inp = gs_scf_inp
        new.ex_scf_inp = ex_scf_inp

        # Save paramenters for the generation of input files at runtime.
        new.relax_kwargs_gs = relax_kwargs_gs
        new.relax_kwargs_ex = relax_kwargs_ex

        new.meta=meta

        # Relaxation for the Ag configuration.
        new.gs_relax_task = new.register_relax_task(gs_scf_inp.new_with_vars(relax_kwargs_gs))
        new.ex_relax_task = new.register_relax_task(ex_scf_inp.new_with_vars(relax_kwargs_ex))


        new.json_data = {}

        return new

    def on_all_ok(self):
        """
        This method is called when all the works in the flow have reached S_OK.

        """
        self.json_data["meta"] = self.meta

        # with self.gs_relax_task.open_gsr() as gsr:
        self.json_data["gs_relax_filepath"] = self.gs_relax_task.gsr_path

        # with self.ex_relax_task.open_gsr() as gsr:
        self.json_data["ex_relax_filepath"] = self.ex_relax_task.gsr_path

        # Write json file in the outdir of the work
        self.write_json_in_outdir("lumi_relaxations.json", self.json_data)

        return super().on_all_ok()




class LumiWorkFromRelax(Work):
    """
    Same as LumiWork, without the relaxations. Typically used after a LumiWork_relaxations work.
    The two relaxed structures (in ground and excited state) are given as input. No creation at run-time

    """

    @classmethod
    def from_scf_inputs(cls, gs_scf_inp, ex_scf_inp, gs_structure, ex_structure, ndivsm=0, nb_extra=10,
                        tolwfr=1e-12, meta=None ,manager=None):
        """
        Args:
            gs_scf_inp: |AbinitInput| representing a GS SCF run for the ground-state.
            exc_scf_inp: |AbinitInput| representing a GS SCF run for the excited-state.
            gs_structure: object representing the relaxed ground state structure
            ex_structure: object representing the excited ground state structure
            ndivsm: Activates the computation of band structure if different from zero.
                if > 0, it's the number of divisions for the smallest segment of the path (Abinit variable).
                if < 0, it's interpreted as the pymatgen `line_density` parameter in which the number of points
                in the segment is proportional to its length. Typical value: -20.
                This option is the recommended one if the k-path contains two high symmetry k-points that are very close
                as ndivsm > 0 may produce a very large number of wavevectors.
            nb_extra: Number of extra bands added to the input nband when computing band structures (ndivsm != 0).
            tolwfr: Tolerance of the residuals used for the NSCF band structure calculations.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
        """
        new = cls(manager=manager)

        new.meta=meta
        new.gs_structure=gs_structure
        new.ex_structure=ex_structure

        # Templates for GS SCF calculations.
        new.gs_scf_inp = gs_scf_inp
        new.ex_scf_inp = ex_scf_inp

        new.ndivsm = int(ndivsm)
        new.tolwfr = tolwfr
        new.nb_extra = int(nb_extra)

        # JSON-compatible dictionary storing the most important results of the Work. (relaxations and 4 pts results
        # are separated)
        # Results will be stored in the `lumi_4pts.json` or 'lumi_relax.json file in the outdata directory of the Work
        # so that one can easily implement additional post-processing tools.
        new.json_data = {}

        # Build GS SCF input for the Ag configuration:
        # use same structure as Ag but with ground occupation factors
        ag_scf_inp = new.gs_scf_inp.new_with_structure(gs_structure)
        new.ag_scf_task = new.register_scf_task(ag_scf_inp)

        # Build GS SCF input for the Ag* configuration:
        # use same structure as Ag but with excited occupation factors.
        agstar_scf_inp = new.ex_scf_inp.new_with_structure(gs_structure)
        new.agstar_scf_task = new.register_scf_task(agstar_scf_inp)

        # Build ex SCF input for the Aestar configuration:
        aestar_scf_inp = new.ex_scf_inp.new_with_structure(ex_structure)
        new.aestar_scf_task = new.register_scf_task(aestar_scf_inp)


        # Build GS SCF input for the Ag* configuration:
        # use same structure as Ag but with excited occupation factors.
        ae_scf_inp = new.gs_scf_inp.new_with_structure(ex_structure)
        new.ae_scf_task = new.register_scf_task(ae_scf_inp)

        if new.ndivsm != 0:
            # Compute band structure for Ag configuration.
            new.ag_scf_task.add_ebands_task_to_work(new, ndivsm=new.ndivsm,
                                                         tolwfr=new.tolwfr, nb_extra=new.nb_extra)
            # Compute band structure for Ag* configuration.
            new.agstar_scf_task.add_ebands_task_to_work(new, ndivsm=new.ndivsm,
                                                             tolwfr=new.tolwfr, nb_extra=new.nb_extra)
            # Compute band structure for aestar configuration.
            new.aestar_scf_task.add_ebands_task_to_work(new, ndivsm=new.ndivsm,
                                                         tolwfr=new.tolwfr, nb_extra=new.nb_extra)
            # Compute band structure for Ae configuration.
            new.ae_scf_task.add_ebands_task_to_work(new, ndivsm=new.ndivsm,
                                                             tolwfr=new.tolwfr, nb_extra=new.nb_extra)

        return new

    def on_all_ok(self):

        self.json_data["meta"] = self.meta


        # Get Ag total energy.
        # with self.ag_scf_task.open_gsr() as gsr:
        self.json_data["Ag_gsr_filepath"] = self.ag_scf_task.gsr_path

        # Get Agstar total energy.
        # with self.ex_relax_task.open_gsr() as gsr:
        self.json_data["Agstar_gsr_filepath"] = self.agstar_scf_task.gsr_path

        # Get Aestar total energy.
        # with self.aestar_scf_task.open_gsr() as gsr:
        self.json_data["Aestar_gsr_filepath"] = self.aestar_scf_task.gsr_path

        # Get Aestar total energy.
        # with self.ae_scf_task.open_gsr() as gsr:
        self.json_data["Ae_gsr_filepath"] = self.ae_scf_task.gsr_path

        # Write json file in the outdir of the work
        self.write_json_in_outdir("lumi.json", self.json_data)

        # Build deltascf results 
        delta_scf = DeltaSCF.from_four_points_file([self.ag_scf_task.gsr_path,
                                                        self.agstar_scf_task.gsr_path,
                                                        self.aestar_scf_task.gsr_path,
                                                        self.ae_scf_task.gsr_path])

        # Create dict with all post-processed results
        d = delta_scf.get_dict_results()
        # save d in json format.
        self.write_json_in_outdir("Delta_SCF.json", d)

        return super().on_all_ok()

