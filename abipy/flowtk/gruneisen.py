# coding: utf-8
"""
Work for computing the Grüneisen parameters with finite differences of DFPT phonons.

WARNING: This code must be tested more carefully.
"""
import numpy as np

from abipy.core.structure import Structure
from .works import Work, PhononWork


class GruneisenWork(Work):
    """
    This work computes the Grüneisen parameters (derivative of frequencies wrt volume)
    using finite differences and the phonons obtained with the DFPT part of Abinit.
    The Anaddb input file needed to compute Grüneisen parameters will be generated
    in the outdata directory of the flow.

    It is necessary to run three DFPT phonon calculations.
    One is calculated at the equilibrium volume and the remaining two are calculated
    at the slightly larger volume and smaller volume than the equilibrium volume.
    The unitcells at these volumes have to be fully relaxed under the constraint of each volume.
    This Work automates the entire procedure starting from an input for GS calculations.
    """

    @classmethod
    def from_gs_input(cls, gs_inp, voldelta, ngqpt, tolerance=None, with_becs=False,
                      ddk_tolerance=None, workdir=None, manager=None):
        """
        Build the work from an |AbinitInput| representing a GS calculations.

        Args:
            gs_inp: |AbinitInput| representing a GS calculation in the initial unit cell.
            voldelta: Absolute increment for unit cell volume. The three volumes are:
                [v0 - voldelta, v0, v0 + voldelta] where v0 is taken from gs_inp.structure.
            ngqpt: three integers defining the q-mesh for phonon calculations.
            tolerance: dict {"varname": value} with the tolerance to be used in the phonon run.
                None to use AbiPy default.
            with_becs: Activate calculation of Electric field and Born effective charges.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run if with_becs.
                None to use AbiPy default.
        """
        new = cls(workdir=workdir, manager=manager)
        new.ngqpt = np.reshape(ngqpt, (3,))
        new.with_becs = with_becs
        new.ddk_tolerance = ddk_tolerance
        new.tolerance = tolerance

        if any(gs_inp["ngkpt"] % new.ngqpt != 0):
            raise ValueError("Kmesh and Qmesh must be commensurate.\nGot ngkpt: `%s`\nand ngqpt: `%s`" % (
                             str(gs_inp["ngkpt"]), str(new.ngqpt)))

        # Build three tasks for structural optimization at constant volume.
        v0 = gs_inp.structure.volume
        if voldelta <= 0:
            raise ValueError("voldelta must be > 0 but got %s" % voldelta)
        volumes = [v0 - voldelta, v0, v0 + voldelta]
        if any(v <= 0 for v in volumes):
            raise ValueError("volumes must be > 0 but got %s" % str(volumes))

        # Keep a copy of the GS input that will be used to generate the Phonon Works
        new.gs_inp = gs_inp.deepcopy()

        new.relax_tasks = []
        for vol in volumes:
            # Build new structure
            new_lattice = gs_inp.structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, gs_inp.structure.species, gs_inp.structure.frac_coords)
            new_input = gs_inp.new_with_structure(new_structure)
            # Set variables for structural optimization at constant volume.
            new_input.pop_tolerances()
            new_input.set_vars(optcell=3, ionmov=3, tolvrs=1e-10, toldff=1.e-6)
            new_input.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)
            t = new.register_relax_task(new_input)
            new.relax_tasks.append(t)

        return new

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks have reached status S_OK.
        """
        self.add_phonon_works_and_build()
        return super().on_all_ok()

    def add_phonon_works_and_build(self):
        """
        Get relaxed structures from the tasks, build Phonons works with new structures.
        Add works to the flow and build new directories.
        """
        ddb_paths, struct_middle, middle_idx = [], None, None
        relaxed_structs = []

        for i, task in enumerate(self.relax_tasks):
            relaxed_structure = task.get_final_structure()
            relaxed_structs.append(relaxed_structure)
            if i == len(self.relax_tasks) // 2:
                middle_idx, struct_middle = i, relaxed_structure

            # work to compute phonons with new structure.
            gsinp_vol = self.gs_inp.new_with_structure(relaxed_structure)
            work = PhononWork.from_scf_input(gsinp_vol, self.ngqpt, is_ngqpt=True, tolerance=self.tolerance,
                                             with_becs=self.with_becs, ddk_tolerance=self.ddk_tolerance)
            # Add it to the flow.
            self.flow.register_work(work)

            ddb_paths.append(work.outdir.path_in("out_DDB"))

        # Write Anaddb input file to compute Gruneisen parameters in flow.outdata.
        from abipy.abio.inputs import AnaddbInput
        anaddb_inp = AnaddbInput.phbands_and_dos(struct_middle, self.ngqpt, nqsmall=20, ndivsm=20,
                                                 chneut=1 if self.with_becs else 0,
                                                 dipdip=1 if self.with_becs else 0,
                                                 lo_to_splitting=self.with_becs,
                                                 comment="Anaddb input file for Grunesein parameters")
        # Add DDB files for Gruns
        anaddb_inp["gruns_nddbs"] = len(ddb_paths)
        anaddb_inp["gruns_ddbs"] = "\n" + "\n".join('"%s"' % p for p in ddb_paths)
        in_path = self.flow.outdir.path_in("anaddb_gruns.in")
        anaddb_inp.write(in_path)

        files_file = []
        app = files_file.append
        app(in_path)                                        # 1) Path of the input file
        app(self.flow.outdir.path_in("anaddb_gruns.out"))   # 2) Path of the output file
        app(ddb_paths[middle_idx])                          # 3) Input derivative database (not used if Gruns)
        for i in range(4):
            app("FOOBAR")

        with open(self.flow.outdir.path_in("anaddb_gruns.files"), "wt") as fh:
            fh.write("\n".join(files_file))

        #task = AbinitTask.temp_shell_task(anaddb_inp, workdir=work.outdir, manager=self.manager)
        #task.start_and_wait(autoparal=False)

        #with_ebands = False
        #if with_ebands:
        #    bands_work = Work(manager=self.manager)
        #    for i, structure in enumerate(relaxed_structs):
        #        nscf_inp = self.gs_inp.new_with_structure(structure)
        #        nscf_inp.pop_tolerances()
        #        nscf_inp.set_kpath(ndivsm, kptbounds=None, iscf=-2)
        #        nscf_inp["tolwfr"] = 1e-18
        #        work.register_nscf_task(nscf_inp, deps={self.relax_tasks[i]: "DEN"})

        #    # Add it to the flow.
        #    self.flow.register_work(work)

        # Allocate new works and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()
