# coding: utf-8
"""
Flow for the computation of Grunesein parameters with Abinit DFPT and finite differences.

WARNING: This code must be tested more carefully.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.core.structure import Structure

from .flows import Flow
from .works import Work


class GruneseinenFlow(Flow):
    """
    This work compute the Grüneisen parameters with the DFPT part of Abinit and
    finite differences in anaddb.

    It is necessary to run three phonon calculations with DFPT.
    One is calculated at the equilibrium volume and the remaining two are calculated
    at the slightly larger volume and smaller volume than the equilibrium volume.
    The unitcells at these volumes have to be fully relaxed under the constraint of each volume.
    """
    @classmethod
    def from_gs_input(cls, gsinp, voldelta, workdir=None, manager=None):
        """
        Build the work from an :class:`AbinitInput` object representing a GS calculations.

	Args:
	    gsinp:
		:class:`AbinitInput` object representing a GS calculation in the initial unit cell.
	    voldelta:
                Absolute increment for unit cell volume. The three volumes are:
                     [v0 - voldelta, v0, v0 + voldelta] where v0 is taken from gsinp.structure.

	Return:
        """
        new = cls(workdir=workdir, manager=manager)

	# Save arguments that will be used to call phonopy for creating
        # the supercells with the displacements once the three volumes have been relaxed.
        #new.scdims = np.array(scdims)
        #if new.scdims.shape != (3,):
        #    raise ValueError("Expecting 3 int in scdims but got %s" % str(new.scdims))
        #new.phonopy_kwargs = phonopy_kwargs if phonopy_kwargs is not None else {}
        #new.displ_kwargs = displ_kwargs if displ_kwargs is not None else {}

        # Build three tasks for structural optimization at constant volume.
        v0 = gsinp.structure.volume
        if voldelta <= 0:
            raise ValueError("voldelta must be > 0 but got %s" % voldelta)
        volumes = [v0 - voldelta, v0, v0  + voldelta]
        if any(v <= 0 for v in volumes):
            raise ValueError("volumes must be > 0 but got %s" % str(volumes))

        for vol in volumes:
            # Build new structure
            new_lattice = gsinp.structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, gsinp.structure.species, gsinp.structure.frac_coords)
            new_input = gsinp.new_with_structure(new_structure)
            # Set variables for structural optimization at constant volume.
            new_input.pop_tolerances()
            new_input.set_vars(optcell=3, ionmov=3, tolvrs=1e-10, toldff=1.e-6)
            new_input.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)
            new.register_relax_task(new_input)

        return new


class GruneseinenWork(Work):

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK.
        """
        self.add_phonopy_works_and_build()
        return super(GruneisenWork, self).on_all_ok()

    def add_phonopy_works_and_build(self):
        """
        Get relaxed structures from the tasks, build Phonopy works with supercells
        constructed from the new structures, add them to the flow and build new directories.
        """
        for i, task in enumerate(self):
            relaxed_structure = task.get_final_structure()
            gsinp = task.input.new_with_structure(relaxed_structure)

            work = PhonopyWork.from_gs_input(gsinp, self.scdims,
                                             phonopy_kwargs=self.phonopy_kwargs,
                                             displ_kwargs=self.displ_kwargs)

            self.flow.register_work(work)
            # Tell the work to copy the results to e.g. `flow/outdir/w0/minus`
            dst = os.path.join(self.pos_str, {0: "minus", 1: "orig", 2: "plus"}[i])
            work.cpdata2dst = self.flow.outdir.path_in(dst)

        # Allocate new works and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()


#class RelaxAndAddPhGammaWork(RelaxWork):
#    """
#    Work for structural relaxations. The first task relaxes the atomic position
#    while keeping the unit cell parameters fixed. The second task uses the final
#    structure to perform a structural relaxation in which both the atomic positions
#    and the lattice parameters are optimized.
#    """
#
#    def on_all_ok(self):
#        """
#        Here I extend the implementation of super in order to create a new workflow
#        for phonons with the optimized structural parameters.
#        """
#        results = super(RelaxAndAddPhGammaWork, self).on_all_ok()
#
#        # Get the relaxed structure.
#        relax_task = self[1]
#        final_structure = relax_task.get_final_structure()
#
#        # Use new structure in GS + DFPT runs and change some values.
#        scf_input = relax_task.input.deepcopy()
#        scf_input.set_structure(final_structure)
#
#        # Remove input variables that can enter into conflict with DFPT.
#        scf_input.pop_tolerances()
#        scf_input.pop_par_vars()
#        scf_input.pop_irdvars()
#        scf_input.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])
#        scf_input.set_vars(tolwfr=1e-20, nstep=80, nbdbuf=4)
#        #nval = scf_input.num_valence_electrons
#
#        # Build GS work and Phonon Work
#        work = PhononDojoWork.from_scf_input(scf_input, self.dojo_qpt)
#        for task in work[1:]:
#            task.set_vars(prtwf=-1)
#
#        # Monkey patch work
#        work.dojo_kppa = self.dojo_kppa
#        work.dojo_qpt = self.dojo_qpt
#        work.ecut = self.ecut
#        work.dojo_pawecutdg = self.dojo_pawecutdg
#        work.dojo_include_soc = self.dojo_include_soc
#        work._dojo_trial = "phgamma" if not self.dojo_include_soc else "phgamma_soc"
#        work._pseudo = self.dojo_pseudo
#
#        self.flow.register_work(work)
#        # Allocate new tasks and update the pickle database.
#        self.flow.allocate()
#        self.flow.build_and_pickle_dump()
#
#        return results
