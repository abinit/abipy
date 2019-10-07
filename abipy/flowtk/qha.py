# coding: utf-8
"""
Tools and workflows for calculations within the quasi-harmonic approximation.

WARNING: This code is still under development.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.flowtk.works import Work, RelaxWork, PhononWork, MergeDdb
from abipy.flowtk.flows import Flow

import logging
logger = logging.getLogger(__name__)


class RelaxAndPhononWork(Work, MergeDdb):
    """
    This work perform a structural relaxation and use the relaxed structure to
    compute phonon frequencies with DFPT.
    """
    _relax_done, _phonons_done = False, False

    @classmethod
    def from_gsinp(cls, gsinp, ngqpt, optcell, ionmov):
        work = cls()
        work.initial_gsinp, work.ngqpt = gsinp, ngqpt
        # Create input for relaxation and register relaxation task.
        relaxinp = gsinp.deepcopy()
        relaxinp.pop_tolerances()
        relaxinp.set_vars(tolvrs=1e-10, toldff=1.e-6)
        relaxinp.set_vars(optcell=optcell, ionmov=ionmov)
        if optcell is not None and optcell != 0 :
            relaxinp.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)

        work.relax_task = work.register_relax_task(relaxinp)
        return work

    def add_phonon_tasks_and_build(self):
        """
        Called by on_all_ok when the structural relaxation is completed.
        Here the phonon tasks with the relaxed structure and constructed
        and added to the work.
        """
        # Get relaxed structure and build tasks for phonon calculations.
        relaxed_structure = self.relax_task.get_final_structure()
        scf_input = self.initial_gsinp.new_with_structure(relaxed_structure)
        scf_task = self.register_scf_task(scf_input)

        qpoints = scf_input.abiget_ibz(ngkpt=self.ngqpt, shiftk=[0, 0, 0]).points
        for qpt in qpoints:
            for ph_inp in scf_task.input.make_ph_inputs_qpoint(qpt, tolerance=None):
                self.register_phonon_task(ph_inp, deps={scf_task: "WFK"})

        # Allocate new works and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def on_all_ok(self):
        """
        This method is called when all the task reach S_OK.
        The first time we enter here, we construct the phonon tasks.
        The second call runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        if not self._relax_done:
            # Build phonon tasks.
            self._relax_done = True
            self.add_phonon_tasks_and_build()
            self.finalized = False

        elif not self._phonons_done:
            # Merge DDB files.
            self._phonons_done = True
            self.merge_ddb_files()
            self.finalized = True

        return super(RelaxAndPhononWork, self).on_all_ok()


class QhaFlow(Flow):

    @classmethod
    def from_gsinp(cls, workdir, gsinp, volumes, ngqpt, manager=None):
        """

        Args:
            workdir:
            gsinp:
            volumes:
            ngqpt:
            manager:
        """
        ngqpt = np.reshape(ngqpt, 3)
        flow = cls(workdir=workdir, manager=manager)

        # Construct len(volumes) works. Each work performs the structure relaxation
        # at fixed volume followed by DFPT calculation with the relaxed structure.
        for vol in volumes:
            # Build GS input file for new structure with rescaled volume.
            new_lattice = gsinp.structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, gsinp.structure.species, gsinp.structure.frac_coords)
            new_input = gsinp.new_with_structure(new_structure)

            # Register work.
            work = RelaxAndPhononWork.from_gsinp(new_input, ngqpt, optcell=3, ionmov=3)
            flow.register_work(work)

        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It transfer the most important results to flow/outdir and perform
        some basic post-processing of the results to facilitate further analysis.
        """
        #data = {}
        #for work in self:
        #    ddb_path = work.outdir.has_abiext("DDB")
        #    with work.relax_task.open_gsr() as gsr
        #        data.append(dict(
        #            volume=gsr.structure.volume,
        #            etotal=gsr.etotal,
        #            structure=gsr.structure,
        #            pressure_GPa=gsr.pressure,
        #            stress_tensor=gsr.stress_tensor,
        #            forces=gsr.forces,
        #            ddb_path=ddb_path,
        #            gsr_path=gsr.path,
        #        ))

        #with open(self.outdir.path_in("qha.sjon"), "wt") as fh:
        #    json.dump(fh, data)

        return super(QhaFlow, self).finalize()
