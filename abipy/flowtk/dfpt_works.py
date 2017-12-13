# coding: utf-8
"""Work subclasses related to DFTP."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from .works import Work


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
        Build NscfDdksWork from a scf_task.

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
