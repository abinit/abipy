# coding: utf-8
"""
Works and Flows for GW calculations with Abinit.
"""
import os
import numpy as np

from . import wrappers
from .nodes import Node
from .tasks import NscfTask
from .works import Work


class ScreeningWork(Work):
    """
    This work parallelizes the calculation of the q-points of the screening.
    It also provides the callback `on_all_ok` that calls the mrgscr tool
    to merge all the partial SCR files produced.

    The final SCR file is made available in the `outdata` directory of the work.
    To use this SCR file as dependency of other tasks, use the syntax:

        deps={scr_work: "SCR"}

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ScreeningWork
    """

    @classmethod
    def from_nscf_task(cls, nscf_task, scr_input, manager=None):
        """
        Construct a `ScreeningWork` from a |NscfTask| object.

        Args:
            nscf_task: |NscfTask| object.
            scr_input: |AbinitInput| object representing a SCREENING calculation.
            manager: |TaskManager| object.
        """
        if not isinstance(nscf_task, NscfTask):
            raise TypeError("task `%s` does not inherit from NscfTask" % nscf_task)

        new = cls(manager=manager)

        # Get list of q-points for the dielectric matrix
        scr_ibz = scr_input.abiget_scr_ibz()

        # Now we can register the task for the different q-points
        for qpoint in scr_ibz.points:
            new.register_scr_task(scr_input.new_with_vars(nqptdm=1, qptdm=qpoint),
                                  deps={nscf_task: "WFK"})

        return new

    @classmethod
    def from_wfk_filepath(cls, wfk_filepath, scr_input, manager=None):
        """
        Construct a `ScreeningWork` from a WFK filepath and a screening input.

        WARNING: the parameters reported in `scr_input` e.g. k-mesh, nband, ecut etc
        must be consisten with those used to generate the WFK file.
        No consistencty check is done at this level. You have been warned!

        Args:
            wfk_filepath: Path to the WFK file.
            scr_input: |AbinitInput| object representing a SCREENING calculation.
            manager: |TaskManager| object.
        """
        new = cls(manager=manager)
        wkf_node = Node.as_node(wfk_filepath)

        # Get list of q-points for the dielectric matrix
        scr_ibz = scr_input.abiget_scr_ibz()

        # Now we can register the task for the different q-points
        for qpoint in scr_ibz.points:
            new.register_scr_task(scr_input.new_with_vars(nqptdm=1, qptdm=qpoint),
                                  deps={wfk_node: "WFK"})

        return new

    def merge_scr_files(self, remove_scrfiles=True, verbose=0):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Work`.
        If remove_scrfiles is True, the partial SCR files are removed after the merge.
        """
        scr_files = list(filter(None, [task.outdir.has_abiext("SCR") for task in self]))

        self.history.info("Will call mrgscr to merge %s SCR files:\n" % len(scr_files))
        assert len(scr_files) == len(self)

        mrgscr = wrappers.Mrgscr(manager=self[0].manager, verbose=verbose)
        final_scr = mrgscr.merge_qpoints(self.outdir.path, scr_files, out_prefix="out")

        if remove_scrfiles:
            for scr_file in scr_files:
                try:
                    os.remove(scr_file)
                except IOError:
                    pass

        return final_scr

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Work`.
        """
        final_scr = self.merge_scr_files()
        return self.Results(node=self, returncode=0, message="mrgscr done", final_scr=final_scr)
