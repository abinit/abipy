#!/usr/bin/env python
from __future__ import division, print_function

import os
import time
import collections
import yaml
import numpy as np

from abipy.abilab import AbinitFlow, Mrgscr
from pymatgen.io.abinitio.tasks import ScfTask, PhononTask
from pymatgen.io.abinitio.workflows import Workflow, BandStructureWorkflow, IterativeWorkflow, PhononWorkflow
from pymatgen.io.abinitio.abiinspect import yaml_read_kpoints, yaml_read_irred_perts

import logging
logger = logging.getLogger(__name__)


class QptdmWorkflow(Workflow):
    """
    This workflow parallelizes the calculation of the q-points of the screening. 
    It also provides the callback `on_all_ok` that calls mrgscr to merge 
    all the partial screening files produced.
    """
    def fill(self, wfk_file, scr_input):
        """
        Create the SCR tasks and register them in self.

        Args:
            wfk_file:
                Path to the ABINIT WFK file to use for the computation of the screening.
            scr_input:
                Input for the screening calculation.

        Return:
            self
        """
        assert len(self) == 0
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)

        # Build a temporary workflow in the tmpdir that will use a shell manager 
        # to run ABINIT in order to get the list of q-points for the screening.
        shell_manager = self.manager.to_shell_manager(mpi_ncpus=1)

        w = Workflow(workdir=self.tmpdir.path_join("_qptdm_run"), manager=shell_manager)

        fake_input = scr_input.deepcopy()
        fake_task = w.register(fake_input)
        w.allocate()
        w.build()

        # Create the symbolic link and add the magic value 
        # nqpdm = -1 to the input to get the list of q-points.
        fake_task.inlink_file(wfk_file)
        fake_task.strategy.add_extra_abivars({"nqptdm": -1})
        w.start_and_wait()

        # Parse the section with the q-points
        try:
            qpoints = yaml_read_kpoints(fake_task.log_file.path, doc_tag="!Qptdms")
            #print(qpoints)
        except:
            raise
        finally:
            w.rmtree()

        # Now we can register the task for the different q-points 
        for qpoint in qpoints:
            qptdm_input = scr_input.deepcopy()
            qptdm_input.set_variables(
                nqptdm=1,
                qptdm=qpoint
            )
            self.register(qptdm_input, manager=self.manager)

        return self.allocate()

    def merge_scrfiles(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        """
        scr_files = filter(None, [task.outdir.has_abiext("SCR") for task in self])

        logger.debug("will call mrgscr to merge %s:\n" % str(scr_files))
        assert len(scr_files) == len(self)

        mrgscr = Mrgscr(verbose=1)
        mrgscr.set_mpi_runner("mpirun")
        final_scr = mrgscr.merge_qpoints(scr_files, out_prefix="out", cwd=self.outdir.path)

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        """
        final_scr = self.merge_scrfiles()

        results = dict(
            returncode=0,
            message="mrgscr done",
            final_scr=final_scr)

        return results


def cbk_qptdm_workflow(flow, work, cbk_data):
    scr_input = cbk_data["input"]
    # Use the WFK file produced by the second 
    # Task in the first Workflow (NSCF step).
    nscf_task = flow[0][1]
    wfk_file = nscf_task.outdir.has_abiext("WFK")

    work.set_manager(flow.manager)
    work.fill(wfk_file, scr_input)
    #work.add_deps(cbk.deps)
    work.connect_signals()
    work.build()

    return work


def g0w0_flow_with_qptdm(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):
    """
    Build an `AbinitFlow` for one-shot G0W0 calculations.
    The computation of the q-points for the screening is parallelized with qptdm
    i.e. we run independent calculation for each q-point and then we merge 
    the final results.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).
        scr_input:
            Input for the SCR run.
        sigma_input:
            Input for the SIGMA run.

    Returns:
        `AbinitFlow`
    """                                                      
    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)

    # Register the first workflow (GS + NSCF calculation)
    bands_work = flow.register_work(BandStructureWorkflow(scf_input, nscf_input))

    assert not bands_work.scf_task.depends_on(bands_work.scf_task)
    assert bands_work.nscf_task.depends_on(bands_work.scf_task)

    # Register the callback that will be executed the workflow for the SCR with qptdm.
    scr_work = flow.register_cbk(cbk=cbk_qptdm_workflow, cbk_data={"input": scr_input},
                                 deps={bands_work.nscf_task: "WFK"}, work_class=QptdmWorkflow
                                )
                             
    assert scr_work.depends_on(bands_work.nscf_task)
    assert not scr_work.depends_on(bands_work.scf_task)

    # The last workflow contains a single SIGMA task that will use 
    # the data produced in the previous two workflows.
    sigma_task = flow.register_task(sigma_input, deps={bands_work.nscf_task: "WFK", scr_work: "SCR"})

    flow.allocate()
    assert sigma_task.depends_on(bands_work.nscf_task)
    assert not sigma_task.depends_on(bands_work.scf_task)
    assert sigma_task.depends_on(scr_work)

    flow.show_dependencies()
    #print("sigma_work.deps", sigma_work.deps)
    print("sigma_task.deps", sigma_task.deps)

    return flow
