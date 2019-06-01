# coding: utf-8
"""
Flows for electron-phonon calculations (high-level interface)
"""
from __future__ import unicode_literals, division, print_function

from .works import Work, PhononWork, PhononWfkqWork
from .flows import Flow


class GkqPathFlow(Flow):
    """
    Flow to compute the gkq matrix elements for a list of q-points (usually a q-path).
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ngqpt, qpath_list, with_becs=False, manager=None):
        """
        Build the flow from an input file representing a GS calculation.
        """
        flow = cls(workdir=workdir, manager=manager)

        # First work with GS run.
        scf_task = flow.register_scf_task(scf_input)[0]

        # Second work to compute phonons on nqgpt q-mesh.
        work_qmesh = PhononWork.from_scf_task(scf_task, qpoints=ngqpt, is_ngqpt=True, with_becs=with_becs)
        flow.register_work(work_qmesh)

        # Compute WFK/WFKQ and phonons for qpt in qpath_list. Don't include becs because already computed previously.
        work_qpath = PhononWfkqWork.from_scf_task(scf_task, qpath_list, ph_tolerance=None, tolwfr=1.0e-22, nband=None,
                      with_becs=False, ddk_tolerance=None, shiftq=(0, 0, 0), is_ngqpt=False, remove_wfkq=False,
                      manager=manager)
        flow.register_work(work_qpath)

        def make_eph_input(scf_inp, ngqpt, qpt):
            """
            Build input file to compute GKQ.nc file from GS SCF input.
            The calculation requires GS wavefunctions WFK, WFQ a DDB and a DVDB file
            """
            return scf_inp.new_with_vars(
                optdriver=7,
                eph_task=-2,
                nqpt=1,
                qpt=qpt,
                ddb_ngqpt=ngqpt,  # q-mesh associated to the DDB file.
                prtphdos=0,
            )

        # Compute matrix elements fully ab-initio for qpt in qpath_list.
        eph_work = Work()
        qseen = set()
        for task in work_qpath.phonon_tasks:
            qpt = tuple(task.input["qpt"])
            if qpt in qseen: continue
            qseen.add(qpt)
            t = eph_work.register_eph_task(make_eph_input(scf_input, ngqpt, qpt), deps=task.deps)
            t.add_deps({work_qmesh: "DDB", work_qpath: "DVDB"})
        flow.register_work(eph_work)

        # Here we buld another work to compute gkq with interpolated potentials.
        inteph_work = Work()
        qseen = set()
        for task in work_qpath.phonon_tasks:
            qpt = tuple(task.input["qpt"])
            if qpt in qseen: continue
            qseen.add(qpt)
            t = inteph_work.register_eph_task(make_eph_input(scf_input, ngqpt, qpt), deps=task.deps)
            t.add_deps({work_qmesh: ["DDB", "DVDB"]})
        flow.register_work(inteph_work)

        return flow
