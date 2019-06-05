# coding: utf-8
"""
Flows for electron-phonon calculations (high-level interface)
"""
from __future__ import unicode_literals, division, print_function

from abipy.core.kpoints import kpath_from_bounds_and_ndivsm

from .works import Work, PhononWork, PhononWfkqWork
from .flows import Flow


class GkqPathFlow(Flow):
    """
    This is an "advanced" flow to compute the gkq e-ph matrix elements for a list of
    q-points (usually a q-path). This is an advanced flow used to analyze the behaviour
    of the matrix elements as function of q and/or compare with the interpolated results.
    The matrix elements are stored in the GKQ.nc files that can be analyzed with the 
    the objects provided by the abipy.eph.gkq.module.
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ngqpt, qbounds, ndivsm=5, with_becs=False, 
                       test_ft_interpolation=True, manager=None):
        """
        Build the flow from an input file representing a GS calculation.

            workdir: Path to the working directory.
            scf_input: Input for the SCF run
            ngqpt: 3 integers defining the q-mesh 
            qbounds: List of boundaries defining the q-path for the computation of the GKQ files.
                The q-path is automatically generated using ndivsm and the reciprocal-space metric. 
                If ndivsm is negative, the  routine assumes that qbounds contains the full list of q-points
                and no post-processing is performed.
            with_becs: Activate calculation of Electric field and Born effective charges.
            test_ft_interpolation: True to add an extra Work in which in the GKQ files are computed 
                using the interpolated DFPT potentials using the q-mesh defined by ngqpt.
            manager: |TaskManager| object.
        """
        flow = cls(workdir=workdir, manager=manager)

        # First work with GS run.
        scf_task = flow.register_scf_task(scf_input)[0]

        # Second work to compute phonons on nqgpt q-mesh.
        work_qmesh = PhononWork.from_scf_task(scf_task, qpoints=ngqpt, is_ngqpt=True, with_becs=with_becs)
        flow.register_work(work_qmesh)

        if ndivsm > 0
            # Generate list of q-points from qbounds and ndivsm.
            qpath_list = kpath_from_bounds_and_ndivsm(qbounds, ndivsm, scf_input.structure)
        else:
            # Use input list of q-pooints.
            qpath_list = np.reshape(qbounds, (-1, 3))

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

        # Compute matrix elements fully ab-initio for each q-point.
        eph_work = Work()
        qseen = set()
        for task in work_qpath.phonon_tasks:
            qpt = tuple(task.input["qpt"])
            if qpt in qseen: continue
            qseen.add(qpt)
            t = eph_work.register_eph_task(make_eph_input(scf_input, ngqpt, qpt), deps=task.deps)
            t.add_deps({work_qmesh: "DDB", work_qpath: "DVDB"})
        flow.register_work(eph_work)

        # Here we build another work to compute gkq with interpolated potentials along the q-path.
        # Note eph_use_ftinterp 1 to force the interpolation of the DFPT potentials.
        if test_ft_interpolation:
            inteph_work = Work()
            qseen = set()
            for task in work_qpath.phonon_tasks:
                qpt = tuple(task.input["qpt"])
                if qpt in qseen: continue
                qseen.add(qpt)
                eph_inp = make_eph_input(scf_input, ngqpt, qpt)
                eph_inp["eph_use_ftinterp"] = 1
                t = inteph_work.register_eph_task(eph_inp, deps=task.deps)
                t.add_deps({work_qmesh: ["DDB", "DVDB"]})
            flow.register_work(inteph_work)

        return flow
