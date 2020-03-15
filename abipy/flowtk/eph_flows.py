# coding: utf-8
"""
Flows for electron-phonon calculations (high-level interface)
"""
import numpy as np

from abipy.core.kpoints import kpath_from_bounds_and_ndivsm
from .works import Work, PhononWork, PhononWfkqWork
from .flows import Flow


class EphPotFlow(Flow):
    r"""
    This flow computes the e-ph scattering potentials on a q-mesh defined by ngqpt
    and a list of q-points (usually a q-path) specified by the user.
    The DFPT potentials on the q-mesh are merged in the DVDB located in the outdata
    of the second work while the DFPT potentials on the q-path are merged in the DVDB
    located in the outdata of the third work.
    These DVDB files are then passed to the EPH code to compute the average over the unit
    cell of the periodic part of the scattering potentials as a function of q.
    Results are stored in the V1QAVG.nc files of the outdata of the tasks in the fourth work.
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ngqpt, qbounds,
                       ndivsm=5, with_becs=True, ddk_tolerance=None, prepgkk=0, manager=None):
        """
        Build the flow from an input file representing a GS calculation.

        Args:
            workdir: Working directory.
            scf_input: Input for the GS SCF run.
            ngqpt: 3 integers defining the q-mesh.
            qbounds: List of boundaries defining the q-path used for the computation of the GKQ files.
                The q-path is automatically generated using `ndivsm` and the reciprocal-space metric.
                If `ndivsm` is 0, the code assumes that `qbounds` contains the full list of q-points
                and no pre-processing is performed.
            ndivsm: Number of points in the smallest segment of the path defined by `qbounds`.
                Use 0 to pass list of q-points.
            with_becs: Activate calculation of Electric field and Born effective charges.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run if `with_becs`.
            prepgkk: 1 to activate computation of all 3 * natom perts (debugging option).
            manager: |TaskManager| object.
        """
        flow = cls(workdir=workdir, manager=manager)

        # First work with GS run.
        scf_task = flow.register_scf_task(scf_input)[0]

        # Second work to compute phonons on the input nqgpt q-mesh.
        work_qmesh = PhononWork.from_scf_task(scf_task, qpoints=ngqpt, is_ngqpt=True,
                                              with_becs=with_becs, ddk_tolerance=ddk_tolerance)
        flow.register_work(work_qmesh)

        if ndivsm > 0:
            # Generate list of q-points from qbounds and ndivsm.
            qpath_list = kpath_from_bounds_and_ndivsm(qbounds, ndivsm, scf_input.structure)
        elif ndivsm == 0:
            # Use input list of q-points.
            qpath_list = np.reshape(qbounds, (-1, 3))
        else:
            raise ValueError("ndivsm cannot be negative. Received ndivsm: %s" % ndivsm)

        # Third Work: compute WFK/WFQ and phonons for qpt in qpath_list.
        # Don't include BECS because they have been already computed in the previous work.
        work_qpath = PhononWfkqWork.from_scf_task(
                       scf_task, qpath_list, ph_tolerance=None, tolwfr=1.0e-22, nband=None,
                       with_becs=False, ddk_tolerance=None, shiftq=(0, 0, 0), is_ngqpt=False, remove_wfkq=True,
                       prepgkk=prepgkk, manager=manager)

        flow.register_work(work_qpath)

        # Now we compute matrix elements fully ab-initio for each q-point.
        eph_work = Work()

        for eph_task in (-15, 15):
            eph_inp = scf_input.new_with_vars(
                optdriver=7,
                ddb_ngqpt=ngqpt,    # q-mesh associated to the DDB file.
                #dvdb_ngqpt=ngqpt,  # q-mesh associated to the DDVDB file.
                prtphdos=0,
                eph_task=eph_task
            )

            if eph_task == -15:
                # Use DVDB with ab-initio POTS along q-path to produce V1QAVG
                deps = {work_qmesh: "DDB", work_qpath: "DVDB"}
            elif eph_task == 15:
                # Use q-mesh to interpolate along the same q-path as above.
                deps = {work_qmesh: ["DDB", "DVDB"]}
                eph_inp.set_vars(ph_nqpath=len(qpath_list), ph_qpath=qpath_list)

            eph_work.register_eph_task(eph_inp, deps=deps)

        flow.register_work(eph_work)

        return flow


class GkqPathFlow(Flow):
    r"""
    This flow computes the gkq e-ph matrix elements <k+q|\Delta V_q|k> for a list of q-points (usually a q-path).
    The results are stored in the GKQ.nc files for the different q-points. These files can be used to analyze the behaviour
    of the e-ph matrix elements as a function of qpts with the the objects provided by the abipy.eph.gkq module.
    It is also possible to compute the e-ph matrix elements using the interpolated DFPT potentials
    if test_ft_interpolation is set to True.
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ngqpt, qbounds,
                       ndivsm=5, with_becs=True, ddk_tolerance=None,
                       test_ft_interpolation=False, prepgkk=0, manager=None):
        """
        Build the flow from an input file representing a GS calculation.

        Args:
            workdir: Working directory.
            scf_input: Input for the GS SCF run.
            ngqpt: 3 integers defining the q-mesh.
            qbounds: List of boundaries defining the q-path used for the computation of the GKQ files.
                The q-path is automatically generated using `ndivsm` and the reciprocal-space metric.
                If `ndivsm` is 0, the code assumes that `qbounds` contains the full list of q-points
                and no pre-processing is performed.
            ndivsm: Number of points in the smallest segment of the path defined by `qbounds`.
                Use 0 to pass list of q-points.
            with_becs: Activate calculation of Electric field and Born effective charges.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run if `with_becs`.
            test_ft_interpolation: True to add an extra Work in which the GKQ files are computed
                using the interpolated DFPT potentials and the q-mesh defined by `ngqpt`.
                The quality of the interpolation depends on the convergence of the BECS, epsinf and `ngqpt`.
            prepgkk: 1 to activate computation of all 3 * natom perts (debugging option).
            manager: |TaskManager| object.
        """
        flow = cls(workdir=workdir, manager=manager)

        # First work with GS run.
        scf_task = flow.register_scf_task(scf_input)[0]

        # Second work to compute phonons on the input nqgpt q-mesh.
        work_qmesh = PhononWork.from_scf_task(scf_task, qpoints=ngqpt, is_ngqpt=True,
                                              with_becs=with_becs, ddk_tolerance=ddk_tolerance)
        flow.register_work(work_qmesh)

        if ndivsm > 0:
            # Generate list of q-points from qbounds and ndivsm.
            qpath_list = kpath_from_bounds_and_ndivsm(qbounds, ndivsm, scf_input.structure)
        elif ndivsm == 0:
            # Use input list of q-points.
            qpath_list = np.reshape(qbounds, (-1, 3))
        else:
            raise ValueError("ndivsm cannot be negative. Received ndivsm: %s" % ndivsm)

        # Third Work. Compute WFK/WFQ and phonons for qpt in qpath_list.
        # Don't include BECS because they have been already computed in the previous work.
        work_qpath = PhononWfkqWork.from_scf_task(
                       scf_task, qpath_list, ph_tolerance=None, tolwfr=1.0e-22, nband=None,
                       with_becs=False, ddk_tolerance=None, shiftq=(0, 0, 0), is_ngqpt=False, remove_wfkq=False,
                       prepgkk=prepgkk, manager=manager)

        flow.register_work(work_qpath)

        def make_eph_input(scf_inp, ngqpt, qpt):
            """
            Build input file to compute GKQ.nc file from GS SCF input.
            The calculation requires GS wavefunctions WFK, WFQ, a DDB file and a DVDB file
            """
            return scf_inp.new_with_vars(
                optdriver=7,
                eph_task=-2,
                nqpt=1,
                qpt=qpt,
                ddb_ngqpt=ngqpt,  # q-mesh associated to the DDB file.
                prtphdos=0,
            )

        # Now we compute matrix elements fully ab-initio for each q-point.
        eph_work = Work()

        qseen = set()
        for task in work_qpath.phonon_tasks:
            qpt = tuple(task.input["qpt"])
            if qpt in qseen: continue
            qseen.add(qpt)
            t = eph_work.register_eph_task(make_eph_input(scf_input, ngqpt, qpt), deps=task.deps)
            t.add_deps({work_qmesh: "DDB", work_qpath: "DVDB"})

        flow.register_work(eph_work)

        # Here we build another work to compute the gkq matrix elements
        # with interpolated potentials along the q-path.
        # The potentials are interpolated using the input ngqpt q-mesh.
        if test_ft_interpolation:
            inteph_work = Work()
            qseen = set()
            for task in work_qpath.phonon_tasks:
                qpt = tuple(task.input["qpt"])
                if qpt in qseen: continue
                qseen.add(qpt)
                eph_inp = make_eph_input(scf_input, ngqpt, qpt)
                # Note eph_use_ftinterp 1 to force the interpolation of the DFPT potentials with eph_task -2.
                eph_inp["eph_use_ftinterp"] = 1
                t = inteph_work.register_eph_task(eph_inp, deps=task.deps)
                t.add_deps({work_qmesh: ["DDB", "DVDB"]})
            flow.register_work(inteph_work)

        return flow
