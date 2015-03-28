# coding: utf-8
"""DDB File."""
from __future__ import print_function, division, unicode_literals

import os
import tempfile
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.io.abinitio.tasks import AnaddbTask, TaskManager
from abipy.core.mixins import TextFile, Has_Structure
from abipy.core.symmetries import SpaceGroup
from abipy.core.structure import Structure
from abipy.abio.inputs import AnaddbInput

import logging
logger = logging.getLogger(__name__)


class TaskException(Exception):
    """
    Exceptions raised when we try to execute :class:`AnaddbTask` in the :class:`DdbFile` methods

    A `TaskException` has a reference to the task and to the :class:`EventsReport` that contains
    the error messages of the run.
    """
    def __init__(self, *args, **kwargs):
        self.task, self.report = kwargs.pop("task"), kwargs.pop("report")
        super(TaskException, self).__init__(*args, **kwargs)

    def __str__(self):
        lines = ["\nworkdir = %s" % self.task.workdir]
        app = lines.append

        if self.report.errors: 
            app("Found %d errors" % len(self.report.errors))
            lines += [str(err) for err in self.report.errors]

        return "\n".join(lines)


class DdbFile(TextFile, Has_Structure):
    """
    This object provides an interface to the DDB file produced by ABINIT
    as well as methods to compute phonon band structures, phonon DOS, thermodinamical properties ...
    """
    @lazy_property
    def structure(self):
        structure = Structure.from_abivars(**self.header)
        # Add Spacegroup (needed in guessed_ngkpt)
        # FIXME: has_timerev is always True
        spgid, has_timerev, h = 0, True, self.header
        structure.set_spacegroup(SpaceGroup(spgid, h.symrel, h.tnons, h.symafm, has_timerev))
        return structure

    @lazy_property
    def header(self):
        """
        Dictionary with the values reported in the header section. 
        Use ddb.header.ecut to access its values
        """
        return self._parse_header()

    def _parse_header(self):
        """Parse the header sections. Returns :class:`AttrDict` dictionary."""
        #ixc         7
        #kpt  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        #     0.25000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        self.seek(0)
        keyvals = []
        for i, line in enumerate(self):
            line = line.strip()
            if not line: continue
            if "Version" in line:
                # +DDB, Version number    100401
                version = int(line.split()[-1])

            if line == "Description of the potentials (KB energies)":
                # Skip section with psps info.
                break

            # header starts here
            if i >= 6:
                # Python does not support exp format with D 
                line = line.replace("D+", "E+").replace("D-", "E-")
                tokens = line.split()
                key = None
                try:
                    float(tokens[0])
                    parse = float if "." in tokens[0] else int
                    keyvals[-1][1].extend(map(parse, tokens))
                except ValueError:
                    # We have a new key
                    key = tokens.pop(0)
                    parse = float if "." in tokens[0] else int
                    keyvals.append((key, map(parse, tokens)))

        h = AttrDict(version=version)
        for key, value in keyvals:
            if len(value) == 1: value = value[0]
            h[key] = value

        # Convert to array. Note that znucl is converted into integer
        # to avoid problems with pymatgen routines that expect integral Z
        # This of course will break any code for alchemical mixing.
        arrays = {
            "kpt": dict(shape=(h.nkpt, 3), dtype=np.double),
            "rprim": dict(shape=(3, 3), dtype=np.double),
            "symrel": dict(shape=(h.nsym, 3, 3), dtype=np.int),
            "tnons": dict(shape=(h.nsym, 3), dtype=np.double),
            "xred":  dict(shape=(h.natom, 3), dtype=np.double),
            "znucl": dict(shape=(-1,), dtype=np.int),
        }

        for k, ainfo in arrays.items():
            h[k] = np.reshape(np.array(h[k], dtype=ainfo["dtype"]), ainfo["shape"])

        # Transpose symrel because Abinit write matrices by colums.
        h.symrel = np.array([s.T for s in h.symrel])
        
        return h

    @lazy_property
    def qpoints(self):
        """`ndarray` with the list of q-points in reduced coordinates."""
        return self._read_qpoints()

    def _read_qpoints(self):
        """Read the list q-points from the DDB file. Returns `ndarray`"""
        # 2nd derivatives (non-stat.)  - # elements :      36
        # qpt  2.50000000E-01  0.00000000E+00  0.00000000E+00   1.0

        # Since there are multiple occurrences of qpt in the DDB file
        # we use seen to remove duplicates.
        self.seek(0)
        tokens, seen = [], set()

        for line in self:
            line = line.strip()
            if line.startswith("qpt") and line not in seen:
                seen.add(line)
                tokens.append(line.replace("qpt", ""))

        qpoints, weights = [], []
        for tok in tokens:
            nums = list(map(float, tok.split()))
            qpoints.append(nums[:3])
            weights.append(nums[3])

        return np.reshape(qpoints, (-1,3))

    @lazy_property
    def guessed_ngqpt(self):
        """
        This function tries to figure out the value of ngqpt from the list of 
        points reported in the DDB file.

        .. warning::
            
            The mesh may not be correct if the DDB file contains points belonging 
            to different meshes and/or the Q-mesh is shifted.
        """
        # Build the union of the stars of the q-points.
        all_qpoints = np.empty((len(self.qpoints) * len(self.structure.spacegroup), 3))
        count = 0
        for qpoint in self.qpoints:
            for op in self.structure.spacegroup:
                all_qpoints[count] = op.rotate_k(qpoint, wrap_tows=False)
                count += 1

        # Replace zeros with np.inf
        for q in all_qpoints: 
            q[q == 0] = np.inf

        # Compute the minimum of the fractional coordinates along the 3 directions and invert
        #print(all_qpoints)
        smalls = np.abs(all_qpoints).min(axis=0)
        smalls[smalls == 0] = 1
        ngqpt = np.rint(1 / smalls)
        ngqpt[ngqpt == 0] = 1
        #print("smalls: ", smalls, "ngqpt", ngqpt)

        return np.array(ngqpt, dtype=np.int)

    @property
    def params(self):
        """Dictionary with the parameters that are usually tested for convergence."""
        return {k: v for k, v in self.header.items() if k in ("nkpt", "nsppol", "ecut", "tsmear", "ixc")}

    def calc_phmodes_at_qpoint(self, qpoint=None, asr=2, chneut=1, dipdip=1, 
                               workdir=None, manager=None, verbose=0, ret_task=False):
        """
        Execute anaddb to compute phonon modes at the given q-point.

        Args:
            qpoint: Reduced coordinates of the qpoint where phonon modes are computed
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            workdir: Working directory. If None, a temporary directory is created.
            manager: :class:`TaskManager` object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information

        Return:
            :class:`PhononBands` object.
        """
        if qpoint is None:
            qpoint = self.qpoints[0] 
            if len(self.qpoints) != 1:
                raise ValueError("%s contains %s qpoints and the choice is ambiguous.\n" 
                                 "Please specify the qpoint in calc_phmodes_at_qpoint" % (self, len(self.qpoints)))

        inp = AnaddbInput.modes_at_qpoint(self.structure, qpoint, asr=asr, chneut=chneut, dipdip=dipdip)

        if manager is None: manager = TaskManager.from_user_config()
        if workdir is None: workdir = tempfile.mkdtemp()
        if verbose: 
            print("workdir:", workdir)
            print("ANADDB INPUT:\n", inp)

        task = AnaddbTask(inp, self.filepath, workdir=workdir, manager=manager.to_shell_manager(mpi_procs=1))

        if ret_task:
            return task

        # Run the task here
        task.start_and_wait(autoparal=False)
        report = task.get_event_report()
        if not report.run_completed:
            raise TaskException(task=task, report=report)

        with task.open_phbst() as ncfile:
            return ncfile.phbands

    #def calc_phbands(self, ngqpt=None, ndivsm=20, asr=2, chneut=1, dipdip=1, workdir=None, manager=None, verbose=0, **kwargs):
    #def calc_phdos(self, ngqpt=None, nqsmall=10, asr=2, chneut=1, dipdip=1, dos_method="tetra" workdir=None, manager=None, verbose=0, **kwargs):

    def calc_phbands_and_dos(self, ngqpt=None, ndivsm=20, nqsmall=10, asr=2, chneut=1, dipdip=1, dos_method="tetra",
                             workdir=None, manager=None, verbose=0, plot=True, ret_task=False):
        """
        Execute anaddb to compute the phonon band structure and the phonon DOS

        Args:
            ngqpt: Number of divisions for the q-mesh in the DDB file. Auto-detected if None (default)
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            workdir: Working directory. If None, a temporary directory is created.
            manager: :class:`TaskManager` object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information
        """
        if ngqpt is None: ngqpt = self.guessed_ngqpt

        inp = AnaddbInput.phbands_and_dos(
            self.structure, ngqpt=ngqpt, ndivsm=ndivsm, nqsmall=nqsmall, 
            q1shft=(0,0,0), qptbounds=None, asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method)

        if manager is None: manager = TaskManager.from_user_config()
        if workdir is None: workdir = tempfile.mkdtemp()
        if verbose: 
            print("workdir:", workdir)
            print("ANADDB INPUT:\n", inp)

        task = AnaddbTask(inp, self.filepath, workdir=workdir, manager=manager.to_shell_manager(mpi_procs=1))

        if ret_task:
            return task

        # Run the task here.
        task.start_and_wait(autoparal=False)

        report = task.get_event_report()
        if not report.run_completed:
            raise TaskException(task=task, report=report)

        with task.open_phbst() as phbst_ncfile, task.open_phdos() as phdos_ncfile:
            phbands, phdos = phbst_ncfile.phbands, phdos_ncfile.phdos
            if plot:
                phbands.plot_with_phdos(phdos, title="Phonon bands and DOS of %s" % self.structure.formula)

            return phbands, phdos

    #def calc_thermo(self, nqsmall, ngqpt=None, workdir=None, manager=None, verbose=0):
    #    """
    #    Execute anaddb to compute thermodinamical properties.

    #    Args:
    #        ngqpt: Number of divisions for the q-mesh in the DDB file. Auto-detected if it is None
    #        asr, chneut, dipdp: Anaddb input variable. See official documentation.
    #        workdir: Working directory. If None, a temporary directory is created.
    #        manager: :class:`TaskManager` object. If None, the object is initialized from the configuration file
    #        verbose: verbosity level. Set it to a value > 0 to get more information
    #        kwargs: Additional variables you may want to pass to Anaddb.
    #    """
    #    if ngqpt is None: ngqpt = self.guessed_ngqpt

    #    inp = AnaddbInput.thermo(self.structure, ngqpt, nqsmall, q1shft=(0, 0, 0), nchan=1250, nwchan=5, thmtol=0.5,
    #           ntemper=199, temperinc=5, tempermin=5., asr=2, chneut=1, dipdip=1, ngrids=10)

    #    if manager is None: manager = TaskManager.from_user_config()
    #    if workdir is None: workdir = tempfile.mkdtemp()
    #    if verbose: 
    #        print("workdir:", workdir)
    #        print("ANADDB INPUT:\n", inp)

    #    task = AnaddbTask(inp, self.filepath, workdir=workdir, manager=manager.to_shell_manager(mpi_procs=1))
    #    task.start_and_wait(autoparal=False)

    #    report = task.get_event_report()
    #    if not report.run_completed:
    #        raise TaskException(task=task, report=report)


#class DdbConverger(object):
#    def __init__(self, ddb, workdir=workdir, manager=None, num_cores=None):
#        self.ddb = DdbFile(ddb)
#        from monty.dev import get_ncpus
#        self.num_cores = get_num_cpus() if num_cores is None else num_cores
#        
#        self.manager = TaskManager.from_user_config() if manager is None else manager
#        self.workdir = tempfile.mkdtemp() if workdir is None else workdir
#
#    def __enter__(self):
#        return self
#
#    def __exit__(self, exc_type, exc_val, exc_tb):
#        """Activated at the end of the with statement. It automatically closes the file."""
#        self.ddb.close()
#
#    def converge_phdos(self, nqsmall_slice):
#        doses = []
#        for nqs in nsmall_range:
#           #self.ddb.calc_phbands_and_dos(ngqpt=None, ndivsm=20, nqsmall=10, workdir=None, manager=self.manager)
#           doses.append(phdos)
#
#        # Compare last three phonon DOSes.
#        # Be careful here because the DOS may be defined on different frequency meshes
#        last_mesh = doses[-1].mesh
#        converged = False
#         for dos in doses[:-1]:
#            dos = dos.spline_on_mesh(last_mesh)
#            diffs.append(dos - doses[-1])
#
#        if converged:
#            return collections.namedtuple("results", "phdos ngsmall plotter")
#                (phdos=phdos, nqsmall=nqsmall, plotter=plotter)
#        else:
#            raise self.Error("Cannot converge the DOS wrt nqsmall")
