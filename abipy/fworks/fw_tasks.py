# coding: utf-8
"""
Task classes for Fireworks.
"""
from __future__ import print_function, division, unicode_literals

try:
    from fireworks.core.firework import Firework, FireTaskBase, FWAction
    from fireworks.utilities.fw_utilities import explicit_serialize
    from fireworks.utilities.fw_serializers import serialize_fw
except ImportError:
    FireTaskBase, FWAction, Firework = 3 * [object]
    explicit_serialize = lambda x: x
    serialize_fw = lambda x: x

import os
import inspect
import subprocess
import logging
import collections
import time
from pymatgen.io.abinitio.utils import Directory, File
from pymatgen.io.abinitio import events
from pymatgen.io.abinitio.nodes import Status
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.serializers.json_coders import pmg_serialize


logger = logging.getLogger(__name__)


@explicit_serialize
class AbiFireTask(FireTaskBase):

    # List of `AbinitEvent` subclasses that are tested in the check_status method.
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = [
    ]

    S_INIT = Status.from_string("Initialized")
    S_LOCKED = Status.from_string("Locked")
    S_READY = Status.from_string("Ready")
    S_SUB = Status.from_string("Submitted")
    S_RUN = Status.from_string("Running")
    S_DONE = Status.from_string("Done")
    S_ABICRITICAL = Status.from_string("AbiCritical")
    S_QCRITICAL = Status.from_string("QCritical")
    S_UNCONVERGED = Status.from_string("Unconverged")
    S_ERROR = Status.from_string("Error")
    S_OK = Status.from_string("Completed")

    ALL_STATUS = [
        S_INIT,
        S_LOCKED,
        S_READY,
        S_SUB,
        S_RUN,
        S_DONE,
        S_ABICRITICAL,
        S_QCRITICAL,
        S_UNCONVERGED,
        S_ERROR,
        S_OK,
    ]

    def __init__(self, abiinput):
        """
        Basic __init__, subclasses are supposed to define the same input parameters, add their own and call super for
        the basic ones. The input parameter should be stored as attributes of the instance for serialization and
        for inspection.
        """
        self.abiinput = abiinput

    @serialize_fw
    def to_dict(self):
        d = {}
        for arg in inspect.getargspec(self.__init__).args:
            if arg != "self":
                val = self.__getattribute__(arg)
                if hasattr(val, "as_dict"):
                    val = val.as_dict()
                elif isinstance(val, (tuple, list)):
                    val = [v.as_dict() if hasattr(v, "as_dict") else v for v in val]
                d[arg] = val

        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        kwargs = {k: dec.process_decoded(v) for k, v in d.items()
                  if k in inspect.getargspec(cls.__init__).args}
        return cls(**kwargs)

    #from Task
    def set_workdir(self, workdir):
        """Set the working directory."""

        self.workdir = os.path.abspath(workdir)

        # Files required for the execution.
        self.input_file = File(os.path.join(self.workdir, "run.abi"))
        self.output_file = File(os.path.join(self.workdir, "run.abo"))
        self.files_file = File(os.path.join(self.workdir, "run.files"))
        self.log_file = File(os.path.join(self.workdir, "run.log"))
        self.stderr_file = File(os.path.join(self.workdir, "run.err"))

        # Directories with input|output|temporary data.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

    # from Task
    def build(self):
        """
        Creates the working directory and the input files of the :class:`Task`.
        It does not overwrite files if they already exist.
        """
        # Create dirs for input, output and tmp data.
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        # Write files file and input file.
        if not self.files_file.exists:
            self.files_file.write(self.filesfile_string)

        self.input_file.write(str(self.abiinput))

    #from Task
    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata", "out"), pj("tmpdata", "tmp"))
    del Prefix, pj

    #from AbintTask
    @property
    def filesfile_string(self):
        """String with the list of files and prefixes needed to execute ABINIT."""
        lines = []
        app = lines.append
        pj = os.path.join

        app(self.input_file.path)                 # Path to the input file
        app(self.output_file.path)                # Path to the output file
        app(pj(self.workdir, self.prefix.idata))  # Prefix for input data
        app(pj(self.workdir, self.prefix.odata))  # Prefix for output data
        app(pj(self.workdir, self.prefix.tdata))  # Prefix for temporary data

        # Paths to the pseudopotential files.
        # Note that here the pseudos **must** be sorted according to znucl.
        # Here we reorder the pseudos if the order is wrong.
        ord_pseudos = []
        znucl = self.abiinput.structure.to_abivars()["znucl"]

        for z in znucl:
            for p in self.abiinput.pseudos:
                if p.Z == z:
                    ord_pseudos.append(p)
                    break
            else:
                raise ValueError("Cannot find pseudo with znucl %s in pseudos:\n%s" % (z, self.pseudos))

        for pseudo in ord_pseudos:
            app(pseudo.path)

        return "\n".join(lines)

    def run_abinit(self, fw_spec):

        with open(self.files_file.path, 'r') as stdin, open(self.log_file.path, 'w') as stdout, \
            open(self.stderr_file.path, 'w') as stderr:

            p = subprocess.Popen(['mpirun', 'abinit'], stdin=stdin, stdout=stdout, stderr=stderr)

        (stdoutdata, stderrdata) = p.communicate()
        self.returncode = p.returncode

    def get_event_report(self):
        """
        Analyzes the main output file for possible Errors or Warnings.

        Returns:
            :class:`EventReport` instance or None if the main output file does not exist.
        """

        if not self.log_file.exists:
            return None

        parser = events.EventsParser()
        try:
            report = parser.parse(self.log_file.path)
            return report

        except parser.Error as exc:
            # Return a report with an error entry with info on the exception.
            logger.critical("%s: Exception while parsing ABINIT events:\n %s" % (self.log_file, str(exc)))
            return parser.report_exception(self.log_file.path, exc)

    def task_analysis(self, fw_spec):

        status, msg = self.check_final_status()

        if self.status != self.S_OK:
            raise AbinitRuntimeError(self)

        return FWAction(stored_data=dict(**self.report.as_dict()))

    def run_task(self, fw_spec):
        self.set_workdir(os.path.abspath('.'))
        self.build()
        self.run_abinit(fw_spec)
        return self.task_analysis(fw_spec)

    def set_status(self, status, msg=None):
        self.status = status
        return status, msg

    def check_final_status(self):
        """
        This function checks the status of the task by inspecting the output and the
        error files produced by the application. Based on abipy task checkstatus().
        """
        # 2) see if an error occured at starting the job
        # 3) see if there is output
        # 4) see if abinit reports problems
        # 5) see if err file exists and is empty
        # 9) the only way of landing here is if there is a output file but no err files...

        # 2) Check the returncode of the process (the process of submitting the job) first.
        if self.returncode != 0:
            # The job was not submitted properly
            return self.set_status(self.S_QCRITICAL, msg="return code %s" % self.returncode)

        # Analyze the stderr file for Fortran runtime errors.
        err_msg = None
        if self.stderr_file.exists:
            err_msg = self.stderr_file.read()

        # Start to check ABINIT status if the output file has been created.
        if self.output_file.exists:
            try:
                self.report = self.get_event_report()
            except Exception as exc:
                msg = "%s exception while parsing event_report:\n%s" % (self, exc)
                logger.critical(msg)
                return self.set_status(self.S_ABICRITICAL, msg=msg)

            if self.report.run_completed:

                # Check if the calculation converged.
                not_ok = self.report.filter_types(self.CRITICAL_EVENTS)
                if not_ok:
                    return self.set_status(self.S_UNCONVERGED)
                else:
                    return self.set_status(self.S_OK)

            # Calculation still running or errors?
            if self.report.errors or self.report.bugs:
                # Abinit reported problems
                if self.report.errors:
                    logger.debug('Found errors in report')
                    for error in self.report.errors:
                        logger.debug(str(error))
                        try:
                            self.abi_errors.append(error)
                        except AttributeError:
                            self.abi_errors = [error]

                # The job is unfixable due to ABINIT errors
                logger.debug("%s: Found Errors or Bugs in ABINIT main output!" % self)
                msg = "\n".join(map(repr, self.report.errors + self.report.bugs))
                return self.set_status(self.S_ABICRITICAL, msg=msg)


        # 9) if we still haven't returned there is no indication of any error and the job can only still be running
        # but we should actually never land here, or we have delays in the file system ....
        # print('the job still seems to be running maybe it is hanging without producing output... ')

        # Check time of last modification.
        if self.output_file.exists and \
           (time.time() - self.output_file.get_stat().st_mtime > self.manager.policy.frozen_timeout):
            msg = "Task seems to be frozen, last change more than %s [s] ago" % self.manager.policy.frozen_timeout
            return self.set_status(self.S_ERROR, msg)

        return self.set_status(self.S_RUN)

    # from GsTask
    @property
    def gsr_path(self):
        """Absolute path of the GSR file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._gsr_path
        except AttributeError:
            path = self.outdir.has_abiext("GSR")
            if path: self._gsr_path = path
            return path

    def open_gsr(self):
        """
        Open the GSR file located in the in self.outdir.
        Returns :class:`GsrFile` object, None if file could not be found or file is not readable.
        """
        gsr_path = self.gsr_path
        if not gsr_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a GSR file in %s" % (self, self.outdir))
            return None

        # Open the GSR file.
        from abipy.electrons.gsr import GsrFile
        try:
            return GsrFile(gsr_path)
        except Exception as exc:
            logger.critical("Exception while reading GSR file at %s:\n%s" % (gsr_path, str(exc)))
            return None


class ErrorCode(object):
    UNRECOVERABLE = 'Unrecoverable'
    UNCLASSIFIED = 'Unclassified'
    UNCONVERGED = 'Unconverged'
    INITIALIZATION = 'Initialization'


class AbinitRuntimeError(Exception):
    """
    Exception raised for errors during Abinit calculation
    Initialized with a task, uses it to prepare a suitable error message
    """
    def __init__(self, task):
        self.task = task

    @pmg_serialize
    def to_dict(self):
        report = self.task.report
        d = {}
        d['num_errors'] = report.num_errors
        d['num_warnings'] = report.num_warnings
        if report.num_errors:
            errors = []
            for error in report.errors:
                errors.append({'name': error.name, 'message': error.message})
            d['errors'] = errors
        if report.num_warnings:
            warnings = []
            for warning in report.warnings:
                warnings.append({'name': warning.name, 'message': warning.message})
            d['warnings'] = warnings

        return d