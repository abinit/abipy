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
from pymatgen.io.abinitio.utils import Directory, File
from pymatgen.io.abinitio import events
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.serializers.json_coders import pmg_serialize


logger = logging.getLogger(__name__)


@explicit_serialize
class AbiFireTask(FireTaskBase):

    # List of `AbinitEvent` subclasses that are tested in the check_status method.
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = [
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
        return p.returncode

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
        if self.output_file.exists:
            try:
                self.report = self.get_event_report()
            except Exception as exc:
                msg = "%s exception while parsing event_report:\n%s" % (self, exc)
                logger.critical(msg)
                return {'error': msg}

            if self.report.run_completed:
                # Check if the calculation converged.
                not_ok = self.report.filter_types(self.CRITICAL_EVENTS)
                if not_ok:
                    return 1
                else:
                    return 0

            # Calculation still running or errors?
            if self.report.errors or self.report.bugs:
                return 1


    def run_task(self, fw_spec):
        self.set_workdir(os.path.abspath('.'))
        self.build()
        self.run_abinit(fw_spec)
        retcode = self.task_analysis(fw_spec)
        if retcode:
            raise AbinitRuntimeError(self)
        return FWAction(stored_data=dict(**self.report.as_dict()))


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