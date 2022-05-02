# coding: utf-8
"""Pseudopotential Generators."""
from __future__ import annotations

import abc
import os
import tempfile
import collections
import shutil
import time

from typing import Union, Optional # Any,
from monty.os.path import which
from monty.termcolor import cprint
from abipy.flowtk.pseudos import Pseudo
from abipy.ppcodes.oncv_parser import OncvParser

import logging
logger = logging.getLogger(__name__)


# Possible status of the PseudoGenerator.

_STATUS2STR = collections.OrderedDict([
    (1, "Initialized"),    # PseudoGenerator has been initialized
    (2, "Running"),        # PseudoGenerator is running.
    (3, "Done"),           # Calculation done, This does not imply that results are OK
    (4, "Error"),          # PP generator error.
    (5, "Completed"),      # Execution completed successfully.
])


class Status(int):
    """
    An integer representing the status of the 'PseudoGenerator`.
    """
    def __repr__(self) -> str:
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self) -> str:
        """String representation."""
        return _STATUS2STR[self]

    @classmethod
    def as_status(cls, obj: Union[Status, str]) -> Status:
        """Convert obj into Status."""
        if isinstance(obj, cls):
            return obj
        else:
            # Assume string
            return cls.from_string(obj)

    @classmethod
    def from_string(cls, s: str) -> Status:
        """Return an instance from its string representation."""
        for num, text in _STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError(f"Wrong string: `{s}`")


class _PseudoGenerator(metaclass=abc.ABCMeta):
    """
    This object receives a string with the input file and generates a pseudopotential.
    It calls the pp generator in a subprocess to produce the results in a temporary directory.
    It also provides an interface to validate/analyze/plot the results
    produced by the pseudopotential code.

    Concrete classes must:

        1) call super().__init__() in their constructor.
        2) the object should have the input file stored in self.input_str

    Attributes:

        workdir: Working directory (output results are produced in workdir)
        status: Flag defining the status of the ps generator.
        retcode: Return code of the code
        parser: Output parser. None if results are not available because
            the calculations is still running or errors
        pseudo:
            :class:`Pseudo` object. None if not available
    """

    # Possible status
    S_INIT = Status.from_string("Initialized")
    S_RUN = Status.from_string("Running")
    S_DONE = Status.from_string("Done")
    S_ERROR = Status.from_string("Error")
    S_OK = Status.from_string("Completed")

    ALL_STATUS = [
        S_INIT,
        S_RUN,
        S_DONE,
        S_ERROR,
        S_OK,
    ]

    # Basenames for stdin/stdout/stderr.
    stdin_basename = "run.in"

    stdout_basename = "run.out"

    stderr_basename = "run.err"

    def __init__(self, workdir: Optional[str] = None) -> None:
        # Set the initial status.
        self.set_status(self.S_INIT)
        self._parser = None

        if workdir is not None:
            workdir = os.path.abspath(workdir)
            if os.path.exists(workdir):
                raise RuntimeError(f"workdir `{workdir}` already exists")
            self.workdir = workdir

        else:
            # Build a temporary directory
            self.workdir = tempfile.mkdtemp(prefix=self.__class__.__name__)

    def __repr__(self) -> str:
        return "<%s at %s>" % (self.__class__.__name__, self.workdir)

    def __str__(self) -> str:
        #print("Using oncvpsp exec:\n\t", self.executable)
        return "<%s at %s, status=%s>" % (self.__class__.__name__, self.workdir, self.status)

    @property
    def stdin_path(self) -> str:
        """Absolute path of the standard input."""
        return os.path.join(self.workdir, self.stdin_basename)

    @property
    def stdout_path(self) -> str:
        """Absolute path of the standard output."""
        return os.path.join(self.workdir, self.stdout_basename)

    @property
    def stderr_path(self) -> str:
        """Absolute path of the standard error."""
        return os.path.join(self.workdir, self.stderr_basename)

    @property
    def status(self) -> Status:
        """The status of the job."""
        return self._status

    @property
    def retcode(self) -> Union[int, None]:
        """
        Return code of the subprocess. None if not available because e.g. the job has not been started yet.
        """
        try:
            return self._retcode
        except AttributeError:
            return None

    @property
    def parser(self):
        return self._parser

    #@property
    #def pseudo(self) -> Union[Pseudo, None]:
    #    """Pseudo object or None if not available"""
    #    try:
    #        return self._pseudo
    #    except AttributeError:
    #        return None

    @property
    def executable(self) -> str:
        """Name of the executable."""
        return self._executable

    @property
    def input_str(self) -> str:
        """String with the input file."""
        return self._input_str

    def start(self) -> int:
        """"
        Run the calculation in a subprocess (non-blocking interface)
        Return 1 if calculation started, 0 otherwise.
        """
        if self.status >= self.S_RUN:
            return 0

        with open(self.stdin_path, "w") as fh:
            fh.write(self.input_str)

        # Start the calculation in a subprocess and return.
        args = [self.executable, "<", self.stdin_path, ">", self.stdout_path, "2>", self.stderr_path]
        self.cmd_str = " ".join(args)

        from subprocess import Popen, PIPE
        self.process = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.workdir)
        self.set_status(self.S_RUN, info_msg="Start on %s" % time.asctime)

        return 1

    def poll(self) -> int:
        """
        Check if child process has terminated. Set and return returncode attribute.
        """
        self._retcode = self.process.poll()

        if self._retcode is not None:
            self.set_status(self.S_DONE)

        return self._retcode

    def wait(self) -> int:
        """
        Wait for child process to terminate. Set and return returncode attribute.
        """
        self._retcode = self.process.wait()
        self.set_status(self.S_DONE)

        return self._retcode

    def kill(self) -> None:
        """Kill the child."""
        self.process.kill()
        self.set_status(self.S_ERROR)
        self.errors.append("Process has beed killed by host code.")
        self._retcode = self.process.returncode

    def set_status(self, status, info_msg=None):
        """
        Set the status.

        Args:
            status: Status object or string representation of the status
            info_msg: string with human-readable message used in the case of errors (optional)
        """
        assert status in _STATUS2STR
        self._status = status

        if status == self.S_DONE:
            self.check_status()

        return status

    def get_stdin(self) -> str:
        return self.input_str

    def get_stdout(self) -> str:
        """
        Returns a string with the stdout of the calculation.
        """
        with open(self.stdout_path, "rt") as out:
            return out.read()

    def get_stderr(self) -> str:
        """
        Return string with the stderr of the calculation.
        """
        with open(self.stderr_path, "rt") as err:
            return err.read()

    def rmtree(self) -> int:
        """
        Remove the temporary directory. Return exit status
        """
        try:
            shutil.rmtree(self.workdir)
            return 0
        except Exception:
            return 1

    ### ABC PROTOCOL ###

    #@abc.abstractproperty
    #def parser(self):
    #    return _

    @abc.abstractmethod
    def check_status(self):
        """
        This function checks the status of the task by inspecting the output and the
        error files produced by the application
        """


class OncvGenerator(_PseudoGenerator):
    """
    This object receives an input file for oncvpsp, a string
    that defines the type of calculation (scalar-relativistic, ...)
    runs the code in a temporary directory and provides methods
    to validate/analyze/plot the final results.

    Attributes:

        retcode: Retcode of oncvpsp
    """

    @classmethod
    def from_file(cls, path: str, calc_type: str, workdir: Optional[str] = None) -> OncvGenerator:
        """
        Build the object from a file containing the input parameters.
        """
        with open(path, "rt") as fh:
            input_str = fh.read()
            return cls(input_str, calc_type, workdir=workdir)

    def __init__(self, input_str: str, calc_type: str, workdir: Optional[str] = None):
        super().__init__(workdir=workdir)

        self._input_str = input_str
        self.calc_type = calc_type

        self._executable = {
            "non-relativistic": which("oncvpspnr.x"),
            "scalar-relativistic": which("oncvpsp.x"),
            "fully-relativistic": which("oncvpspr.x"),
        }[calc_type]

        if self._executable is None:
            msg = "Cannot find oncvpsp executable in $PATH. Use `export PATH=dir_with_oncvps_executable:$PATH`"
            raise RuntimeError(msg)

    def check_status(self):
        """
        Check the status of the run, set and return self.status attribute.
        """
        if self._status == self.S_OK:
            return self._status

        parser = self._parser = OncvParser(self.stdout_path)

        try:
            parser.scan()
        except parser.Error as exc:
            cprint(str(exc), color="red")
            self._status = self.S_ERROR
            return self._status

        logger.info("run_completed:", parser.run_completed)
        if self.status == self.S_DONE and not parser.run_completed:
            logger.info("Run is not completed!")
            self._status = self.S_ERROR

        if parser.run_completed:
            logger.info("setting status to S_OK")
            self._status = self.S_OK

            psp8_filepath = os.path.join(self.workdir, parser.atsym + ".psp8")

            # Write psp8 file.
            psp8_str = parser.get_psp8_str()
            if psp8_str is not None:
                with open(psp8_filepath, "wt") as fh:
                    fh.write(psp8_str)

            # Add UPF string if present.
            upf_str = parser.get_upf_str()
            if upf_str is not None:
                with open(psp8_filepath.replace(".psp8", ".upf"), "wt") as fh:
                    fh.write(upf_str)

            # Initialize self.pseudo from file.
            pseudo = Pseudo.from_file(psp8_filepath)

            # Add md5 checksum to dojo_report
            if pseudo.has_dojo_report:
                pseudo.dojo_report["md5"] = p.compute_md5()
                pseudo.write_dojo_report(report=p.dojo_report)

        if parser.errors:
            logger.warning("setting status to S_ERROR")
            self._status = self.S_ERROR

        return self._status
