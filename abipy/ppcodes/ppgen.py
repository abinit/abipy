# coding: utf-8
"""Pseudopotential Generators."""
from __future__ import annotations

import abc
import os
import tempfile
import collections
import shutil
import time

from typing import Any, Union, Optional
from itertools import product
from monty.os.path import which
from abipy.flowtk.pseudos import Pseudo
from abipy.ppcodes.oncvpsp import OncvOutputParser

import logging
logger = logging.getLogger(__name__)


# Possible status of the PseudoGenerator.
_STATUS2STR = collections.OrderedDict([
    (1, "Initialized"),    # PseudoGenerator has been initialized
    (2, "Running"),        # PseudoGenerator is running.
    (3, "Done"),           # Calculation done, This does not imply that results are ok
    (4, "Error"),          # Generator error.
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
    def from_string(cls, s):
        """Return a :class:`Status` instance from its string representation."""
        for num, text in _STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)


class PseudoGenerator:
    """
    This object receives a string with the input file and generates a pseudopotential.
    It calls the pp generator in a subprocess to produce the results in a temporary directory.
    It also provides an interface to validate/analyze/plot the results
    produced by the pseudopotential code. Concrete classes must:

        1) call super().__init__() in their constructor.
        2) the object should have the input file stored in self.input_str

    Attributes:
        workdir: Working directory (output results are produced in workdir)
        status: Flag defining the status of the ps generator.
        retcode: Return code of the code
        errors: List of strings with errors.
        warnings: List of strings with warnings.
        parser: Output parser. None if results are not available because
            the calculations is still running or errors
        results: Dictionary with the most important results. None if results are not available because
            the calculations is still running or errors
        pseudo:
            :class:`Pseudo` object. None if not available
    """
    __metaclass__ = abc.ABCMeta

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

    stdin_basename = "run.in"
    stdout_basename = "run.out"
    stderr_basename = "run.err"

    def __init__(self, workdir: Optional[str] = None) -> None:
        # Set the initial status.
        self.set_status(self.S_INIT)
        self.errors, self.warnings = [], []

        if workdir is not None:
            workdir = os.path.abspath(workdir)
            if os.path.exists(workdir):
                raise RuntimeError("workdir %s already exists" % workdir)
            self.workdir = workdir
        else:
            # Build a temporary directory
            self.workdir = tempfile.mkdtemp(prefix=self.__class__.__name__)

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
        Return code of the subprocess. None if not available because e.g. the job
        has not been started yet.
        """
        try:
            return self._retcode
        except AttributeError:
            return None

    @property
    def pseudo(self):
        """:class:`Pseudo` object."""
        try:
            return self._pseudo
        except AttributeError:
            return None

    @property
    def executable(self) -> str:
        """Name of the executable."""
        return self._executable

    @property
    def input_str(self) -> str:
        """String with the input file."""
        return self._input_str

    def __repr__(self) -> str:
        return "<%s at %s>" % (self.__class__.__name__, os.path.basename(self.workdir))

    def __str__(self) -> str:
        #return "<%s at %s, status=%s>" % (self.__class__.__name__, os.path.basename(self.workdir), self.status)
        return "<%s at %s, status=%s>" % (self.__class__.__name__, self.workdir, self.status)

    def start(self) -> int:
        """"
        Run the calculation in a sub-process (non-blocking interface)
        Return 1 if calculation started, 0 otherwise.
        """
        print(self.executable)
        if self.status >= self.S_RUN:
            return 0

        logger.info("Running in %s:" % self.workdir)
        with open(self.stdin_path, "w") as fh:
            fh.write(self.input_str)

        # Start the calculation in a subprocess and return.
        args = [self.executable, "<", self.stdin_path, ">", self.stdout_path, "2>", self.stderr_path]
        self.cmd_str = " ".join(args)

        from subprocess import Popen, PIPE
        self.process = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.workdir)
        self.set_status(self.S_RUN, info_msg="Start on %s" % time.asctime)

        return 1

    def poll(self):
        """Check if child process has terminated. Set and return returncode attribute."""
        self._retcode = self.process.poll()

        if self._retcode is not None:
            self.set_status(self.S_DONE)

        return self._retcode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._retcode = self.process.wait()
        self.set_status(self.S_DONE)

        return self._retcode

    def kill(self):
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

        #changed = True
        #if hasattr(self, "_status"):
        #    changed = (status != self._status)

        self._status = status

        if status == self.S_DONE:
            self.check_status()

        #if status == self.S_OK:
        #    self.on_ok()

        return status

    @abc.abstractmethod
    def check_status(self):
        """
        This function checks the status of the task by inspecting the output and the
        error files produced by the application
        """

    def get_stdin(self):
        return self.input_str

    def get_stdout(self) -> str:
        """Returns a string with the stdout of the calculation."""
        if not os.path.exists(self.stdout_path):
            return "Stdout file does not exist"

        with open(self.stdout_path) as out:
            return out.read()

    def get_stderr(self) -> str:
        """Returns a string with the stderr of the calculation."""
        if not os.path.exists(self.stdout_path):
            return "Stderr file does not exist"

        with open(self.stderr_path) as err:
            return err.read()

    def rmtree(self) -> int:
        """Remove the temporary directory. Return exit status"""
        try:
            shutil.rmtree(self.workdir)
            return 0
        except Exception:
            return 1

    #def on_ok(self):
    #    """
    #    Method called when calculation reaches S_OK
    #    Perform operations to finalize the run. Subclasses should provide their own implementation.
    #    """

    @abc.abstractmethod
    def plot_results(self, **kwargs):
        """
        Plot the results with matplotlib.
        """

    def parse_output(self):
        """
        Uses OutputParser to parse the output file.the output file.the output file.the output file.
        """
        parser = self.OutputParser(self.stdout_path)
        try:
            parser.scan()
        except parser.Error:
            time.sleep(1)
            try:
                parser.scan()
            except parser.Error:
                raise

    @property
    def results(self):
        return getattr(self, "_results", None)

    @property
    def plotter(self):
        """
        Return a plotter object that knows how to plot the results
        of the pseudopotential generator
        """
        return getattr(self, "_plotter", None)


class OncvGenerator(PseudoGenerator):
    """
    This object receives an input file for oncvpsp, a string
    that defines the type of calculation (scalar-relativistic, ...)
    runs the code in a temporary directory and provides methods
    to validate/analyze/plot the final results.

    Attributes:

        retcode: Retcode of oncvpsp
    """
    OutputParser = OncvOutputParser

    @classmethod
    def from_file(cls, path: str, calc_type: str, workdir: Optional[str] = None) -> OncvGenerator:
        """
        Build the object from a file containing input parameters.
        """
        with open(path, "rt") as fh:
            input_str = fh.read()
            return cls(input_str, calc_type, workdir=workdir)

    def __init__(self, input_str: str, calc_type: str, workdir: Optional[str] = None):
        super(OncvGenerator, self).__init__(workdir=workdir)
        self._input_str = input_str
        self.calc_type = calc_type

        calctype2exec = {
            "non-relativistic": which("oncvpspnr.x"),
            "scalar-relativistic": which("oncvpsp.x"),
            "fully-relativistic": which("oncvpspr.x")}

        self._executable = calctype2exec[calc_type]
        if self.executable is None:
            msg = "Cannot find executable for oncvpsp is PATH. Use `export PATH=dir_with_executable:$PATH`"
            raise RuntimeError(msg)

        self.format = 'psp8'

    def check_status(self):
        """Check the status of the run, set and return self.status attribute."""
        if self.status == self.S_OK:
            return self._status

        parser = self.OutputParser(self.stdout_path)

        try:
            parser.scan()
        except parser.Error:
            self._status = self.S_ERROR
            return self._status

        logger.info("run_completed:", parser.run_completed)
        if self.status == self.S_DONE and not parser.run_completed:
            logger.info("Run is not completed!")
            self._status = self.S_ERROR

        if parser.run_completed:
            logger.info("setting status to S_OK")
            self._status = self.S_OK
            #########################################
            # Here we initialize results and plotter.
            #########################################
            if parser.warnings:
                self.errors.extend(parser.warnings)

            try:
                self._results = parser.get_results()
            except parser.Error:
                # File may not be completed.
                time.sleep(2)
                try:
                    self._results = parser.get_results()
                except Exception:
                    raise

            self._plotter = parser.make_plotter()

            # Write Abinit pseudopotential.
            filepath = os.path.join(self.workdir, parser.atsym + "." + self.format)
            self.filepath = filepath
            #if os.path.exists(filepath):
            #    raise RuntimeError("File %s already exists" % filepath)

            if self.format == 'psp8':

                # Write psp8 file.
                psp8_str = parser.get_psp8_str()
                if psp8_str is not None:
                    with open(filepath, "wt") as fh:
                        fh.write(psp8_str)

                # Add upf string (if present).
                upf_str = parser.get_upf_str()
                if upf_str is not None:
                    with open(filepath.replace(".psp8", ".upf"), "wt") as fh:
                        fh.write(upf_str)

                # Initialize self.pseudo from file.
                self._pseudo = p = Pseudo.from_file(filepath)

                # Add md5 checksum to dojo_report
                if p.has_dojo_report:
                    p.dojo_report["md5"] = p.compute_md5()
                    p.write_dojo_report(report=p.dojo_report)

        if parser.errors:
            logger.warning("setting status to S_ERROR")
            self._status = self.S_ERROR
            self.errors.extend(parser.errors)
            print(self.errors)

        return self._status

    def plot_results(self, **kwargs):
        """Plot the results with matplotlib."""
        #if not self.status == self.S_OK:
        #    logger.warning("Cannot plot results. ppgen status is %s" % self.status)
        #    return

        # Call the output parser to get the results.
        parser = self.OutputParser(self.stdout_path)
        parser.scan()

        # Build the plotter and plot data according to **kwargs
        plotter = parser.make_plotter()
        plotter.plot_atanlogder_econv()


class OncvMultiGenerator:
    """
    This object receives a template input file and generates multi
    pseudos by changing particular parameters.
    """
    def __init__(self, filepath: str, calc_type: str = "scalar-relativistic") -> None:
        """
        Args:
            filepath: File with the input file
        """
        self.filepath = os.path.abspath(filepath)
        self.calc_type = calc_type

        with open(filepath, "r") as fh:
            self.template_lines = fh.readlines()

    def change_icmod3(self, fcfact_list=(3, 4, 5), rcfact_list=(1.3, 1.35, 1.4, 1.45, 1.5, 1.55)):
        """
        Change the value of fcfact and rcfact in the template. Generate the new pseudos
        and create new directories with the pseudopotentials in the current working directory.

        Return: List of `Pseudo` objects

        Old version with icmod == 1.

        # icmod fcfact
        1 0.085

        New version with icmod == 3.
        # icmod, fcfact (rcfact)
            3    5.0  1.3
        """
        magic = "# icmod fcfact"
        for i, line in enumerate(self.template_lines):
            if line.strip() == magic: break
        else:
            raise ValueError("Cannot find magic line `%s` in template:\n%s" % (magic, "\n".join(self.template_lines)))

        # Extract the parameters from the line.
        pos = i + 1
        line = self.template_lines[pos]

        tokens = line.split()
        icmod = int(tokens[0])

        #if len(tokens) != 3:
        #    raise ValueError("Expecting line with 3 numbers but got:\n%s" % line)
        #icmod, old_fcfact, old_rcfact = int(tokens[0]), float(tokens[1]), float(tokens[2])
        #if icmod != 3:
        #    raise ValueError("Expecting icmod == 3 but got %s" % icmod)

        base_name = os.path.basename(self.filepath).replace(".in", "")
        ppgens = []
        for fcfact, rcfact in product(fcfact_list, rcfact_list):
            new_input = self.template_lines[:]
            new_input[pos] = "%i %s %s\n" % (3, fcfact, rcfact)
            input_str = "".join(new_input)
            #print(input_str)
            ppgen = OncvGenerator(input_str, calc_type=self.calc_type)

            name = base_name + "_fcfact%3.2f_rcfact%3.2f" % (fcfact, rcfact)
            ppgen.name = name
            ppgen.stdin_basename = name + ".in"
            ppgen.stdout_basename = name + ".out"

            # Attach fcfact and rcfact to ppgen
            ppgen.fcfact, ppgen.rcfact = fcfact, rcfact

            if not ppgen.start() == 1:
                raise RuntimeError("ppgen.start() failed!")
            ppgens.append(ppgen)

        for ppgen in ppgens:
            retcode = ppgen.wait()
            ppgen.check_status()

        # Ignore errored calculations.
        ok_ppgens = [gen for gen in ppgens if gen.status == gen.S_OK]
        print("%i/%i generations completed with S_OK" % (len(ok_ppgens), len(ppgens)))

        ok_pseudos = []
        for ppgen in ok_ppgens:
            # Copy files to dest
            pseudo = ppgen.pseudo
            #dest = os.path.basename(self.filepath) + "_fcfact%3.2f_rcfact%3.2f" % (ppgen.fcfact, ppgen.rcfact)
            dest = os.path.split(self.filepath)[0]
            shutil.copy(os.path.join(ppgen.workdir,ppgen.stdin_basename), dest)
            shutil.copy(os.path.join(ppgen.workdir,ppgen.stdout_basename), dest)

            # Reduce the number of ecuts in the DOJO_REPORT
            # Re-parse the output and use devel=True to overwrite initial psp8 file
            psp8_path = os.path.join(dest, ppgen.name + ".psp8")
            out_path = os.path.join(dest, ppgen.name + ".out")

            parser = OncvOutputParser(out_path)
            parser.scan()

            # Rewrite pseudo file in devel mode.
            with open(psp8_path, "w") as fh:
                fh.write(parser.get_pseudo_str(devel=True))

            # Build new pseudo.
            p = Pseudo.from_file(psp8_path)
            ok_pseudos.append(p)

        return ok_pseudos
