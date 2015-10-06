from __future__ import print_function, division #, unicode_literals

import sys
import os
import warnings
import subprocess
import numpy as np

from os.path import basename, dirname, join, relpath, abspath
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import OrderedDict
from copy import deepcopy

from abipy.core import release, Structure, Density
from abipy.profile import abipy_env
from .utils import parse_ewc
from .abinitinput import AbinitInput


__all__ = [
    'Launcher', 
    'MassLauncher',
]

# =========================================================================== #

class LauncherArgParser(ArgumentParser):
    """
    Base parser used in Launcher to parse the arguments provided by the user
    when the function 'execute' is called.  The parser consists of a top level
    parser responsible for parsing global options such as verbosity level,
    version, etc, and sub-parsers for the different sub-commands.

    Every class that inherits from Launcher, MassLauncher should define
    a parser that  inherits from this base class, and use register_subparser
    to extend or customize the subparsers.  The list of commands supported
    by the instance is stored in self.commands so that we know whether
    a particular sub-parser can handle the argument or if it should delegate
    the superclass.
    """
    def __init__(self, *args, **kwargs):

        #formatter_class=RawDescriptionHelpFormatter)
        ArgumentParser.__init__(self, *args, **kwargs)
        
        # Create the top-level parse responsible for parsing global options
        #   such as verbosity level, version ....
        # NOTE: except for verbose and version, any other option for the top
        #       level parser should be stored in a variable named --top-optname
        #       so that we don't pollute the namespace of the subcommands

        self.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
          help='verbose, can be supplied multiple times to increase verbosity')  

        self.add_argument('--version', action='version', version="abipy " + release.version)

        # Create the parser for the sub-commands
        self.subparsers = self.add_subparsers(dest='command', help='sub-command help')

        # Workaround. see http://stackoverflow.com/questions/8757338/sub-classing-the-argparse-argument-parser
        # I know I shouldn't do this, but I don't want to wrap the parser in the object.
        self.subparsers._parser_class = ArgumentParser

        p_make = self.subparsers.add_parser('make', help='Make files and directories')
        p_make.add_argument('-f', '--force', action='store_true', default=False, help='Force file creation')
        self.register_subparser("make", p_make)

        p_submit = self.subparsers.add_parser('submit', help='Submit the calculation to a batch server.')
        #p_submit.add_argument('bar', type=int, help='bar help')
        self.register_subparser("submit", p_submit)

        p_run = self.subparsers.add_parser('run', help='Run the calculation from the shell.')
        #p_run.add_argument('-n', '--py-nthreads', metavar='NUM', type=int, default=1, help='Number of jobs.')
        self.register_subparser("run", p_run)

        p_report = self.subparsers.add_parser('report', help='Tell if the calculation completed or not.')
        #p_report.add_argument('bar', type=int, help='bar help')
        self.register_subparser("report", p_report)

        p_inspect = self.subparsers.add_parser('inspect', help='Inspect files using EDITOR.')
        p_inspect.add_argument('what_inspect', metavar='character(s)', default = "o", 
          help='Files to inspect: i for the input, o for the output, l for log, j for the job file. f for files file\n' +
          'Characters can be concatenated. Use "ol", for example, to inspect both the output and the log file.')
        self.register_subparser("inspect", p_inspect)

        p_clean = self.subparsers.add_parser('clean', help='Remove log, data files and temporary files.')
        #p_clean.add_argument('-f', '--force', action='store_true', help='Force')
        self.register_subparser("clean", p_clean)

        p_destroy = self.subparsers.add_parser('destroy', help='Remove all files, including input and outputs.')
        p_destroy.add_argument('-f', '--force', action='store_true', help='Force')
        self.register_subparser("destroy", p_destroy)

        p_show = self.subparsers.add_parser('show', help='Print the calculation name and return True.')
        #p_show.add_argument('-f', '--force', action='store_true', help='Force')
        self.register_subparser("show", p_show)

        p_visu = self.subparsers.add_parser('visualize', #aliases=['visu'], 
                                            help='Visualize data with an external program e.g. Xcrysden.')
        p_visu.add_argument('what_visualize', metavar='STRING', default = "crystal", help=' Type of data to visualize')

        self.register_subparser("visualize", p_visu)

    def myparse_args(self, args=None, namespace=None):
        """
        Wrap the parse_args method of ArgumentParsers
        :return: options, args, kwargs

        where options is the default output of parse_args and kwargs is a dictionary option_name -> value
        """
        if args is None:
            self.parse_args(args=["--help"], namespace=namespace)

        if '__iter__' not in dir(args): args = [args]

        # Call the "true" parse_args
        options = self.parse_args(args=args, namespace=namespace)

        args = list()
        kwargs = deepcopy(vars(options))

        return options, args, kwargs

    @property
    def commands(self):
        "The commands registered in the parser"
        return self._cmd2subparser.keys()

    def can_handle(self, command):
        "True if the parser can handle command"
        return command in self.commands

    def iter_cmdsubparser(self):
       "Iterate over (command_string, subparser)"
       for tup in self._cmd2subparser.items(): yield tup

    def register_subparser(self, command, subparser, solve_conflict=False):
        """
        Register the subparser associate to a given command. 

        :arg solve_conflict: By defaut it's not possible to override an existent subparser associated 
                             to the same command. Use solve_conflict if subparser should replace the old one.
        """
        if not hasattr(self, "_cmd2subparser"): self._cmd2subparser = OrderedDict()
        if command in self._cmd2subparser and not solve_conflict:
            raise ValueError("Cannot overwrite subparser for command %s. Use solve_conflict=True" % command)
        self._cmd2subparser[command] = subparser

    def unregister_subparser(self, command):
        """Unregister the subparser associated to the given command. Return the subparser removed"""
        return self._cmd2subparser.pop(command)


# =========================================================================== #

class LauncherError(Exception): 
    """base class for the exceptions raised by Launcher."""

class Launcher(AbinitInput):
    """
    A more powerful version of :class:`~abipy.htc.AbinitInput`.
    Allows to run a calculation, either from a script or from command line.

    .. code-block:: python

        >> calc = Launcher('Mycalc/root',
        ..                 jobtype=None,
        ..                 pseudodir='/path/to/pseudos/',
        ..                 pseudos=['14si.pspnc'],
        ..                 bindir='/path/to/binaries/')
        >>
        >> calc.read('myinput.in')
        >>
        >> # Write the files.
        >> calc.make()
        >>
        >> # Run the calculation with abinit.
        >> calc.run()
        >>
        >> # Inquire about the calculation status.
        >> status = calc.report()
        >>
        >> if status == 'Completed':
        ..     # Remove log and data files.
        ..     calc.clean(force=True)

    You can perform all these actions from the command line, using the function 'execute'.
    """
    Error = LauncherError

    # Parser class and instance are stored as class attributes.
    ArgParser =  LauncherArgParser
    argparser = LauncherArgParser()

    def output_files(self):
        """Return all output files produced, in alphabetical order."""
        base = self.output_name
        files = [base] + [ base + a for a in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' ]
        return filter(os.path.exists, files)

    def last_output(self):
        """Return the last output file produced."""
        files = self.output_files()
        if not files:
            return None
        return files[-1]

    def idat_files(self):
        """Return all the input data files."""
        files = list()
        for file in os.listdir(dirname(self.idat_root)):
            if file.startswith(basename(self.idat_root)):
                files.append(join(dirname(self.idat_root), file))
        return files

    def odat_files(self):
        """Return all the output data files produced."""
        files = list()
        for file in os.listdir(dirname(self.odat_root)):
            if file.startswith(basename(self.odat_root)):
                files.append(join(dirname(self.odat_root), file))
        return files

    def tmp_files(self):
        """Return all the input data files produced."""
        files = list()
        for file in os.listdir(dirname(self.tmp_root)):
            if file.startswith(basename(self.tmp_root)):
                files.append(join(dirname(self.tmp_root), file))
        return files

    def junk_files(self):
        """Return all the junk files produced."""
        files = list()
        for file in os.listdir(self.jobfile.dirname):

            if (file.startswith('fort.') or
                file.endswith('.dat') or
                file in ('_GWDIAG',)):

                files.append(join(self.jobfile.dirname, file))

        return files

    def read_mainlog_ewc(self):
       """
       Read errors, warnings and comments from the main output and the log file.
                                                                                        
       :return: Two namedtuple instances: main, log. 
                The lists of strings with the corresponding messages are
                available in main.errors, main.warnings, main.comments, log.errors etc.
       """
       from pymatgen.io.abinit.events import EventsParser
       parser = EventsParser()
       main_events = parser.parse(self.last_output())
       log_events = parser.parse(self.log_name)

       return main_events, log_events

    def make(self, *args, **kwargs):
        """
        Write the files.

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.get('verbose'):
            print('Writing ' + self.name)

        self.write(*args, **kwargs)

    def run(self, *args, **kwargs):
        """
        Run the calculation by executing the job file from the shell.

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.

        .. warning::
             This method must be thread safe since we may want to run several indipendent
             calculations with different python threads. 
        """
        if kwargs.get('verbose'):
            print('Running ' + self.name + '\n')

        subprocess.call((self.jobfile.shell, self.jobfile.absname))

    def submit(self, *args, **kwargs):
        """
        Submit the calculation to a batch server.

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.

        """
        if kwargs.get('verbose'):
            print('Submitting ' + self.name)

        curdir = abspath(os.curdir)
        os.chdir(self.jobfile.absdir)
        subprocess.call((self.jobfile.submission_command, self.jobfile.basename))
        os.chdir(curdir)

    def clean(self, *args, **kwargs):
        """
        Remove log file, data files and tmp files.

        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.get('verbose'):
            print('Cleaning ' + self.name)

        destroy = [self.log_name]
        for files in (self.odat_files(), self.tmp_files(), self.junk_files()):
            destroy.extend(files)

        if destroy:
            self._remove_files(destroy, kwargs.get('force', False))

    def destroy(self, *args, **kwargs):
        """
        Remove all calculation files and directories, if empty.

        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.

        """
        if kwargs.get('verbose'):
            print('Destroying ' + self.name)

        destroy = [self.input_name, self.log_name, self.job_name, self.files_name]

        for files in (self.output_files(), self.idat_files(),
                      self.odat_files(), self.tmp_files(),  self.junk_files()):
            destroy.extend(files)

        if destroy:
            self._remove_files(destroy, kwargs.get('force', False))
        self._remove_directory_tree(self.dirname)
        
    def show(self, form=str, *args, **kwargs):
        """Print the calculation name and return True."""
        print(form(self.name))
        return True

    def inspect(self, *args, **kwargs):
        """
        Inspect the input/(last) output/ log produced by the run.

        :arg what_inspect: 
                   "i" for the input, "o" for the output, "l" for log, "j" for the job file. "f" for the files file
                   characters can be concatenated. what="ol", for example, will inspect both the output and the log file.

        The environment variable EDITOR defines the application to use (default vi).
        """
        from ..tools import Editor
        editor = Editor()

        what_inspect = kwargs.get("what_inspect", "o")

        filenames = []
        if "i" in what_inspect: filenames.append(self.input_name)
        if "o" in what_inspect: filenames.append(self.last_output())
        if "l" in what_inspect: filenames.append(self.log_name)
        if "j" in what_inspect: filenames.append(self.job_name)
        if "f" in what_inspect: filenames.append(self.files_name)
        editor.edit_files(filenames, ask_for_exit=True)

    def visualize(self, *args, **kwargs):
        # TODO Here I have to decide if this method should be defined
        # in a subclass of Launcher e.g GSLauncher or in the base class
        from .utils import find_file
        what_visualize = kwargs.get("what_visualize", "crystal")

        visualizer = abipy_env.get_uservar("visualizer", kwargs)

        # Find the correct output file
        out_files = self.odat_files() 
        gsfname = find_file(out_files, "GSR")
                                                                          
        if gsfname is None:
            raise RuntimeError("Cannot find GSR file among " % out_files)



        if what_visualize == "crystal":

            structure = Structure.from_file(gsfname)
            structure.visualize(visualizer)()

        elif what_visualize == "density":
            raise NotImplementedError("den_fname?")

            density = Density.from_file(den_fname)
            density.visualize(visualizer)()

        elif what_visualize in ["electrons", "fermisurface",]:
            from ..electrons import ElectronBands
            energies = ElectronBands.from_file(gsfname)

            if what_visualize == "electrons": energies.plot()
            if what_visualize == "fermisurface":
                raise RuntimeError("No hanlder found for fermisurface")
                #visu = energies.visualize(self, visualizer, structure)
                #visu()
        else:
            raise RuntimeError("No handler found for %s" % what_visualize)

    # TODO
    #def __str__(self):
    #   string = ""
    #   return string

    def report(self, *args, **kwargs):
        """
        Print information on the calculation status and return a status.

        Keyword arguments:
            verbose: (0)
                0 : do not print anything
                > 0 : print status
                > 1 : print number of errors, warnings and comments

        """
        output = self.last_output()

        from ..tools import StringColorizer
        str_colorizer = StringColorizer(sys.stdout)
                                                                       
        status2txtcolor = {
          "Completed"  : lambda string : str_colorizer(string, "green"),
          "Unfinished" : lambda string : str_colorizer(string, "blue"),
          "Unstarted"  : lambda string : str_colorizer(string, "cyan"),
        }

        def color(status): return status2txtcolor[status](status)

        verbose = kwargs.get('verbose', 0)

        if output and self._iscomplete(output):
            status = 'Completed'
            msg = relpath(output) + ' : ' + color(status)

            if verbose:
    
                # Does not work!
                pass

                ## Read the number of errors, warnings and comments 
                ##for the (last) main output and the log file.
                #main, log =  self.read_mainlog_ewc()

                #main_info = main.tostream(sys.stdout)
                #log_info  = log.tostream(sys.stdout)

                #msg += "\n  " + "\n  ".join([main_info, log_info])

        elif os.path.exists(self.log_name):
            status = 'Unfinished'
            msg = self.name + ' : ' + color(status)

        else:
            status = 'Unstarted'
            msg = self.name + ' : ' + color(status)

        if verbose:
            print(msg)
            if status == 'Completed':
                pass
                # Does not work!
                #for w in main.warnings: print(w)
                #if verbose > 1:
                #    for w in log.warnings: print(w)

        return status

    def execute(self, *args, **kwargs):
        """
        Execute an action from the command line.

            *  make        --  Write the files.
            *  submit      --  Submit the calculation to a batch server.
            *  run         --  Run the calculation from the shell.
            *  report      --  Tell if the calculation completed or not.
            *  inspect     --  Open files in EDITOR
            *  clean       --  Remove log, data files and temporary files.
            *  show        --  Signify that the calculation exists.
            *  destroy     --  Remove all files, including input and outputs.
            *  visualize   --  Visualize data.

        Suppose this is the content of 'myscript.py':

        .. code-block:: python

            > # myscript.py
            calc = Launcher('Mycalc/root', jobtype=None,
                            pseudodir='Data', pseudos=['14si.pspnc'],
                            executable='abinit')
            calc.read('myinput.in')
            calc.execute()

        Then, from the shell, one can perform the following::

        .. code-block:: python

            > # bash
            > python myscript.py make
            Writing Mycalc/root
            >
            > python myscript.py show
            Mycalc/root
            >
            > python myscript.py run
            Running Mycalc/root
            >
            > python myscript.py report
            Mycalc/root.out : Completed
            >
            > python myscript.py clean
            Cleaning Mycalc/root

        Typically, one would use the command 'submit' instead of 'run'.
        """
        if not args:
            args = sys.argv[1:]

        options, args, kwargs = self.argparser.myparse_args(args)

        if self.argparser.can_handle(options.command):
            getattr(self, options.command)(*args, **kwargs)

        else:
            raise RuntimeError("Don't know how to handle command %s. This should not happen!" % options.command)

    @staticmethod
    def _remove_files(files, force=False):
        """Remove a list of file, asking confirmation."""

        files = filter(os.path.exists, files)

        if files and not force:

            print("About to remove the following files:")
            for file in files:
                print(file)

            proceed = raw_input("Do you want to proceed? (y/n) ")
            if not proceed.lower().startswith('y'):
                return

        for file in files:
            try:
                os.remove(file)
            except Exception as exc:
                warnings.warn(str(exc))

    @staticmethod
    def _remove_directory_tree(topdir):
        """Remove a directory hierarchy, if it contains no files."""
        dirs = [topdir]
        for d in dirs:
            for f in os.listdir(d):
                sub = os.path.join(d, f)
                if os.path.isdir(sub):
                    dirs.append(sub)
                else:
                    return

        for d in reversed(dirs):
            try:
                os.rmdir(d)
            except OSError:
                warnings.warn("Directory tree partially removed: " + topdir)

    @staticmethod
    def _iscomplete(output_file):
        "Return True if an abinit output file is complete."
        with open(output_file, 'read') as f:
            lines = f.readlines()
            lines.reverse()

            for i in range(10):
                try:
                    line = lines[i]
                except:
                    return False

                if 'Calculation completed.' in line:
                    return True
        return False

        # Does not work !
        #from pymatgen.io.abinit.utils import abinit_output_iscomplete
        #return abinit_output_iscomplete(output_file)

# =========================================================================== #

class MassLauncherArgParser(LauncherArgParser):
    """
    top level parser and subparsers used by MassLauncher.
    Handle all the options of Launcher and add the option -c to select the calculations.
    """
    def __init__(self, *args, **kwargs):

        LauncherArgParser.__init__(self, *args, **kwargs)

        for (cmd, subparser) in self.iter_cmdsubparser():
            # Add new option 
            subparser.add_argument('-c', '--calc', dest='only', nargs='*', type=str, 
                                   help="Execute command only for the selected calculations.")

            # Add py_nthreads arg to the commands that support threads.
            if cmd in ["run",]:
                subparser.add_argument('-n', '--py_nthreads', nargs='?', type=int, default=1,
                                       help="The number of threads (to run simultaneously).")


            # Register the curried subparser so that MassLauncher will take over in execute.
            self.register_subparser(cmd, subparser, solve_conflict=True) 

class MassLauncher(object):
    """
    To launch several nearly-identical launchers.
    Acts like a list of launcher.

    .. code-block:: python

        >> # Let's create four calculations           # The rootnames are
        >> calcs = MassLauncher(4, 'Mycalcs/calc',    #   Mycalcs/calc1
        >>                      jobtype='PBS',        #   Mycalcs/calc2
        >>                      pseudodir='Data',     #   Mycalcs/calc3
        >>                      executable='abinit')  #   Mycalcs/calc4
        >>
        >> # == Common properties ==
        >> calcs.read('common.in')
        >>
        >> calcs.ecut = 10.
        >> calcs.tolwfr = 1e-8
        >> calcs.nstep = 0
        >> calcs.iscf = 7
        >>
        >> unit_cell = {'ntypat' : 1, 'znucl' : [14], 'natom' : 2, 'typat' : [1, 1],
        >>              'rprim' : [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]],
        >>              'acell' : 3*[10.261], 'xred' : [[.0, .0, .0], [.25,.25,.25]]}
        >> calcs.set_variables(unit_cell)
        >>
        >> calcs.set_pseudos('14si.pspnc')
        >>
        >> calcs.set_jobname('MyTest')
        >> calcs.set_nodes(1)
        >> calcs.set_ppn(12)
        >> calcs.set_memory('1gb')
        >> calcs.set_runtime(48)
        >>
        >> # == Specific properties ==
        >> ecut = 10.
        >> for calc in calcs:
        >>     calc.ecut = ecut
        >>     ecut += 5.
        >>
        >> # Write them all.
        >> calcs.make()
    """

    ArgParser = MassLauncherArgParser
    argparser = MassLauncherArgParser()

    def __init__(self, n=0, name='Calc', *args, **kwargs):

        self._setattr(launchers = list())
        self._setattr(ids = list())

        if n > 0:
            index_format = '0=' + str(int(np.log10(n)) + 1)
            for i in range(1, n+1):
                index = '{i:{f}}'.format(i=i, f=index_format)
                calc_name = name + index
                launcher = Launcher(calc_name, *args, **kwargs)
                self.add_launcher(launcher, index=index)

    def __getitem__(self, i): return self.launchers[i]

    def __setitem__(self, i, calc): self.launchers[i] = calc

    def __delitem__(self, i): del self.launchers[i]

    def __iter__(self): return iter(self.launchers)

    def __len__(self): return len(self.launchers)

    def _setattr(self, **kwargs):
        self.__dict__.update(kwargs)

    #def __getattr__(self, name):
    #    """Return a function that passes the arguments to all launchers."""
    #    def f(*args, **kwargs):
    #        return [ getattr(c, name)(*args, **kwargs) for c in self ]
    #    return f

    def _distributed(self, func_name):
        """Return a function that passes the arguments to all launchers."""
        def f(*args, **kwargs):
            return [ getattr(c, func_name)(*args, **kwargs) for c in self ]
        f.__doc__ =  getattr(self[0], func_name).__doc__
        return f

    def _make_distributed(f):
        """Make a function distributed to all launchers."""
        def g(self, *args, **kwargs):
            return self._distributed(f.__name__)(*args, **kwargs)
        g.__doc__ = getattr(Launcher, f.__name__).__doc__
        return g

    def __setattr__(self, name, value):
        """Set an attribute to all launchers."""
        return [ setattr(calc, name, value) for calc in self ]

    def add_launcher(self, launcher, index=None):
        """Add a Launcher instance (or derivatives) to the list."""
        if index is None:
            index = str(len(self.launchers) + 1)
        # maybe needs an OrderedDict here
        self.launchers.append(launcher)
        self.ids.append(index)

        # Create missing property setter
        if len(self.launchers) == 1:
            for prop in self.launchers[0].properties():
                setter = 'set_' + prop
                if not setter in dir(self):
                    self.__dict__[setter] = self._distributed(setter)

    def properties(self):
        """Return the list of properties with a `set_` function."""
        funcs = filter(lambda s: s.startswith('set_'), dir(self))
        return [ f.split('set_', 1)[-1] for f in funcs ]

    def only(self, only=None):
        """Return a list of indexes and a list of calc which are included in 'only'."""
        if only:
            filtered = list()
            for i, calc in zip(self.ids, self):
                if str(i) in map(str, only):
                    filtered.append((i, calc))
        else:
            filtered= zip(self.ids, self)
        return filtered

    @_make_distributed
    def set_pseudodir(self): return

    @_make_distributed
    def set_pseudos(self): return

    @_make_distributed
    def read(self): return

    @_make_distributed
    def set_variables(self): return

    @_make_distributed
    def set_comment(self): return

    @_make_distributed
    def link_idat(self): return

    @_make_distributed
    def link_odat(self): return

    @_make_distributed
    def link_io(self): return

    def execute(self, *args, **kwargs):
        """
        Execute an action given from the command line.

            *  make        --  Write the files.
            *  submit      --  Submit the calculation to a batch server.
            *  run         --  Run the calculation from the shell.
            *  report      --  Tell if the calculation completed or not.
            *  clean       --  Remove log, data files and temporary files.
            *  show        --  Signify that the calculation exists.
            *  destroy     --  Remove all files, including input and outputs.

        Command line optional arguments:
            -c [id1 [,id2 [id3, ... ]]] :

        Select a subset of calculations.

        With the previous example, one could issue::

            > python myscript.py make
            Writing Mycalc/calc1
            Writing Mycalc/calc2
            Writing Mycalc/calc3
            Writing Mycalc/calc4
            >
            > python myscript.py show -c 1 2
            Mycalc/calc1
            Mycalc/calc2
            >
            > python myscript.py run -c 2 3
            Running Mycalc/calc2
            Running Mycalc/calc3
            >
            > python myscript.py report
            Mycalc/calc1 : Unstarted
            Mycalc/calc2.out : Completed
            Mycalc/calc3.out : Completed
            Mycalc/calc4 : Unstarted
        """
        if not args:
            args = sys.argv[1:]

        options, args, kwargs = self.argparser.myparse_args(args)
        kwargs.update(vars(options))

        if self.argparser.can_handle(options.command):

            try:
                nthreads = options.py_nthreads
            except AttributeError:
                nthreads = 1
                
            print("About to run command", options.command," with nthreads", nthreads)

            if nthreads == 1:

                if options.command in dir(self):
                    getattr(self, options.command)(*args, **kwargs)
                else:
                    for index, calc in self.only(options.only):
                    #for i, calc in zip(self.ids, self):
                    #    if options.only and str(i) not in map(str, options.only):
                    #        continue
                        getattr(calc, options.command)(*args, **kwargs)

            else:
                # Threaded version.
                from threading import Thread
                from Queue import Queue

                def worker():
                    while True:
                        func, args, kwargs = q.get()
                        func(*args, **kwargs)
                        q.task_done()
                                                             
                q = Queue()
                for i in range(nthreads):
                    t = Thread(target=worker)
                    t.setDaemon(True)
                    t.start()
                                                             
                #for (i, calc) in zip(self.ids, self):
                #    if options.only and i not in options.only: continue
                for index, calc in self.only(options.only):
                    func = getattr(calc, options.command)
                    q.put((func, args, kwargs))

                # Block until all tasks are done. 
                q.join()  

        else:
            raise RuntimeError("Don't know how to handle command %s. This should not happen!" % options.command)

    @_make_distributed
    def report(self): return

    @_make_distributed
    def odat_files(self): return

    @_make_distributed
    def last_output(self): return

    #@_make_distributed
    def show(self, form=None, *args, **kwargs):
        only = kwargs.get('only')
        if form is None:
            form = str
        
        #@form
        def tmpform(s):
            return str(index) + ' ' + s

        newform = lambda s: form(tmpform(s))

        for index, calc in self.only(only):
            calc.show(form=newform, *args, **kwargs)

        return

    @_make_distributed
    def make(self): return

    @_make_distributed
    def run(self): return

    @_make_distributed
    def submit(self): return

    @_make_distributed
    def clean(self): return

    @_make_distributed
    def destroy(self): return
