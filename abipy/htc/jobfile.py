from __future__ import print_function, division #, unicode_literals

from os import makedirs
from os.path import basename, dirname, join, abspath, exists, realpath
from .utils import listify

__all__ = [
    'JobFile', 
    'PBSJobFile', 
    'SGEJobFile', 
    'SlurmJobFile',
    'MoabJobFile',
]

# =========================================================================== #


class JobFile(object):
    """
    The job file is organized as follow::

        #!/bin/csh                           # 1) Shell specification.
                                             #
        #PBS -N jobname                      # 2) Submission commands.
        #PBS -l walltime=48:00:00            #    Depends on the subclass used.
        #PBS -l nodes=4:ppn=12               
                                             
        set MPIRUN="mpiexec"                 # 3) Declarations.
        set EXECUTABLE=/path/to/executable   #    These are also properties
        set INPUT=calculation.in             #    of the jobfile object.
        set LOG=calculation.log              
                                             
        module load intel-compilers          # 4) Modules.
        module load MPI/Intel/mvapich2       
                                             
        cd ${PBS_O_WORKDIR}                  # 5) Lines before execution.
        limit stacksize unlimited            
                                             
        $MPIRUN $EXECUTABLE < $INPUT > $LOG  # 6) Execution line.
                                             
        echo "Job done!"                     # 7) Lines after execution.
        date                                 

    .. attributes:

        shell:
            The shell binary to be used. E.g. '/bin/csh'.
            Default is '/bin/bash'.
        mpirun:
            The mpi runner. E.g. 'mpiexec -npernode 6'.
            Default is none.
        executable:
            The binary to be executed. Default is abinit.
        bindir:
            The directory in which to look for binaries. Default is none.
        input:
            The input file to feed in the executable as the standard input.
            Mandatory.
        log:
            The file into which the standard output is redirected.
            Default is 'log'.
        stderr:
            The file into which the standard error is redirected.
            Default is 'stderr'.
        modules:
            The modules which will be loaded with 'module load'.
            Default is none.
        lines_before:
            Lines before the main execution.
            Default is none.
        lines_after:
            Lines after the main execution.
            Default is none.
        other_lines:
            Other lines your job submission script would like to have.
            Must be preceded by the approbriate tag (#!, #PBS).
            Default is none.
        submission_command:
            The command which should be used to launch the job.
            E.g. 'qsub', 'bqsub', 'sbatch'.
            Default depends on the job type.
    """
    _executable = 'abinit'
    _mpirun = ''
    _modules = list()
    _other_lines = list()
    _lines_before = list()
    _lines_after = list()

    def __init__(self, name='job.sh', **kwargs):

        # Name
        self.name = name
        self.absdir = realpath(dirname(self.name))
        self.absname = realpath(self.name)

        # Shell
        self._shell = '/bin/bash'

        # Execution lines
        self.input = ''
        self.log = 'log'
        self.stderr = 'stderr'
        self.executable = 'abinit'
        self.bindir = ''
        self.mpirun = ''

        # Modules
        self.modules = list()

        # Other lines
        self.other_lines = list()
        self.lines_before = list()
        self.lines_after = list()

        # Command used to submit the job
        self.submission_command = 'qsub'

        # Set attributes
        for (arg, val) in kwargs.items():
            try:
                getattr(self, 'set_' + arg)(val)
            except:
                pass

    def __str__(self):
        lines = []
        def app(line):
            if '__iter__' in dir(line):
                lines.extend(line)
            else:
                lines.append(line)

        # Shell line
        app('#!' + self.shell)
        app('')

        # Submission instructions
        app(self._get_command_lines())

        # Other submission inscrutions
        app(self.other_lines)
        app('')

        # Declarations
        for (key, value) in [('MPIRUN', self.mpirun),
                             ('EXECUTABLE', self.executable),
                             ('INPUT', self.input),
                             ('LOG', self.log),
                             ('STDERR', self.stderr)]:
            app(self._declare(key, value))
        app('')

        # Modules
        for module in self.modules:
            app('module load ' + module)
        app('')

        # Lines before execution
        app(self.lines_before)
        app('')

        # Execution lines
        if 'csh' in self.shell:
            execline = "($MPIRUN $EXECUTABLE < $INPUT > $LOG) >& $STDERR"
        else:
            execline = "$MPIRUN $EXECUTABLE < $INPUT > $LOG 2> $STDERR"
        app(execline)
        app('')

        # Lines after execution
        app(self.lines_after)
        app('')

        return "\n".join(lines)

    def write(self, name=None):
        """Write the file."""
        if name is None:
            name = self.name

        if self.dirname and not exists(self.dirname):
            makedirs(self.dirname)

        with open(name, 'w') as f:
            f.write(str(self))

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        return list()

    def _declare(self, key, val):
        """Return a lines setting a variable."""
        if 'csh' in self.shell:
            declare = 'set '
        else:
            declare = ''
        return declare + key + '=' + val

    def _set_property(self, name, *args, **kwargs):
        """Set a property through the corresponding set_ function."""
        return getattr(self, 'set_' + key)(*args, **kwargs)

    @property
    def dirname(self):
        """The directory containing the file."""
        return dirname(self.name)

    @property
    def path(self):
        """The path of the file."""
        return abspath(self.name)

    @property
    def basename(self):
        """The base name of the file."""
        return basename(self.name)

    @classmethod
    def properties(cls):
        """Return the list of properties with a set function."""
        funcs = filter(lambda s: s.startswith('set_'), dir(cls))
        return [ f.split('set_', 1)[-1] for f in funcs ]

    @property
    def executable(self):
        return join(self.bindir, self._executable)

    @executable.setter
    def executable(self, executable):
        self._executable = basename(executable)
        if basename(executable) != executable:
            self.set_bindir(dirname(executable))

    @property
    def mpirun(self):
        return '"' + self._mpirun.strip('"').strip("'") + '"'

    @mpirun.setter
    def mpirun(self, mpirun):
        self._mpirun = str(mpirun)

    @property
    def shell(self):
        return self._shell

    @shell.setter
    def shell(self, shell):
        if shell == basename(shell):
            self._shell = join('/bin', shell)
        else:
            self._shell = abspath(shell)

    @property
    def modules(self):
        return self._modules

    @modules.setter
    def modules(self, modules):
        self._modules = listify(modules)

    @property
    def other_lines(self):
        return self._other_lines

    @other_lines.setter
    def other_lines(self, lines):
        self._other_lines = listify(lines)

    @property
    def lines_before(self):
        return self._lines_before

    @lines_before.setter
    def lines_before(self, lines):
        self._lines_before = listify(lines)

    @property
    def lines_after(self):
        return self._lines_after

    @lines_after.setter
    def lines_after(self, lines):
        self._lines_after = listify(lines)

    def set_shell(self, shell):
        """
        Sets the shell type.  The argument can either be an absolute path,
        or just the shell type e.g. bash, csh, tcsh, in which case
        the executable is assumed to be located in /bin/.
        The shell also determine how a variable is declared.
        """
        self.shell = shell

    def set_mpirun(self, mpirun):
        """
        Set the mpi runner to execute the program.
        E.g. 'mpiexec -npernode 6', 'mpirun -np 12', ''.
        """
        self.mpirun = mpirun

    def set_bindir(self, bindir):
        """Set the directory for binaries (abinit, mrgscr...)."""
        self.bindir = realpath(bindir)

    def set_executable(self, executable):
        """Set the executable to use."""
        self.executable = executable

    def set_input(self, input):
        """Set the input file for the main executable."""
        self.input = input

    def set_log(self, log):
        """Set the log file to collect standard output of the executable."""
        self.log = log

    def set_stderr(self, stderr):
        """Set the log file to collect standard output of the executable."""
        self.stderr = stderr

    def set_modules(self, *modules):
        """Set one or many modules to be loaded."""
        self.modules = modules

    def set_other_lines(self, *lines):
        """Set other command lines for the batch submission system."""
        self.other_lines = lines

    def set_lines_before(self, *lines):
        """Set one or many lines to be executed before the main execution."""
        self.lines_before = lines

    def set_lines_after(self, *lines):
        """Set one or many lines to be executed after the main execution."""
        self.lines_after = lines

    def set_submission_command(self, command):
        """
        Sets the command used for job submission,
        e.g. qsub, bqsub, sbatch, ...
        """
        self.submission_command = command

# =========================================================================== #


class PBSJobFile(JobFile):
    """
    Portable Batch System.

    .. attributes:

        jobname:
            Name of the job.
        runtime:
            Maximum time for the job.
        nodes:
            Number of nodes on which to run the job.
        ppn:
            Number of processors per node.
        memory:
            Memory per node. E.g. '48G'.
        queue:
            The queue to which the job is submitted.
        mail:
            The mail to which a notification will be sent.
        mail_options:
            The conditions under which a mail will be sent.
            E.G. 'abe'.
        submission_command:
            default is 'qsub'.

    See man qsub for more info.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#PBS "

    def __init__(self, **kwargs):

        kwargs.setdefault('submission_command', 'qsub')
        JobFile.__init__(self, **kwargs)

    nodes = None
    def set_nodes(self, val):
        self.nodes = val

    ppn = None
    def set_ppn(self, val):
        self.ppn = val

    memory = None
    def set_memory(self, val):
        self.memory = val

    runtime = None
    def set_runtime(self, val):
        """Either set the numer of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.runtime = val

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    queue = None
    def set_queue(self, val):
        self.queue = val

    mail = None
    def set_mail(self, val):
        self.mail = val

    mail_options = None
    def set_mail_options(self, val):
        self.mail_options = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line) # + '\n')

        if self.jobname:
            add('-N ' + str(self.jobname))

        if self.runtime:
            add('-l walltime={0}:{1}:{2}'.format(*self.runtime))

        if self.nodes and self.ppn:
            add('-l nodes=' + str(self.nodes) + ':ppn=' + str(self.ppn))

        if self.memory:
            add('-l mem=' + str(self.memory))

        if self.queue:
            add('-q ' + self.queue)

        if self.mail:
            add('-M ' + self.mail)

        if self.mail_options:
            add('-m ' + self.mail_options)

        return lines

# =========================================================================== #


class SGEJobFile(JobFile):
    """
    Sun Grid Engine.

    .. attributes:

        jobname:
            Name of the job.
        runtime:
            Maximum time for the job.
        nproc:
            Number of processors.
        queue:
            The queue to which the job is submitted.
        environment:
            The parallel environment under which the job is ran.
        memory:
            The requested memory, in M.
        mail:
            The mail to which a notification will be sent.
        mail_options:
            The conditions under which a mail will be sent.
            E.G. 'abe'.
        submission_command:
            default is 'qsub'.

    See man qsub for more info.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#$ "

    def __init__(self, **kwargs):

        kwargs.setdefault('submission_command', 'qsub')
        JobFile.__init__(self, **kwargs)

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    runtime = None
    def set_runtime(self, val):
        """Either set the numer of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.runtime = val

    nproc = None
    def set_nproc(self, val):
        self.nproc = val

    queue = None
    def set_queue(self, val):
        self.queue = val

    environment = None
    def set_environment(self, val):
        self.environment = val

    memory = None
    def set_memory(self, val):
        self.memory = val

    mail = None
    def set_mail(self, val):
        self.mail = val

    mail_options = None
    def set_mail_options(self, val):
        self.mail_options = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line) # + '\n')

        if self.jobname:
            add('-N ' + str(self.jobname))

        if self.runtime:
            add('-l h_rt={0}:{1}:{2}'.format(*self.runtime))

        if self.environment and self.nproc:
            line = '-pe ' + self.environment + ' ' + str(self.nproc)
            if self.memory:
                line += ' -l mem=' + str(self.memory)
            add(line)

        if self.queue:
            add('-q ' + self.queue)

        if self.mail:
            add('-M ' + self.mail)

        if self.mail_options:
            add('-m ' + self.mail_options)

        return lines

# =========================================================================== #


class SlurmJobFile(JobFile):
    """
    Simple Linux Utility for Resource Management.

    .. Attributes:

        jobname:
            Name of the job.
        time:
            Maximum time for the job.
        ntasks:
            The number of processes.
        cpus_per_task:
            The number of cpus per process.
        mem_per_cpu:
            The memory per cpu.
        partition:
            The partition...
        mail_user:
            The mail to which a notification is sent.
        mail_type:
            The conditions unde which to send a mail.
        submission_command:
            default is 'sbatch'.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#SBATCH "

    def __init__(self, **kwargs):

        kwargs.setdefault('submission_command', 'sbatch')
        JobFile.__init__(self, **kwargs)

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    time = None
    def set_time(self, val):
        """Either set the number of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.time = val

    def set_runtime(self, val): 
        self.set_time(val)

    ntasks = None
    def set_ntasks(self, val):
        self.ntasks = val

    ntasks_per_node = None
    def set_ntasks_per_node(self, val):
        self.ntasks_per_node = val

    cpus_per_task = None
    def set_cpus_per_task(self, val):
        self.cpus_per_task = val

    mem_per_cpu = None
    def set_mem_per_cpu(self, val):
        self.mem_per_cpu = val

    partition = None
    def set_partition(self, val):
        self.partition = val

    mail_user = None
    def set_mail_user(self, val):
        self.mail_user = val

    mail_type = None
    def set_mail_type(self, val):
        self.mail_type = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line) # + '\n')

        if self.jobname:
            add('--job-name=' + str(self.jobname))

        if self.time:
            add('--time={0}:{1}:{2}\n'.format(*self.time))

        if self.ntasks:
            add('--ntasks=' + str(self.ntasks))

        if self.partition:
            add('--partition=' + self.partition)

        if self.ntasks_per_node:
            add('--ntasks-per-node=' + str(self.ntasks_per_node))

        if self.cpus_per_task:
            add('--cpus-per-task=' + str(self.cpus_per_task))

        if self.mem_per_cpu:
            add('--mem-per-cpu=' + str(self.mem_per_cpu))

        if self.mail_user:
            add('--mail-user=' + self.mail_user)

        if self.mail_type:
            add('--mail-type=' + self.mail_type)

        return lines

# =========================================================================== #


class MoabJobFile(JobFile):
    """
    Moab Workload Manager

    .. Attributes:
        start_after:
          Declares the time after which the job is eligible for execution.
          Syntax: (brackets delimit optional items with the default being       
          current date/time): [CC][YY][MM][DD]hhmm[.SS]
        account:
            Defines the account associated with the job.
        hold:
            Put a user hold on the job at submission time.
        combine:
            Combine stdout and stderr into the same output file.
        resources:
            Defines the resources that are required by the job.
        mail:
            Defines the set of conditions (a=abort,b=begin,e=end) when the
            server will send a mail message about the job to the user.
        jobname:
            Gives a user specified name to the job.
        priority:
            Assigns a user priority value to a job.
        queue:
            Run the job in the specified queue (pdebug, pbatch, etc.). A host
            may also be specified if it is not the local host.
        rerun:
            Automatically rerun the job is there is a system failure.
        env:
            Specifically adds a list of environment variables that are exported
            to the job.
        allenv:
            Declares that all environment variables in the msub environment are
            exported to the batch job.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#MSUB "

    def __init__(self, **kwargs):

        kwargs.setdefault('submission_command', 'srun')
        JobFile.__init__(self, **kwargs)

    start_after = None
    def set_start_after(self, val):
        self.start_after = val

    account = None
    def set_account(self, val):
        self.account = val

    hold = None
    def set_hold(self, val):
        self.hold = val

    combine = None
    def set_combine(self, val):
        self.combine = val

    resources = dict()
    def set_resources(self, val):
        self.resources = val

    mail = None
    def set_mail(self, val):
        self.mail = val

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    priority = None
    def set_priority(self, val):
        self.priority = val

    queue = None
    def set_queue(self, val):
        self.queue = val

    rerun = None
    def set_rerun(self, val):
        self.rerun = val

    env = None
    def set_env(self, val):
        self.env = val

    allenv = None
    def set_allenv(self, val):
        self.allenv = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line) # + '\n')

        if self.start_after:
            add('-a ' + self.start_after)

        if self.account:
            add('-A ' + self.account)

        if self.hold is True:
            add('-h ')

        if self.combine is True:
            add('-j oe')

        if self.resources:
            for (arg, val) in self.resources.items():
                add('-l ' + arg + '=' + val)

        if self.mail:
            add('-m ' + self.mail)

        if self.jobname:
            add('-N ' + self.jobname)

        if self.priority:
            add('-p ' + self.priority)

        if self.queue:
            add('-q ' + self.queue)

        if self.rerun is True:
            add('-r y')

        if self.env:
            add('-v ' + ','.join(self.env))

        if self.allenv is True:
            add('-V')

        return lines
