from __future__ import print_function, division #, unicode_literals

import subprocess

from os import makedirs, readlink, symlink
from os.path import basename, dirname, exists, join, realpath

from abipy.profile import abipy_env
from .filesfile import FilesFile
from .abinitfiles import AbinitFiles
from .inputfile import InputFile 
from .jobfile import JobFile, PBSJobFile, SGEJobFile, SlurmJobFile

__all__ = ['AbinitInput']

# =========================================================================== #

class AbinitInput(AbinitFiles):
    """
    To Create an Abinit calculation.
    The AbinitInput contains three main internal objects:
    an :class:`~abipy.htc.InputFile`, a :class:`~abipy.htc.FilesFile`,
    and a :class:`~abipy.htc.JobFile`.

    **Arguments**:
        name:
            A string of the form 'directory/rootname' or simply 'directory'.
        jobtype:
            The type of job, e.g. 'PBS', 'SGE', 'Slurm' or None.

    **Keyword arguments for the inputfile**:
        variables:
            A dictionary of input variables for the calculation.

    **Keyword arguments for the filesfile**:
        pseudodir:
            The directory in which to look for the pseudopotential files.
        pseudos:
            A list of pseudopotential files.

    **Keyword arguments for the jobfile**:
        bindir:
            The directory in which to look for binaries.
        executable:
            The binary to be executed.  Default is 'abinit'.
        mpirun:
            The mpi runner.  E.g. 'mpiexec -npernode 6'.
        modules:
            List of modules which will be loaded with 'module load'.
        lines_before:
            List of lines to be executed before the main execution.
        lines_after:
            List of lines to be executed after the main execution.
        other_lines:
            List of command lines for the job submission system.

        Any property of the JobFile or the selected subclass
        (PBSJobFile, SGEJobFile, SlurmJobFile, ...)

    .. code-block:: python

        >> calc = AbinitInput('Mycalc/rootname', jobtype='PBS',
        ..                    pseudodir='Data',
        ..                    bindir='/path/to/binaries/',
        ..                    lines_before=['cd ${PBS_O_WORKDIR}'])
        >>
        >> # Input file properties
        >> calc.ecut = 10.
        >> calc.tolwfr = 1e-8
        >>
        >> calc.kptopt = 1
        >> calc.ngkpt = [2, 2, 2]
        >> calc.nshiftk = 1
        >> calc.shiftk = [0.0, 0.0, 0.0]
        >>
        >> unit_cell = {
        ..     'acell' : 3*[8.6277],
        ..     'rprim' : [[.0, .5, .5],
        ..                [.5, .0, .5],
        ..                [.5, .5, .0]],
        ..     'ntypat' : 2,
        ..     'znucl' : [30, 8],
        ..     'natom' : 2,
        ..     'typat' : [1, 2],
        ..     'xred' : [[.0, .0, .0],
        ..               [.25,.25,.25]]}
        >>
        >> calc.set_variables(unit_cell)
        >>
        >> # Files file properties
        >> calc.set_pseudos('Zn.psp', 'O.psp')
        >>
        >> # Job file properties
        >> calc.set_jobname('MyTest')
        >> calc.set_nodes(2)
        >> calc.set_ppn(12)
        >> calc.set_memory('24gb')
        >> calc.set_runtime(48)
        >> calc.set_mpirun('mpiexec -npernode 6')
        >>
        >> calc.write()

    As you can see, setting an attribute of the AbinitInput object
    will result in setting that input variable in the input file.

    Note also that *all* of the :class:`~abipy.htc.JobFile` functions are available
    though the AbinitInput.
    """

    def __init__(self, name='Calculation/abinit', jobtype='PBS', **kwargs):

        AbinitFiles.__init__(self, name)

        # Initialize the inputfile
        self._setattr(inputfile = InputFile(self.input_name))

        # Initialize the filesfile
        self._setattr(filesfile = FilesFile(self.files_name,
                                            input=self.input_name,
                                            output=self.output_name,
                                            idat_root=self.idat_root,
                                            odat_root=self.odat_root,
                                            tmp_root=self.tmp_root))

        # Initialize the jobfile
        jobargs = dict(name=self.job_name, executable='abinit',
                       input=self.files_name,
                       log=self.log_name) 

        #jobtype = abipy_env.get_uservar("jobtype", kwargs)
        if jobtype is None: jobtype = ''

        if jobtype.lower() == 'pbs':
            self._setattr(jobfile = PBSJobFile(**jobargs))

        elif jobtype.lower() == 'sge':
            self._setattr(jobfile = SGEJobFile(**jobargs))

        elif jobtype.lower() == 'slurm':
            self._setattr(jobfile = SlurmJobFile(**jobargs))

        else:
            self._setattr(jobfile = JobFile(**jobargs))

        # Create setter function for jobfile attributes.
        for prop in self.jobfile.properties():
            function = 'set_' + prop
            self.__dict__[function] = getattr(self.jobfile, function)

        # Pairs of (target, pointer) to be linked.
        self._setattr(_to_link = list())

        # Other arguments
        for key in self.properties():
            val = abipy_env.get_default(key, kwargs)
            if val is not None:
                getattr(self, 'set_' + key)(val)

    def _setattr(self, **kwargs):
        self.__dict__.update(kwargs)

    def __setattr__(self, name, val):
        setattr(self.inputfile, name, val)

    def write(self, *args, **kwargs):
        """
        Write all the files for the calculation.
        """

        force = kwargs.get("force", False)
        if not force and self.inputfile.exists:
            print("Cannot overwrite: %s\tUse -f to force file creation." % self.name)
            return

        # Create directories
        for file in (self.input_name, self.output_name, self.job_name,
                     self.files_name, self.log_name, self.idat_root,
                     self.odat_root, self.tmp_root):
            directory = dirname(file)
            if not exists(directory):
                makedirs(directory)

        # Link files
        for pair in self._to_link:
            self._link(*pair)

        # Write files
        self.inputfile.write()
        self.filesfile.write()
        self.jobfile.write()

    def link_idat(self, file, dtset='auto', datatype='auto'):
        """
        Make a symbolic link for input data.

        Args:
            file:
                The name of the file to be linked.
            dtset:
                The index of the dataset which should read the file.
                By default, it is read from the file name.
                Set to zero for no dataset index.
            datatype:
                The type of datafile, e.g. 'DEN' or 'WFK'.
                By default, it is read from the file name.
        """
        if datatype == 'auto':
            datatype = file.split('_')[-1]

        if dtset == 'auto':
            dtset = file.split('_DS', 1)[-1].split('_')[0]
            try:
                dtset = int(dtset)
            except:
                dtset = 0

        self._to_link.append((file, self.get_idat(datatype, dtset)))

    def link_odat(self, file, dtset='auto', datatype='auto'):
        """
        Make a symbolic link for output data.

        Args:
            file:
                The name of the file to be linked.
            dtset:
                The index of the dataset from which the file belongs.
                By default, it is read from the file name.
                Set to zero for no dataset index.
            datatype:
                The type of datafile, e.g. 'DEN' or 'WFK'.
                By default, it is read from the file name.
        """
        if datatype == 'auto':
            datatype = file.split('_')[-1]

        if dtset == 'auto':
            dtset = file.split('_DS', 1)[-1].split('_')[0]
            try:
                dtset = int(dtset)
            except:
                dtset = 0

        self._to_link.append((file, self.get_odat(datatype, dtset)))

    def link_io(self, idtset, odtset, datatype):
        """
        Make a symbolic link from an output data to an input data.

        Args:
            idtset:
                The dataset index of the input file.
            odtset:
                The dataset index of the output file.
            datatype:
                The type of datafile, e.g. 'DEN' or 'WFK'.
        """
        idat = self.get_idat(datatype, idtset)
        odat = self.get_odat(datatype, odtset)
        self._to_link.append((odat, idat))

    def set_pseudodir(self, pseudodir):
        """Set the directory for the pseudopotentials. Linked to FilesFile."""
        return self.filesfile.set_pseudodir(pseudodir)

    def set_pseudos(self, *pseudopotentials):
        """Set the pseudopotential files. Linked to FilesFile."""
        return self.filesfile.set_pseudos(*pseudopotentials)

    def read(self, *args, **kwargs):
        """Read an input file. Linked to InputFile."""
        return self.inputfile.read(*args, **kwargs)

    def set_variables(self, *args, **kwargs):
        """Set input variables. Linked to InputFile."""
        return self.inputfile.set_variables(*args, **kwargs)

    def set_comment(self, *args, **kwargs):
        """Set a comment in the input file. Linked to InputFile."""
        return self.inputfile.set_comment(*args, **kwargs)

    def properties(self):
        """Return the list of properties with a 'set_' function."""
        funcs = filter(lambda s: s.startswith('set_'), dir(self))
        return [ f.split('set_', 1)[-1] for f in funcs ]

    @staticmethod
    def _link(target, pointer):
        """
        Create the symbolic link  target <-- pointer.
        Overwrite an existing link, but don't overwrite a file.
        """
    
        atarget = realpath(target)
        pointerdir = realpath(dirname(pointer))
        apointer = join(pointerdir, basename(pointer))
    
        if not exists(pointerdir):
            subprocess.call(('mkdir', '-p', pointerdir))
    
        try:
            symlink(atarget, apointer)
        except OSError:
            try:
                oldtarget = realpath(readlink(apointer))
                if oldtarget == atarget:
                    return
                else:
                    subprocess.call(('ln', '-fs', atarget, apointer))
            except OSError:
                raise OSError("Unable to link the files. " +
                              "Maybe {0} exists and is not a link.".format(apointer))

