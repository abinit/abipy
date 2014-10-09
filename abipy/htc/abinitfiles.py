"""
Base class for an ABINIT calculation.
"""
from __future__ import print_function, division #, unicode_literals
from os.path import basename, dirname, join, splitext, realpath

__all__ = ["AbinitFiles"]


# =========================================================================== #

class AbinitFiles(object):
    """
    A manager for the abinit file names.

    A calculation is organised as follow::

        CalcDir/
        |-- CalcRoot.in
        |-- CalcRoot.out
        |-- input_data/
        |   |-- idat_CalcRoot_DS2_DEN
        |-- run/
        |   |-- CalcRoot.files
        |   |-- CalcRoot.sh
        |   |-- CalcRoot.log
        |   |-- out_data/
        |   |   |-- odat_CalcRoot_DS1_DEN
        |   |   |-- odat_CalcRoot_DS1_WFK
        |   |-- tmp_data/
        |   |   |-- tmp_CalcRoot_LOG_0001
        |   |   |
    """

    _directory = {
        'inp'   : '',
        'out'   : '',
        'job'   : 'run',
        'files' : 'run',
        'log'   : 'run',
        'idat'  : 'input_data',
        'odat'  : join('run', 'out_data'),
        'tmp'   : join('run', 'tmp_data')}

    _prefix = {
        'inp'   : '',
        'out'   : '',
        'job'   : '',
        'files' : '',
        'log'   : '',
        'idat'  : 'idat_',
        'odat'  : 'odat_',
        'tmp'   : 'tmp_'}

    _suffix = {
        'inp'   : '.in',
        'out'   : '.out',
        'files' : '.files',
        'log'   : '.log',
        'idat'  : '',
        'odat'  : '',
        'tmp'   : '',
        'job'   : '.sh'}

    def _form_name(self, filetype):
        d = self._directory[filetype]
        pr = self._prefix[filetype]
        sf = self._suffix[filetype]
        return join(self.absdir, d, pr + self.basename + sf)

    def __init__(self, rootname='AbinitCalculation/root'):

        # Take the root name
        if '.' in rootname:
            rootname = splitext(rootname)[0]

        # The calculation must be in a directory.
        if rootname == basename(rootname):
            rootname = join(rootname, 'calc')

        self.set_rootname(rootname)

    def set_rootname(self, name):
        """Set the root name."""
        self.__dict__['rootname'] = name
        self.__dict__['absdir'] = realpath(dirname(name))

    @property
    def dirname(self):
        """The top-level directory."""
        return dirname(self.rootname)

    @property
    def basename(self):
        """The basename of the root name."""
        return basename(self.rootname)

    @property
    def name(self):
        """The root name."""
        return self.rootname

    @property
    def files_name(self):
        """The name of the ".files" file given to abinit."""
        return self._form_name('files')

    @property
    def job_name(self):
        """The submission script name."""
        return self._form_name('job')

    @property
    def input_name(self):
        """The input file name."""
        return self._form_name('inp')

    @property
    def log_name(self):
        """The name of the log file."""
        return self._form_name('log')

    @property
    def idat_root(self):
        """The root name for input data files."""
        return self._form_name('idat')

    @property
    def odat_root(self):
        """The root name for output data files."""
        return self._form_name('odat')

    @property
    def tmp_root(self):
        """The root name for temporaty data files."""
        return self._form_name('tmp')

    @property
    def output_name(self):
        """The output file produced (based on the name only)."""
        return self._form_name('out')

    def get_odat(self, datatype, dtset=0):
        """
        Returns an output data file name.

        Args:
            datatype:
                The type of datafile, e.g. 'DEN' or 'WFK'.

            dtset:
                The dataset index from which to take the data file.
                If 0 (the default), no dataset index is used.
        """
        file = self.odat_root

        if int(dtset) > 0: file += '_DS' + str(dtset)

        file += '_' + datatype.upper().lstrip('_')

        return file

    def get_idat(self, datatype, dtset=0):
        """
        Returns an input data file name.

        Args:
            datatype:
                The type of datafile, e.g. 'DEN' or 'WFK'.

            dtset:
                The dataset index from which to take the data file.
                If 0 (the default), no dataset index is used.
        """
        file = self.idat_root

        if int(dtset) > 0: file += '_DS' + str(dtset)

        file += '_' + datatype.upper().lstrip('_')

        return file

    def get_netcdf(self, datatype, dtset=0):
        """
        Returns a netcdf output data file name.

        Args:
            datatype:
                The type of datafile, 'DEN' or 'WFK'.

            dtset:
                The dataset index from which to take the data file.
                If 0 (the default), no dataset index is used.
        """
        return self.get_odat(datatype, dtset=dtset) + '.nc'
