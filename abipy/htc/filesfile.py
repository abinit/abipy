from __future__ import print_function, division #, unicode_literals
import warnings

from os import makedirs
from os.path import dirname, join, exists, realpath
from copy import deepcopy

from .abinitfiles import AbinitFiles

__all__ = ['FilesFile']

# =========================================================================== #


class FilesFile(object):
    """
    A .files file used as input for abinit.

    **Keyword arguments:**

        input: Input file.
        output: Output file.
        idat_root: Root for input data.
        odat_root: Root for output data.
        tmp_root: Root for temporary files.
        pseudodir: The directory for the pseudopotentials.
        pseudos: List of pseudopotential files.

    .. code-block:: python

        >> filesfile = Filesfile('run/abinit.files',
        ..                       input='abinit.in',
        ..                       output='abinit.out',
        ..                       idat_root='run/data_files/idat_',
        ..                       odat_root='run/data_files/odat_',
        ..                       tmp_root='run/data_files/tmp_',)
        ..                       pseudodir='/path/to/pseudopotential/directory',)
        >> 
        >> filesfile.set_pseudos(['H.psp', 'O.psp'])
        >> filesfile.write()
        >> with open('run/abinit.files', 'r') as f: print f.read()
        ../CalcRoot.in
        ../CalcRoot.out
        data_files/idat_
        data_files/odat_
        tmp_files/tmp_
        /path/to/pseudopotential/directory/H.psp
        /path/to/pseudopotential/directory/O.psp
    """

    def __init__(self, name='calc.files', **kwargs):

        self.name = name

        self.input = 'calc.in'
        self.output = 'calc.out'
        self.idat_root = 'idat_calc'
        self.odat_root = 'odat_calc'
        self.tmp_root = 'tmp_calc'

        self.pseudodir = '.'
        self.pseudos = list()

        for (arg, val) in kwargs.items():
            getattr(self, 'set_' + arg)(val)

    def set_input(self, name):
        """Set the input file."""
        self.input = name

    def set_output(self, name):
        """Set the output file."""
        self.output = name

    def set_idat_root(self, name):
        """Set root for input data."""
        self.idat_root = name

    def set_odat_root(self, name):
        """Set root for output data."""
        self.odat_root = name

    def set_tmp_root(self, name):
        """Set root for temporary files."""
        self.tmp_root = name

    def set_pseudodir(self, directory):
        """Set the directory for the pseudopotentials."""
        self.pseudodir = realpath(directory)

    def set_pseudos(self, *pseudos):
        """
        Sets the pseudopotential files.

        **Arguments:**
            pseudos:
                Pseudopotential files.

        Both of these syntax work::
            >> f.set_pseudos('H.psp', 'O.psp')
            >> f.set_pseudos(['H.psp', 'O.psp'])
        """
        if not pseudos:
            self.pseudos = list()
            return

        elif '__iter__' in dir(pseudos[0]):
            pseudos = tuple(pseudos[0])

        self.pseudos = pseudos

    @property
    def pseudos(self):
        """Pseudopotential files."""
        return [ join(self.pseudodir, pseudo) for pseudo in self._pseudos ]

    @pseudos.setter
    def pseudos(self, vals):
        self._pseudos = deepcopy(vals)

    @property
    def dirname(self):
        """The directory containing the file."""
        return dirname(self.name)

    def check_pseudos(self):
        """Issue a warning for each pseudopotential file not found."""
        for pseudo in self.pseudos:
            if not exists(pseudo):
                warnings.warn('Pseudopotential file not found: ' + pseudo)

    def __str__(self):
        lines = [self.input, self.output, self.idat_root, self.odat_root,
                 self.tmp_root] + self.pseudos
        return '\n'.join(lines) + '\n'

    def write(self):
        """Write the file."""
        self.check_pseudos()

        if self.dirname and not exists(self.dirname):
            makedirs(self.dirname)

        with open(self.name, 'w') as f:
            f.write(str(self))
