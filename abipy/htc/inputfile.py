from __future__ import print_function, division #, unicode_literals

import string
import os.path
import warnings
import numpy as np

from os import makedirs
from os.path import dirname, abspath, exists
from collections import OrderedDict
from copy import deepcopy

from .utils import flatten, listify, is_number, is_iter
from .variable import InputVariable, SpecialInputVariable, _UNITS

__all__ = [
    'InputFile', 
    'VariableBlock',
]


_input_variable_blocks = OrderedDict((
('Datasets', '''
    ndtset jdtset udtset
    '''),
('Basis set', '''
    ecut ecutsm
    '''),
('Bands', '''
    nband nbdbuf
    '''),
('k-point grid', '''
    kptopt nkpt kpt ngkpt kptrlatt
    nshiftk shiftk kptbounds kptns
    '''),
('Models', '''
    ixc ppmodel ppmfreq usepawu upawu jpawu
    '''),
('PAW options', '''
    bxctmindg  dmatpawu   dmatpuopt   dmatudiag  iboxcut
    jpawu  lpawu   lexexch  mqgriddg  ngfftdg
    pawcpxocc   pawcross   pawecutdg   pawfatbnd
    pawlcutd   pawlmix   pawmixdg   pawnhatxc   pawnphi
    pawntheta   pawnzlm   pawoptmix   pawovlp
    pawprtden   pawprtdos   pawprtvol   pawprtwf
    pawspnorb   pawstgylm   pawsushat   pawusecp
    pawxcdev   prtcs   prtefg   prtfc   prtnabla
    ptcharge  quadmom  spnorbscl  usedmatpu   upawu
    useexexch   usepawu   usexcnhat
    '''),
('SCF procedure', '''
    iscf nstep nline tolvrs tolwfr
    toldfe toldff tolimg tolmxf tolrff
    '''),
('KSS generation', '''
    kssform nbandkss
    '''),
('GW procedure', '''
    optdriver gwcalctyp spmeth nkptgw kptgw
    bdgw nqptdm qptdm
    '''),
('GW param', '''
    ecuteps ecutsigx ecutwfn nomegasf
    nfreqim nfreqre freqremax npweps rhoqpmix
    '''),
('GW options', '''
    userre awtr symchi gwpara symsigma gwmem fftgw
    '''),
('Structural optimization', '''
    amu  bmass  delayperm   diismemory   dilatmx   dtion   dynimage
    ecutsm  friction   fxcartfactor  getcell   getxcart   getxred
    goprecon   goprecprm  iatcon   iatfix   iatfixx   iatfixy   iatfixz
    imgmov   ionmov   istatimg  mdtemp   mdwall  natfix   natfixx
    natfixy   natfixz   natcon   nconeq   nimage   nnos   noseinert
    ntime   ntimimage  optcell  pimass   pitransform   prtatlist  qmass
    random_atpos   restartxf  signperm   strfact   strprecon   strtarget
    tolimg   tolmxf  vel   vis  wtatcon
    '''),
('Response function', '''
    bdeigrf elph2_imagden esmear frzfermi
    ieig2rf mkqmem mk1mem prepanl prepgkk
    prtbbb rfasr rfatpol rfddk rfdir rfelfd
    rfmeth rfphon rfstrs rfuser rf1atpol rf1dir
    rf1elfd rf1phon rf2atpol rf2dir rf2elfd
    rf2phon rf3atpol rf3dir rf3elfd rf3phon
    sciss smdelta td_maxene td_mexcit
    '''),
('Wannier 90', '''
    w90iniprj w90prtunk
    '''),
('Parallelisation', '''
    gwpara localrdwf ngroup_rf npband npfft
    npimage   npkpt   npspinor paral_kgb
    paral_rf use_gpu_cuda
    '''),
('Unit cell', '''
    acell angdeg rprim ntypat znucl natom typat xred xcart
    '''),
('Printing', '''
    prtvol enunit
    '''),
('Files', '''
    irdddk irdden ird1den irdqps irdkss irdscr
    irdsuscep irdwfk irdwfq ird1wf getcell
    getddk getden getgam_eig2nkq getkss getocc
    getqps getscr getsuscep getvel getwfk
    getwfq getxcart getxred get1den get1wf
    '''),
))


# =========================================================================== #


class VariableBlock(list):
    """A block of abinit variables."""

    def __init__(self, title, register=''):

        # The block title
        self.title = title

        # A register of all possible input variable.
        if isinstance(register, str):
            self.register = register.split()
        else:
            self.register = list(register)

    def clear(self):
        del self[:]

    def __str__(self):
        lines = ['#== {0} ==#'.format(self.title)]
        for variable in sorted(self):
            svar = str(variable)
            if svar:
                lines.append(svar)
        return '\n'.join(lines)


# =========================================================================== #


class InputFile(object):
    """
    Abinit input file.

    .. code-block:: python

        >> f = InputFile('myfile.in')
        >> f.read('otherfile.in')
        >> 
        >> f.ndtset = 4          # Variables are set with integers,
        >> f.jdtset = [1,2,3,4]  # or lists.
        >> 
        >> f.ecut = 25.          # Here is a floats.
        >> f.ecut = '25.'        # But using strings is always possible.
        >> 
        >> f.tolwfr = 1e-20      # Scientific notation
        >> f.tolwfr = '1d-20'    # is translated like this.
        >> 
        >> f.rprim =  [[0.0,0.5,0.5], [0.5,0.0,0.5], [0.5,0.5,0.0]]  # These three lines
        >> f.rprim =  [ 0.0,0.5,0.5 ,  0.5,0.0,0.5 ,  0.5,0.5,0.0]   # produce exactly
        >> f.rprim ='\\n 0.0 0.5 0.5 \\n 0.5 0.0 0.5 \\n 0.5 0.5 0.0'   # the same result.
        >> 
        >> f.nband4 = 300        # Dataset-specific variable.
        >> 
        >> f.udtset = '2 3'      # We will have to remember that:
        >> f.nband1__s = 100      #   ':' <==> '__s'  (start)
        >> f.nband1__i = 50       #   '+' <==> '__i'  (increment)
        >> f.ecut__a2 = 20.0      #   '?' <==> '__a'  (any)
        >> 
        >> f.istwfk = '*1'       # In some cases, string is the only way!
        >> 
        >> f.tolvrs = None       # Unset a variable. It won't appear in the file.
        >> 
        >> f.fuzzy = 10          # But non-existent variables are written anyway!
        >> 
        >> f.ecut = '100 eV'     # This is OK but not recommended since variables
        >> f.ecut = [100, 'eV']  # are converted to default units when read from file.
        >> 
        >> f.set_comment('''This is a comment.
        ..                  It will be printed at the top of the file.''')
        >> 
        >> f.write()

    See also the function :func:`~abipy.htc.InputFile.set_variables`.
        return self.variables
    """

    _blocks = _input_variable_blocks

    def __init__(self, name='abinit.in'):

        self.__dict__['name'] = str(name)
        self.__dict__['variables'] = dict()
        self.__dict__['_comment'] = str()
        self.__dict__['variables_blocks'] = list()

        for (name, register) in _input_variable_blocks.items():
            self.variables_blocks.append(VariableBlock(name, register))
        self.variables_blocks.append(VariableBlock('Other'))

    def __str__(self):
        lines = list()
        
        # Comments
        if self.comment:
            lines.append(self.comment)
            lines.append('')

        # Clear blocks
        for block in self.variables_blocks:
            block.clear()

        # Sort variables in blocks
        for name, value in self.variables.items():
            variable = SpecialInputVariable(name, value)
            placed = False
            for block in self.variables_blocks:
                if variable.basename in block.register:
                    block.append(variable)
                    placed = True
                    break
            if not placed:
                self.variables_blocks[-1].append(variable)

        # Make the string
        for block in self.variables_blocks:
            if block:
                lines.append(str(block))
                lines.append('')
            block.clear()

        return '\n'.join(lines)

    def __setattr__(self, name, value):
        """
        F.__setattr__('name', value) <==> F.name = value

        Declare a variable in the internal dictionary.
        """
        self.set_variable(name, value)

    def set_name(self, name):
        """Set the name of the file."""
        self.__dict__['name'] = name

    @property
    def path(self):
        """The absolute path."""
        return abspath(self.name)

    @property
    def exists(self):
        """True if self.path exists."""
        return os.path.exists(self.path)

    @property
    def comment(self):
        return self._comment

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        lines = [ '# ' + l.lstrip('#').strip() for l in comment.splitlines() ]
        self.__dict__['_comment'] = '\n'.join(lines)

    def write(self, name=None):
        """Write the inputs to the file."""
        if name is None:
            name = self.name

        if not exists(dirname(self.name)):
            makedirs(dirname(self.name))

        with open(name, 'w') as f:
            f.write(str(self))

    def clear(self):
        """Clear variables."""
        self.variables.clear()
        for block in self.variables_blocks:
            block.clear()

    def read_string(self, bigstring):
        """Initialize all variables from a string."""

        # Split the big string into parts
        parts = list()
        for line in bigstring.splitlines():
            line = line.replace('=', ' ').split('#')[0].strip()
            parts.extend(line.split())

        # Make a list of variable string declaration
        var_list, var_string = list(), ''
        for part in parts:
            if not part:
                continue
            if part[0].isalpha() and part not in _UNITS:
                if var_string:
                    var_list.append(var_string)
                    var_string = ''
            var_string += ' ' + part
        if var_string:
            var_list.append(var_string)

        # Initialize all variables.
        for var_string in var_list:
            variable = SpecialInputVariable.from_str(var_string)
            self.variables[variable.name] = variable.get_value()
                
    @classmethod
    def from_str(cls, bigstring):
        """Initialize from a string."""
        inputfile = cls()
        inputfile.read_string(bigstring)

    def read(self, file):
        """
        Reads the content of an input file and store the variables in the
        internal dictionary with the proper type. Comments are thrown away.
        """
        self.clear()
        with open(file, 'r') as f:
            self.read_string(f.read())

    def set_variable(self, name, value):
        """Set a single variable."""
        self.variables[name] = value

    def set_variables(self, variables=None, dataset=0, **kwargs):
        """
        Sets variables by providing a dictionary, or expanding a dictionary,
        and possibly append them by a dataset index.

        .. code-block:: python

            >> kpoint_grid_shifted = {
            >>     'kptopt' : 1,
            >>     'ngkpt' : 3*[4],
            >>     'nshiftk' : 4,
            >>     'shiftk' : [[0.5,0.5,0.5],
            >>                 [0.5,0.0,0.0],
            >>                 [0.0,0.5,0.0],
            >>                 [0.0,0.0,0.5]],}
            >> 
            >> kpoint_grid_unshifted = {
            >>     'kptopt' : 1,
            >>     'ngkpt' : 3*[4],
            >>     'nshiftk' : 1,
            >>     'shiftk' : [0,0,0],}
            >> 
            >> cell = {
            >>     'ntypat' : 1
            >>     'znucl'  : 6.0
            >>     'natom'  : 2
            >>     'typat'  : [1, 1]
            >>     'xred'   : [[0,0,0],[0.25,0.25,0.25]]
            >>     'acell'  : 3*[6.9]
            >>     'rprim'  : [[0.0,0.5,0.5],
            >>                 [0.5,0.0,0.5],
            >>                 [0.5,0.5,0.0]]}
            >> 
            >> f = InputFile()
            >> f.set_variables(ndtset=3, ecut=4.0, ecutsm=0.5)
            >> 
            >> f.set_variables(cell)    # These two lines
            >> f.set_variables(**cell)  # are equivalent.
            >> 
            >> # Here we append a dataset index at the end of all variables.
            >> f.set_variables(kpoint_grid_shifted, dataset=1)
            >> f.set_variables(kpoint_grid_unshifted, dataset=[2, 3])
            >> 
            >> f.write('myfile.in')  # The name was not set at initialization.
        """
        if variables is None:
            variables = dict()
        variables.update(kwargs)

        if not dataset:
            dataset = ['']

        for ds in listify(dataset):
            for (key, val) in variables.items():
                newkey = key + str(ds)
                self.set_variable(newkey, val)

    def get_variables(self):
        """Return a dictionary of the variables."""
        return deepcopy(self.variables)

    def get_variable(self, variable):
        """Return the value of a variable, or None if it is not set."""
        return self.variables.get(variable)
