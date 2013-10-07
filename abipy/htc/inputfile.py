
from __future__ import print_function, division

import string
import os.path
import warnings
import numpy as np

from os import makedirs
from os.path import dirname, abspath, exists
from collections import OrderedDict
from copy import deepcopy

from .utils import flatten, listify, is_number, is_iter

__all__ = ['InputFile', 'AbinitVariable', 'VariableBlock']

# =========================================================================== #

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



_units = {
        'bohr' : 1.0,
        'angstrom' : 1.8897261328856432,
        'hartree' : 1.0,
        'Ha' : 1.0,
        'eV' : 0.03674932539796232,
        }

_available_units = [
    'bohr',
    'angstrom',
    'hartree',
    'Ha',
    'eV',
    ]


# =========================================================================== #


def convert_number(value):
    """
    Converts some object to a float or a string.
    If the argument is an integer or a float, returns the same object.
    If the argument is a string, tries to convert to an integer,
    then to a float.
    The string '1.0d-03' will be treated the same as '1.0e-03'
    """

    if isinstance(value, float) or isinstance(value, int):
        return value

    elif isinstance(value, str):

        if is_number(value):
            try:
                val = int(value)
            except ValueError:
                val = float(value)
            return val

        else:
            val = value.replace('d', 'e')
            if is_number(val):
                val = float(val)
                return val
            else:
                raise ValueError('convert_number failed')

    else:
        raise ValueError('convert_number failed')

# =========================================================================== #


class AbinitVariable(object):
    """An Abinit variable."""

    def __init__(self, name, value, units=''):

        # This is the "internal" name, using '_s', '_i' and '_a'.
        self.name = name

        # The value, in any format.
        self.value = value

        # The units
        self.units = units

        if is_iter(self.value) \
        and isinstance(self.value[-1], str) \
        and self.value[-1] in _units:
            self.value = list(self.value)
            self.units = self.value.pop(-1)

    def get_value(self):
        """Return the value."""
        if self.units:
            return list(self.value) + [self.units]
        else:
            return self.value

    @staticmethod
    def internal_to_declared(name):
        """
        Make the conversion
            _s --> :
            _i --> +
            _a --> ?
        """
        for old, new in (('_s', ':'), ('_i', '+'), ('_a', '?')):
            name = name.replace(old, new)
        return name

    @staticmethod
    def declared_to_internal(name):
        """
        Make the conversion
             : --> _s
             + --> _i
             ? --> _a
        """
        for old, new in ((':', '_s'), ('+', '_i'), ('?', '_a')):
            name = name.replace(old, new)
        return name

    @property
    def basename(self):
        """Return the name trimmed of any dataset index."""
        basename = self.name
        for r in ('_s', ':', '_i', '+', '_a', '?'):
            basename = basename[:-4] + basename[-4:].replace(r, '')
        return basename.rstrip(string.digits)

    @property
    def declared_name(self):
        """Return the name for declaration in the input file."""
        return self.internal_to_declared(self.name)

    @property
    def dataset(self):
        """Return the (internal) dataset index in string form."""
        return self.name.split(self.basename)[-1]

    def __str__(self):
        """Declaration of the variable in the input file."""

        value = self.value
        if value is None or not str(value):
            return ''
    
        var = self.declared_name
        line = ' ' + var
    
        # By default, do not impose a number of decimal points
        floatdecimal = 0
    
        # For some inputs, impose number of decimal points...
        if any(inp in var for inp in ('xred','xcart','qpt','kpt')):
            #TODO Shouldn't do that
            floatdecimal = 10
    
        # ...but not for those
        if any(inp in var for inp in ('ngkpt','kptrlatt', 'ngqpt',)):
            #TODO Shouldn't do that
            floatdecimal = 0
    
        if isinstance(value, np.ndarray):
            n = 1
            for i in np.shape(value):
                n *= i
            value = np.reshape(value, n)
            value = list(value)
    
        # values in lists
        if isinstance(value, list):
    
            # Reshape a list of lists into a single list
            if all(isinstance(v, list) for v in value):
                line += self.format_list2d(value, floatdecimal)
    
            else:
                # Maximum number of values per line.
                valperline = 3
                if any(inp in var for inp in ['bdgw']):
                    #TODO Shouldn't do that
                    valperline = 2
    
                line += self.format_list(value, valperline, floatdecimal)
    
        # scalar values
        else:
            line += ' ' + str(value)

        # Add units
        if self.units:
            line += ' ' + self.units
    
        return line

    def format_scalar(self, val, floatdecimal=0):
        """
        Format a single numerical value into a string
        with the appropriate number of decimal.
        """
        sval = str(val)
        if sval.isdigit() and floatdecimal == 0:
            return sval
    
        try:
            fval = float(val)
        except:
            return sval
    
        if fval == 0 or (abs(fval) > 1e-3 and abs(fval) < 1e4):
            form = 'f'
            addlen = 5
        else:
            form = 'e'
            addlen = 8
    
        ndec = max(len(str(fval-int(fval)))-2, floatdecimal)
        ndec = min(ndec, 10)
    
        sval = '{v:>{l}.{p}{f}}'.format(v=fval, l=ndec+addlen, p=ndec, f=form)
    
        sval = sval.replace('e', 'd')
    
        return sval

    def format_list2d(self, values, floatdecimal=0):
        """Format a list of lists."""
    
        lvals = flatten(values)
    
        # Determine the representation
        if all(isinstance(v, int) for v in lvals):
            type_all = int
        else:
            try:
                for v in lvals:
                    float(v)
                type_all = float
            except:
                type_all = str
    
        # Determine the format
        width = max(len(str(s)) for s in lvals)
        if type_all == int:
            formatspec = '>{}d'.format(width)
        elif type_all == str:
            formatspec = '>{}'.format(width)
        else:
    
            # Number of decimal
            maxdec = max(len(str(f-int(f)))-2 for f in lvals)
            ndec = min(max(maxdec, floatdecimal), 10)
    
            if all(f == 0 or (abs(f) > 1e-3 and abs(f) < 1e4) for f in lvals):
                formatspec = '>{w}.{p}f'.format(w=ndec+5, p=ndec)
            else:
                formatspec = '>{w}.{p}e'.format(w=ndec+8, p=ndec)
    
        line = '\n'
        for L in values:
            for val in L:
                line += ' {v:{f}}'.format(v=val, f=formatspec)
            line += '\n'
    
        return line.rstrip('\n')

    def format_list(self, values, valperline=3, floatdecimal=0):
        """
        Format a list of values into a string.
        The result might be spread among several lines.
        """
    
        line = ''
    
        # Format the line declaring the value
        for i, val in enumerate(values):
            line += ' ' + self.format_scalar(val, floatdecimal)
            if (i+1) % valperline == 0:
                line += '\n'
    
        # Add a carriage return in case of several lines
        if '\n' in line.rstrip('\n'):
            line = '\n' + line
    
        return line.rstrip('\n')

    @staticmethod
    def string_to_value(sval):
        """
        Interpret a string variable and attempt to return a value of the
        appropriate type.  If all else fails, return the initial string.
        """

        value = None
    
        try:
            for part in sval.split():
    
                if '*' in part:
    
                    # cases like istwfk *1
                    if part[0] == '*':
                        value = None
                        break
    
                    # cases like acell 3*3.52
                    else:
                        n = int(part.split('*')[0])
                        f = convert_number(part.split('*')[1])
                        if value is None:
                            value = []
                        value += n * [f]
                        continue
    
                # Fractions
                if '/' in part:
                    (num, den) = (float(part.split('/')[i]) for i in range(2))
                    part = num / den
    
                # Unit
                if part in _units.keys():
    
                    if value is None:
                        msg = "Could not apply the unit tolken '%s'." %(part)
                        warnings.warn(msg)
                    elif isinstance(value, list):
                        value.append(part)
                    else:
                        value = [value, part]
    
                    # Convert
                    if False:
                        if isinstance(value, list):
                            for i in range(len(value)):
                                value[i] *= _units[part]
                        elif isinstance(value, str):
                            value = None
                            break
                        else:
                            value *= _units[part]
    
                    continue
    
                # Convert
                try:
                    val = convert_number(part)
                except:
                    val = part
    
                if value is None:
                    value = val
                elif isinstance(value, list):
                    value.append(val)
                else:
                    value = [value, val]
        except:
            value = None
    
        if value is None:
            value = sval
    
        return value

    @classmethod
    def from_str(cls, bigstring):
        """Return an instance from a string declaration."""
        parts = bigstring.split()

        # Perform checks on the string
        if len(parts) < 2 or (parts[-1] in _units and len(parts) < 3):
            msg = '\n'.join(['Unable to initialize variable from string:',
                             bigstring, 'not enough tokens.'])
            raise ValueError(msg)
        elif not parts[0].isalpha():
            msg = '\n'.join(['Unable to initialize variable from string:',
                             bigstring, 'no valid variable name found.'])
            raise ValueError(msg)

        # Make the name
        name = parts.pop(0)
        name = cls.declared_to_internal(name)

        # Make the units
        if parts[-1] in _units:
            units = parts.pop(-1)
        else:
            units = None

        value = cls.string_to_value(' '.join(parts))

        return cls(name, value, units)

    @property
    def sorting_name(self):
        """Name for sorting purposes."""
        conversion = zip(string.digits, string.letters) + (
                     zip(['_s', '_i', '_a'], string.letters[10:]))

        dataset = self.dataset
        for this, that in conversion:
            dataset = dataset.replace(this, that)

        return self.basename + '_' + dataset

    def __gt__(self, other):
        return self.sorting_name > other.sorting_name

    def __lt__(self, other):
        return self.sorting_name < other.sorting_name

    def __eq__(self, other):
        return self.sorting_name == other.sorting_name



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
        lines = ['#== {} ==#'.format(self.title)]
        for variable in sorted(self):
            svar = str(variable)
            if svar:
                lines.append(svar)
        return '\n'.join(lines)


# =========================================================================== #


class InputFile(object):
    """
    Abinit input file.

    Example::

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
        >> f.nband1_s = 100      #   ':' <==> '_s'  (start)
        >> f.nband1_i = 50       #   '+' <==> '_i'  (increment)
        >> f.ecut_a2 = 20.0      #   '?' <==> '_a'  (any)
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
        for variable in self.variables.itervalues():
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
            if part[0].isalpha() and part not in _units:
                if var_string:
                    var_list.append(var_string)
                    var_string = ''
            var_string += ' ' + part
        if var_string:
            var_list.append(var_string)

        # Initialize all variables.
        for var_string in var_list:
            variable = AbinitVariable.from_str(var_string)
            self.variables[variable.name] = variable
                
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
        self.variables[name] = AbinitVariable(name, value)

    def set_variables(self, variables=dict(), dataset=0, **kwargs):
        """
        Sets variables by providing a dictionary, or expanding a dictionary,
        and possibly append them by a dataset index.

        Example::

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
        variables.update(kwargs)

        if not dataset:
            dataset = ['']

        for ds in listify(dataset):
            for (key, val) in variables.items():
                newkey = key + str(ds)
                self.set_variable(newkey, val)

    def get_variables(self):
        """Return a dictionary of the variables."""
        variables = dict()
        for name, var in self.variables:
            variables[name] = var.get_value()
        return variables

    def get_variable(self, variable):
        """Return the value of a variable, or None if it is not set."""
        if variable not in self.variables:
            return None
        return self.variables.get(variable).get_value()

