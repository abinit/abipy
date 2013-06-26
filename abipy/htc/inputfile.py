from __future__ import print_function, division

import string
import os.path
import warnings
import numpy as np

from os import makedirs
from os.path import dirname, abspath, exists
from collections import OrderedDict
from copy import deepcopy

__all__ = ['InputFile']

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
    """

    _blocks = _input_variable_blocks

    def __init__(self, name='abinit.in'):

        self.__dict__['name'] = str(name)
        self.__dict__['path'] = abspath(name)
        self.__dict__['variables'] = dict()
        self.__dict__['_comment'] = str()

    def __str__(self):
        return str(self.variables)

    def __setattr__(self, name, val):
        """
        F.__setattr__('name', value) <==> F.name = value

        Stores values inside F.variables dictionary.

        To set some attribute, and not an input variable,
        use F.__dict__['name'] = value
        """
        self.variables[name] = val

    @property
    def exists(self):
        "True if self.path exists."
        return os.path.exists(self.path)

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        comment = str(comment).lstrip('#\n').rstrip('\n')
        comment = '# ' + comment
        comment = comment.replace('\n', '\n# ')
        comment = comment.replace('\n# #', '\n# ')
        comment += '\n'
        self.__dict__['_comment'] = comment

    def write(self, name=None):
        """Write the inputs to the file."""
        if name is None:
            name = self.name

        if not exists(dirname(self.name)):
            makedirs(dirname(self.name))

        variables = deepcopy(self.variables)

        with open(name, 'w') as f:

            # Write comments
            if self._comment:
                f.write(self._comment)

            # Write defined input blocks
            for (block_name, block_var) in self._blocks.items():

                block_dict = dict()
                for var in block_var.split():
                    for var_instance in self._variable_instances(var):
                        val = variables.pop(var_instance, None)
                        if val is not None:
                            block_dict[var_instance] = val

                if block_dict:
                    lines = '#== {} ==\n'.format(block_name)
                    lines += format_variable_dict(block_dict)
                    lines += '\n'
                    f.write(lines)

            # Write other input variables
            if variables:
                lines = '#== {} ==\n'.format('Other')
                lines += format_variable_dict(variables)
                lines += '\n'
                f.write(lines)

    def read(self, file):
        """
        Reads the content of an input file and store the variables in the
        internal dictionary with the proper type. Comments are thrown away.
        """
        self.variables.clear()

        # Split the file into a dictionary of variables
        #   and a string containing their value.
        var = None
        val = ''
        tmp_dict = dict()
        with open(file, 'r') as f:

            for line in f.readlines():

                line = line.replace('=', ' ').split('#')[0].strip()

                for part in line.split():

                    # New input variable
                    if part[0].isalpha() and part not in _units:

                        # Store old variable
                        if var:
                            tmp_dict[var] = val

                        # Reset variable and value
                        var = expand_var(part)
                        val = ''

                    # append to value
                    else:
                        val += ' ' + part

            # Store last variable
            if var and val != '':
                tmp_dict[var] = val

        # Interpret values
        for var, sval in tmp_dict.iteritems():
            self.set_variable(var, string_to_value(sval))

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
            dataset = ''

        for ds in listify(dataset):
            for (key, val) in variables.items():
                newkey = key + str(ds)
                self.set_variable(newkey, val)

    def get_variables(self):
        """Return the variables dict."""
        return self.variables

    def get_variable(self, variable):
        """Return the value of a variable, or None if it is not set."""
        return self.variables.get(variable)

    def set_variable(self, variable, value):
        """Set a single variable."""
        self.variables[variable] = value

    def _variable_instances(self, var):
        """
        Return the list of all instances of a particular input variable defined,
        appended or not by a digit.  Return the keys only, sorted.
        """

        inputlist = list()
        for inp in self.variables.keys():
            if strip_var(inp) == var:
                inputlist.append(inp)

        sorting = list()
        for inp in inputlist:
            digit = inp.split(var)[-1]
            n = 0

            if digit == '':
                n += 9000000

            elif digit.isdigit():
                n += 8000000 + int(digit)

            elif digit.count('_') == 2 and len(digit) == 4:
                n += 7000000

                if digit[:2] == '_s':   n += 10
                elif digit[:2] == '_i': n += 20
                elif digit[:2] == '_a': n += 30

                if digit[2:] == '_s':   n +=  1
                elif digit[2:] == '_i': n +=  2
                elif digit[2:] == '_a': n +=  3

            elif digit.startswith('_') and digit[2:].isdigit():
                n += 6000000

                if digit[:2] == '_s':   n += 100000
                elif digit[:2] == '_i': n += 200000
                elif digit[:2] == '_a': n += 300000

                n +=  int(digit[2:])

            elif digit[-2:].startswith('_') and digit[:-2].isdigit():
                n += 5000000

                if digit[-2:] == '_s':   n += 100000
                elif digit[-2:] == '_i': n += 200000
                elif digit[-2:] == '_a': n += 300000

                n += int(digit[:-2])

            sorting.append((inp,n))

        # Sort the inputs
        sorting = sorted(sorting, cmp=lambda x,y: cmp(x[1],y[1]))
        inputlist = [ imp[0] for imp in sorting ]

        return inputlist


# =========================================================================== #


def strip_var(var):
    """Strip an input of any trailling dataset indication."""
    inp = var.rstrip(string.digits)
    for r in ('_s', ':', '_i', '+', '_a', '?'):
        inp = inp[:-4] + inp[-4:].replace(r, '')
    return inp

def expand_var(inputvar):
    """
    This is a two-way substitution function between
        ':' <==> '_s'
        '+' <==> '_i'
        '?' <==> '_a'
    Also, if the trailing digits are '0' or '00', suppresses it.
    """
    inputvar = str(inputvar)

    alphapart = inputvar.rstrip(string.digits)
    rest = inputvar.split(alphapart)[-1]
    if not rest:
        pass

    elif rest.isdigit() and int(rest) == 0:
        inputvar = alphapart

    elif rest.isdigit():
        pass

    replacelist = (('_s', ':'), ('_i', '+'), ('_a', '?'))
    if any(r[1] in inputvar for r in replacelist):
        for (a, b) in replacelist:
            inputvar = inputvar.replace(b, a)

    elif any(r[0] in inputvar for r in replacelist):
        for (a, b) in replacelist:
            inputvar = inputvar[:-4] + inputvar[-4:].replace(a, b)

    return inputvar

def format_var(var, value):
    """
    Return a string for a single input variable
    to be written in the input file.
    Always ends with a carriage return.
    """

    if value is None or not str(value):
        return ''

    var = expand_var(var)
    line = ' ' + var

    # By default, do not impose a number of decimal points
    floatdecimal = 0

    # For some inputs, impose number of decimal points...
    if any(inp in var for inp in ('xred','xcart','qpt','kpt')):
        floatdecimal = 10

    # ...but not for those
    if any(inp in var for inp in ('ngkpt','kptrlatt', 'ngqpt',)):
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
            line += format_list_of_list(value, floatdecimal)

        else:
            # Maximum number of values per line.
            valperline = 3
            if any(inp in var for inp in ['bdgw']):
                valperline = 2

            line += format_list(value, valperline, floatdecimal)

    # scalar values
    else:
        line += ' ' + str(value) + '\n'

    return line

def format_scalar(val, floatdecimal=0):
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

def format_list_of_list(values, floatdecimal=0):
    """Format a list of lists."""

    lvals = flattens(values)

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

    return line

def format_list(values, valperline=3, floatdecimal=0):
    """
    Format a list of values into a string.  The result might be spread among
    several lines, each line starting with a blank space.
    'valperline' specifies the number of values per line.
    """

    line = ''

    # Format the line declaring the value
    for i, val in enumerate(values):
        line += ' ' + format_scalar(val, floatdecimal)
        if (i+1) % valperline == 0:
            line += '\n'

    # Add a carriage return in case of several lines
    if '\n' in line.rstrip('\n'):
        line = '\n' + line
    # Add a carriage return at the end
    if not line.endswith('\n'):
        line += '\n'

    return line


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

def format_variable_dict(variables):
    """Format a dictionary of input variables into a multi-line string."""
    string = ''
    for (var, val) in variables.iteritems():
        string += format_var(var, val)
    return string


def flattens(lists):
    """Transform a list of lists into a single list."""
    if not isinstance(lists[0], list):
        return lists
    L = list()
    for l in lists:
        L += l
    return L

def listify(obj):
    """Transform any object, iterable or not, to a list."""
    if '__iter__' in dir(obj):
        return list(obj)
    else:
        return [obj]

def is_number(s):
    """Returns True if the argument can be made a float."""
    try:
        float(s)
        return True
    except:
        return False

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

