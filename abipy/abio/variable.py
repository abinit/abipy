from __future__ import print_function, division, absolute_import #, unicode_literals

import string
import warnings
import collections
import numpy as np


__all__ = [
    'InputVariable',
]

_SPECIAL_DATASET_INDICES = (':', '+', '?')

_DATASET_INDICES = ''.join(list(string.digits) + list(_SPECIAL_DATASET_INDICES))

_INTERNAL_DATASET_INDICES = ('__s', '__i', '__a')

_SPECIAL_CONVERSION = zip(_INTERNAL_DATASET_INDICES, _SPECIAL_DATASET_INDICES)

_UNITS = {
    'bohr' : 1.0,
    'angstrom' : 1.8897261328856432,
    'hartree' : 1.0,
    'Ha' : 1.0,
    'eV' : 0.03674932539796232,
}


#def convert_number(value):
#    """
#    Converts some object to a float or a string.
#    If the argument is an integer or a float, returns the same object.
#    If the argument is a string, tries to convert to an integer,
#    then to a float.
#    The string '1.0d-03' will be treated the same as '1.0e-03'
#    """
#    if isinstance(value, (float, int)):
#        return value
#
#    elif isinstance(value, str):
#
#        if is_number(value):
#            try:
#                val = int(value)
#            except ValueError:
#                val = float(value)
#            return val
#
#        else:
#            val = value.replace('d', 'e')
#            if is_number(val):
#                val = float(val)
#                return val
#            else:
#                raise ValueError('convert_number failed')
#
#    else:
#        raise ValueError('convert_number failed')


class InputVariable(object):
    """
    An Abinit input variable.
    """
    def __init__(self, name, value, units='', valperline=3):

        self._name = name
        self.value = value
        self._units = units

        # Maximum number of values per line.
        self.valperline = valperline
        if name in ['bdgw']:
            #TODO Shouldn't do that
            self.valperline = 2

        if (is_iter(self.value) and isinstance(self.value[-1], str) and self.value[-1] in _UNITS):
            self.value = list(self.value)
            self._units = self.value.pop(-1)

    def get_value(self):
        """Return the value."""
        if self.units:
            return list(self.value) + [self.units]
        else:
            return self.value

    @property
    def name(self):
        return self._name

    @property
    def basename(self):
        """Return the name trimmed of any dataset index."""
        basename = self.name
        return basename.rstrip(_DATASET_INDICES)

    @property
    def dataset(self):
        """Return the dataset index in string form."""
        return self.name.split(self.basename)[-1]

    @property
    def units(self):
        """Return the units."""
        return self._units

    def __str__(self):
        """Declaration of the variable in the input file."""
        value = self.value
        if value is None or not str(value):
            return ''

        var = self.name
        line = ' ' + var

        # By default, do not impose a number of decimal points
        floatdecimal = 0

        # For some inputs, impose number of decimal points...
        if any(inp in var for inp in ('xred', 'xcart', 'rprim', 'qpt', 'kpt')):
            #TODO Shouldn't do that
            floatdecimal = 16

        # ...but not for those
        if any(inp in var for inp in ('ngkpt', 'kptrlatt', 'ngqpt', 'ng2qpt')):
            #TODO Shouldn't do that
            floatdecimal = 0

        if isinstance(value, np.ndarray):
            n = 1
            for i in np.shape(value):
                n *= i
            value = np.reshape(value, n)
            value = list(value)

        # values in lists
        if isinstance(value, (list, tuple)):

            # Reshape a list of lists into a single list
            if all(isinstance(v, (list, tuple)) for v in value):
                line += self.format_list2d(value, floatdecimal)

            else:
                # Maximum number of values per line.
                #valperline = 3
                #if any(inp in var for inp in ['bdgw']):
                #    #TODO Shouldn't do that
                #    valperline = 2

                line += self.format_list(value, floatdecimal)

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
        if sval.lstrip('-').lstrip('+').isdigit() and floatdecimal == 0:
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
            formatspec = '>{0}d'.format(width)
        elif type_all == str:
            formatspec = '>{0}'.format(width)
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

    def format_list(self, values, floatdecimal=0):
        """
        Format a list of values into a string.
        The result might be spread among several lines.
        """
        line = ''

        # Format the line declaring the value
        for i, val in enumerate(values):
            line += ' ' + self.format_scalar(val, floatdecimal)
            if self.valperline is not None and (i+1) % self.valperline == 0:
                line += '\n'

        # Add a carriage return in case of several lines
        if '\n' in line.rstrip('\n'):
            line = '\n' + line

        return line.rstrip('\n')

    #@staticmethod
    #def string_to_value(sval):
    #    """
    #    Interpret a string variable and attempt to return a value of the
    #    appropriate type.  If all else fails, return the initial string.
    #    """
    #    value = None

    #    try:
    #        for part in sval.split():

    #            if '*' in part:
    #                # cases like istwfk *1
    #                if part[0] == '*':
    #                    value = None
    #                    break

    #                # cases like acell 3*3.52
    #                else:
    #                    n = int(part.split('*')[0])
    #                    f = convert_number(part.split('*')[1])
    #                    if value is None:
    #                        value = []
    #                    value += n * [f]
    #                    continue

    #            # Fractions
    #            if '/' in part:
    #                (num, den) = (float(part.split('/')[i]) for i in range(2))
    #                part = num / den

    #            # Unit
    #            if part in _UNITS.keys():

    #                if value is None:
    #                    warnings.warn("Could not apply the unit token '%s'." % part)
    #                elif isinstance(value, list):
    #                    value.append(part)
    #                else:
    #                    value = [value, part]

    #                # Convert
    #                if False:
    #                    if isinstance(value, list):
    #                        for i in range(len(value)):
    #                            value[i] *= _UNITS[part]
    #                    elif isinstance(value, str):
    #                        value = None
    #                        break
    #                    else:
    #                        value *= _UNITS[part]

    #                continue

    #            # Convert
    #            try:
    #                val = convert_number(part)
    #            except:
    #                val = part

    #            if value is None:
    #                value = val
    #            elif isinstance(value, list):
    #                value.append(val)
    #            else:
    #                value = [value, val]
    #    except:
    #        value = None

    #    if value is None:
    #        value = sval

    #    return value

    #@classmethod
    #def from_str(cls, bigstring):
    #    """Return an instance from a string declaration."""
    #    parts = bigstring.split()

    #    # Perform checks on the string
    #    if len(parts) < 2 or (parts[-1] in _UNITS and len(parts) < 3):
    #        msg = '\n'.join(['Unable to initialize variable from string:',
    #                         bigstring, 'not enough tokens.'])
    #        raise ValueError(msg)
    #    elif not parts[0].isalpha():
    #        msg = '\n'.join(['Unable to initialize variable from string:',
    #                         bigstring, 'no valid variable name found.'])
    #        raise ValueError(msg)

    #    # Make the name
    #    name = parts.pop(0)
    #    #name = cls.declared_to_internal(name)

    #    # Make the units
    #    if parts[-1] in _UNITS:
    #        units = parts.pop(-1)
    #    else:
    #        units = None

    #    value = cls.string_to_value(' '.join(parts))

    #    return cls(name, value, units)


#def is_number(s):
#    """Returns True if the argument can be made a float."""
#    try:
#        float(s)
#        return True
#    except:
#        return False


def is_iter(obj):
    """Return True if the argument is list-like."""
    return hasattr(obj, '__iter__')


def flatten(iterable):
    """Make an iterable flat, i.e. a 1d iterable object."""
    iterator = iter(iterable)
    array, stack = collections.deque(), collections.deque()
    while True:
        try:
            value = next(iterator)
        except StopIteration:
            if not stack:
                return tuple(array)
            iterator = stack.pop()
        else:
            if not isinstance(value, str) \
               and isinstance(value, collections.Iterable):
                stack.append(iterator)
                iterator = iter(value)
            else:
                array.append(value)
