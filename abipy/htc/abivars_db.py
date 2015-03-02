"""Database with the names of the input variables used in Abinit and in other main programs."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import json
import yaml
import html2text

from six.moves import cPickle as pickle
from monty.string import is_string, list_strings
from collections import OrderedDict
from monty.functools import lazy_property


with open(os.path.join(os.path.dirname(__file__), "abinit_vars.json")) as fh:
    ABI_VARNAMES = json.load(fh)


# Unit names.
ABI_UNITS = [
    'au',
    'Angstr',
    'Angstrom',
    'Angstroms',
    'Bohr',
    'Bohrs',
    'eV',
    'Ha',
    'Hartree',
    'Hartrees',
    'K',
    'Ry',
    'Rydberg',
    'Rydbergs',
    'T',
    'Tesla',
]

# Operators.
ABI_OPS = ['sqrt', 'end', '*', '/']


# NEW VERSION BASED ON Yannick's database.

list_specials = [
    ('AUTO_FROM_PSP','Means that the value is read from the PSP file'),
    ('CUDA','True if CUDA is enabled (compilation)'),
    ('ETSF_IO','True if ETSF_IO is enabled (compilation)'),
    ('FFTW3','True if FFTW3 is enabled (compilation)'),
    ('MPI_IO','True if MPI_IO is enabled (compilation)'),
    ('NPROC','Number of processors used for Abinit'),
    ('PARALLEL','True if the code is compiled with MPI'),
    ('SEQUENTIAL','True if the code is compiled without MPI'),
]


class literal(str): 
    pass

def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_unicode_representer)


class Variable(yaml.YAMLObject):

    vartype = ''          # String containing the type
    characteristic = None # String containing the characteristics
    definition = None     # String containing the mnemonics
    dimensions = None     # Array containing either int, formula or another variable
    defaultval = None     # Either constant number, formula or another variable
    text = None           # Description (str)
    varname = None        # Name of the variable (str)
    commentdefault = None
    commentdims = None
    section = None
    range = None
    requires = None       # str
    excludes = None       # str

    yaml_tag = u'!variable'
    
    def attrs(self):
        return ['vartype','characteristic','definition','dimensions','defaultval','text',
                'varname','section']

    def __init__(self, vartype=None, characteristic=None,
                definition=None,dimensions=None,default=None,
                text=None, varname=None,section=None,range=None,
                commentdefault=None,commentdims=None):

        self.vartype = vartype
        self.characteristic = characteristic
        self.definition = definition
        self.dimensions = dimensions
        self.defaultval = default
        self.text = literal(text)
        self.varname = varname
        self.section = section
        self.commentdefault = commentdefault
        self.commentdims = commentdims
        self.range = range

        # Fix mnemonics with newlines!
        if self.definition is not None:
            self.definition = self.definition.replace("\n", "")

    @classmethod
    def from_array(cls, array):
        return Variable(vartype=array["vartype"],characteristic=array["characteristic"],
                        definition=array["definition"],dimensions=array["dimensions"],
                        default=array["default"],text=array["text"],varname=array["varname"],
                        section=array["section"],range=array["range"],commentdefault=array["commentdefault"],
                        commentdims=array["commentdims"])

    def __str__(self):
        #return "Variable " + str(self.varname) + " (default = " + str(self.defaultval) + ")"
        return html2text.html2text("<h2>Default value:</h2>" + str(self.defaultval) + "<br/><h2>Description</h2>" +str(self.text))

    def __repr__(self):
        """variable name + mnemonics"""
        return self.varname + "  <" + str(self.definition) + ">"

    def _repr_html_(self):
        """For Ipython notebook"""
        html = "<h2>Default value:</h2>" + str(self.defaultval) + "<br/><h2>Description</h2>" + str(self.text)
        return html.replace("[[", "<b>").replace("]]", "</b>")

    @property
    def info(self):
        """String with Extra info on the variable."""
        attrs = [
            "vartype", "characteristic",  "definition", "dimensions", "defaultval", #"text", 
            "varname", "commentdefault", "commentdims", "section", "range",
            "requires", "excludes"
            ]

        def astr(obj):
            return str(obj).replace("[[", "").replace("]]", "")

        d =  {k: astr(getattr(self, k)) for k in attrs} #if getattr(self, k) is not None}
        from six.moves import StringIO
        stream = StringIO()
        json.dump(d, stream, indent=4, sort_keys=True)
        return stream.getvalue()


class ValueWithUnit(yaml.YAMLObject):
    
    value = None
    units = None
    yaml_tag = u'!valuewithunit'
    
    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units
        
    def __str__(self):
        return str(self.value)+" "+str(self.units)  
    
    def __repr__(self):
        return str(self)
    
def valuewithunit_representer(dumper, data):
    
    return dumper.represent_mapping('!valuewithunit',data.__dict__)
    
class Range(yaml.YAMLObject):
    
    start = None
    stop = None
    
    yaml_tag = u'!range'
    
    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop
        
    def isin(self, value):
        isin = True
        if (start is not None):
            isin = isin and (start <= value)
        if (stop is not None):
            isin = isin and (stop > value)
        return str(self)
    
    def __repr__(self):
        if (self.start is not None and self.stop is not None):
            return "["+str(self.start)+" .. "+str(self.stop)+"]"
        if (self.start is not None):
            return "["+str(self.start)+"; ->"
        if (self.stop is not None):
            return "<-;"+str(self.stop)+"]"
        else:
            return None

class ValueWithConditions(yaml.YAMLObject):
    
    yaml_tag = u'!valuewithconditions'
    
    def __repr__(self):
        s = ''
        for key in self.__dict__.keys():
           if key != 'defaultval':
             s += str(self.__dict__[key])+' if '+str(key)+',\n'
        s += str(self.defaultval) + ' otherwise.\n'
        return s
        
    def __str__(self):
        return self.__repr__()

class MultipleValue(yaml.YAMLObject):

    number = None
    value = None
    
    yaml_tag = u'!multiplevalue'

    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number == None:
            return "*"+str(self.value)
        else:
	        return str(self.number)+"*"+str(self.value)

# Public API
__VARS_DATABASE = None


def get_abinit_variables():
    """Returns the database with the description of the ABINIT variables."""
    global __VARS_DATABASE

    if __VARS_DATABASE is None: 
        pickle_file = os.path.join(os.path.dirname(__file__), "abinit_vars.pickle")
        
        if os.path.exists(pickle_file): 
            #print("Reading from pickle")
            with open(pickle_file, "rb") as fh:
                __VARS_DATABASE = pickle.load(fh)

        else:
            #print("Reading database from YAML file and generating pickle version. It may take a while...")
            from abipy import data as abidata
            yaml_file = abidata.var_file('abinit_vars.yml')
            with open(yaml_file, "r") as fh:
                var_list = yaml.load(fh)

            # Build ordered dict with variables in alphabetical order.
            var_list = sorted(var_list, key=lambda v: v.varname)
            __VARS_DATABASE = VariableDatabase([(v.varname, v) for v in var_list])

            # Save object to pickle file so that can we can reload it from pickle instead of yaml (slower)
            with open(pickle_file, "wb") as fh:
                pickle.dump(__VARS_DATABASE, fh)
                os.chmod(pickle_file, 0o444)

    return __VARS_DATABASE
        

class VariableDatabase(OrderedDict):
    """Stores the mapping varname --> variable object."""

    @lazy_property
    def characteristics(self):
        from abipy import data as abidata
        with open(abidata.var_file('characteristics.yml'),'r') as f:
            return yaml.load(f)
    
    @lazy_property
    def sections(self):
        from abipy import data as abidata
        with open(abidata.var_file('sections.yml'),'r') as f:
            return yaml.load(f)

    def apropos(self, varname):
        """Return the list of variables that are related` to the given varname"""
        vars = []
        for v in self.values():
            if (varname in v.text or 
               (v.dimensions is not None and varname in str(v.dimensions)) or
               (v.requires is not None and varname in v.requires) or
               (v.excludes is not None and varname in v.excludes)):
                vars.append(v)

        return vars

    def vars_with_section(self, sections):
        sections = set(list_strings(sections))
        vars = []
        for v in self.values():
            if v.section in sections:
                vars.append(v)
        return vars

    def vars_with_char(self, chars):
        chars = set(list_strings(chars))
        vars = []
        for v in self.values():
            if v.characteristic: print(v.characteristic)
            if v.characteristic in chars:
                vars.append(v)
        return vars

    def json_dumps_varnames(self):
        return json.dumps(list(self.keys()))


def docvar(varname):
    """Return the `Variable` object associated to this name."""
    return get_abinit_variables()[varname]


def abinit_help(varname, info=True, stream=sys.stdout):
    """Print the abinit documentation on the ABINIT input variable `varname`"""
    database = get_abinit_variables()
    if isinstance(varname, Variable): varname = varname.varname
    try:
        var = database[varname]
    except KeyError:
        stream.write("Variable %s not in the database\n" % varname)
        return
        
    text = html2text.html2text("<h2>Default value:</h2>" + str(var.defaultval) + "<br/><h2>Description</h2>" + str(var.text))
    if info: text += var.info
    stream.write(text.replace("[[", "\033[1m").replace("]]", "\033[0m"))
