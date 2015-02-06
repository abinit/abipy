"""Database with the names of the input variables used in Abinit and in other main programs."""
from __future__ import print_function, division, unicode_literals

import os
import json

import cPickle as pickle
#import pickle

import yaml
import html2text
import abipy.data as abidata


with open(os.path.join(os.path.dirname(__file__), "abinit_vars.json")) as fh:
    ABI_VARNAMES = json.load(fh)

# Add internal variables i.e. those variables that are not declared in dtset%
ABI_VARNAMES += [
    "acell",
    "xred",
    "rprim",
    "kptbounds",
    "ndivsm",
    "qpt",
]

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

with open(abidata.var_file('characteristics.yml'),'r') as f:
    list_chars = yaml.load(f);

with open(abidata.var_file('sections.yml'),'r') as f:
    list_sections = yaml.load(f)

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

class literal(str): pass

def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_unicode_representer)


class Variable(yaml.YAMLObject):

    vartype = '' # String containing the type
    characteristic = None # String containing the characteristics
    definition = None # String containing the mnemonics
    dimensions = None # Array containing either int, formula or another variable
    defaultval = None # Either constant number, formula or another variable
    text = None # Description (str)
    varname = None # Name of the variable (str)
    commentdefault = None
    commentdims = None
    section = None
    range = None
    requires = None
    excludes = None

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

    @classmethod
    def from_array(cls,array):
        return Variable(vartype=array["vartype"],characteristic=array["characteristic"],
                        definition=array["definition"],dimensions=array["dimensions"],
                        default=array["default"],text=array["text"],varname=array["varname"],
                        section=array["section"],range=array["range"],commentdefault=array["commentdefault"],
                        commentdims=array["commentdims"])

    def __str__(self):
        return "Variable "+str(self.varname)+" (default = "+str(self.defaultval)+")"


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
        if(start is not None):
            isin = isin and (start <= value)
        if(stop is not None):
            isin = isin and (stop > value)
        return str(self)
    
    def __repr__(self):
        if(self.start is not None and self.stop is not None):
            return "["+str(self.start)+" .. "+str(self.stop)+"]"
        if(self.start is not None):
            return "["+str(self.start)+"; ->"
        if(self.stop is not None):
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
        s+= str(self.defaultval)+' otherwise.\n'
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
    if __VARS_DATABASE is None: __VARS_DATABASE = VariableDatabase()
    return __VARS_DATABASE
        

class VariableDatabase(object):

    def __init__(self):
        pickle_file = os.path.join(os.path.dirname(__file__), "abinit_vars.pickle")
        
        if os.path.exists(pickle_file): 
            print("Reading from pickle")
            with open(pickle_file, "rb") as fh:
                self.all_vars = pickle.load(fh)
        else:
            yaml_file = abidata.var_file('abinit_vars.yml')
            with open(yaml_file, "r") as fh:
                variables = yaml.load(fh)

            self.all_vars = {}
            for var in variables:
                self.all_vars[var.varname] = var

            with open(pickle_file, "wb") as fh:
                pickle.dump(self.all_vars, fh)

    def get_var(self, varname):
        return self.all_vars[varname]


def abinit_help(varname):
    """
    Print the abinit documentation on the ABINIT input variable `varname`
    """
    database = get_abinit_variables()
    try:
        var = database.get_var(varname)
    except KeyError:
        print("Variable %s not in the database" % varname)
        return
        
    text = html2text.html2text("<h2>Default value : </h2>"+str(var.defaultval)+"<br /><h2>Description</h2>"+str(var.text))
    print(text.replace("[[", "\033[1m").replace("]]", "\033[0m"))
