"""Database with the names of the input variables used in Abinit and in other main programs."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import json
import yaml
import html2text

from collections import OrderedDict, defaultdict
from six.moves import cPickle as pickle
from monty.string import is_string, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint


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
    ('AUTO_FROM_PSP', 'Means that the value is read from the PSP file'),
    ('CUDA', 'True if CUDA is enabled (compilation)'),
    ('ETSF_IO', 'True if ETSF_IO is enabled (compilation)'),
    ('FFTW3', 'True if FFTW3 is enabled (compilation)'),
    ('MPI_IO', 'True if MPI_IO is enabled (compilation)'),
    ('NPROC', 'Number of processors used for Abinit'),
    ('PARALLEL', 'True if the code is compiled with MPI'),
    ('SEQUENTIAL', 'True if the code is compiled without MPI'),
]


class literal(str):
    pass


def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_unicode_representer)


class Variable(yaml.YAMLObject):
    yaml_tag = u'!variable'

    def __init__(self, vartype=None, characteristics=None, definition=None, dimensions=None,
                default=None, text=None, varname=None, section=None, range=None,
                commentdefault=None, commentdims=None, requires=None, excludes=None):
        """
        Args:
            vartype: String containing the type
            characteristic: String containing the characteristics
            definition: String containing the mnemonics
            dimensions: Array containing either int, formula or another variable
            default: Either constant number, formula or another variable
            text: Description (str)
            varname: Name of the variable (str)
            commentdefault:
            commentdims:
            section:
            range:
            requires:
            excludes:
        """
        self.vartype = vartype
        self.characteristics = characteristics
        self.definition = definition
        self.dimensions = dimensions
        self.defaultval = default
        #self.text = literal(text)
        self.varname = varname
        self.section = section
        self.range = range
        self.commentdefault = commentdefault
        self.commentdims = commentdims

        try:
            self.text = unicode(text)  # py2
        except:
            self.text = str(text, "utf-8")

        # Fix mnemonics with newlines!
        if self.definition is not None:
            self.definition = self.definition.replace("\n", "")

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    @property
    def name(self):
        """An alias for self.varname"""
        return self.varname

    def __repr__(self):
        """variable name + mnemonics"""
        return self.varname + "  <" + str(self.definition) + ">"

    def __str__(self):
        s = html2text.html2text("<h2>Default value:</h2>" + str(self.defaultval) + "<br/><h2>Description</h2>" + str(self.text))
        return str(s.encode("utf-8", errors="ignore"))

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        html = "<h2>Default value:</h2>" + str(self.defaultval) + "<br/><h2>Description</h2>" + str(self.text)
        return html.replace("[[", "<b>").replace("]]", "</b>")

    @property
    def info(self):
        """String with Extra info on the variable."""
        attrs = [
            "vartype", "characteristics",  "definition", "dimensions", "defaultval", #"text",
            "varname", "commentdefault", "commentdims", "section", #"range",
            "requires", "excludes"
            ]

        def astr(obj):
            return str(obj).replace("[[", "").replace("]]", "")

        d = {k: astr(getattr(self, k)) for k in attrs} #if getattr(self, k) is not None}
        return json.dumps(d, indent=4, sort_keys=True)

    @lazy_property
    def isarray(self):
        """True if the variable is an array."""
        if (self.dimensions is None or
            # FIXME: Improve treatment of ValueWithConditions.
            #isinstance(self.dimensions, (Range, ValueWithConditions)) or
            isinstance(self.dimensions, (Range,)) or
            is_string(self.dimensions)
            ): return False
        return True

    def depends_on_dimension(self, dimname):
        """
        True if variable is an array whose shape depends on dimension name `dimname`.

        Args:
            dimname: String of :class:`Variable` object.
        """
        if isinstance(dimname, Variable): dimname = dimname.varname
        if not self.isarray: return False
        # This test is not very robust and can fail.
        if dimname in str(self.dimensions): return True
        return False

    @property
    def url(self):
        """The url associated to the variable."""
        # TODO: root will change once we move to the new website.
        root = "https://www.abinit.org/sites/default/files/last/input_variables/html_automatically_generated/"
        return root + "%s.html#%s" % (self.section, self.varname)

    def html_link(self, label=None):
        """String with the URL of the web page."""
        label = self.varname if label is None else label
        return '<a href="%s" target="_blank">%s</a>' % (self.url, label)

    def browse(self):
        """Open variable documentation in browser."""
        import webbrowser
        return webbrowser.open(self.url)



class ValueWithUnit(yaml.YAMLObject):
    yaml_tag = u'!valuewithunit'

    value = None
    units = None

    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units

    def __str__(self):
        return str(self.value) + " " + str(self.units)

    def __repr__(self):
        return str(self)


def valuewithunit_representer(dumper, data):
    return dumper.represent_mapping('!valuewithunit', data.__dict__)


class Range(yaml.YAMLObject):
    yaml_tag = u'!range'

    start = None
    stop = None

    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop

    def isin(self, value):
        isin = True
        if self.start is not None:
            isin = isin and self.start <= value
        if self.stop is not None:
            isin = isin and self.stop > value
        return str(self)

    def __repr__(self):
        if self.start is not None and self.stop is not None:
            return "[" + str(self.start) + " .. " + str(self.stop) + "]"
        if self.start is not None:
            return "[" + str(self.start) + "; ->"
        if self.stop is not None:
            return "<-;" + str(self.stop) + "]"
        else:
            return None


class ValueWithConditions(yaml.YAMLObject):
    yaml_tag = u'!valuewithconditions'

    def __repr__(self):
        s = ''
        for key in self.__dict__.keys():
            if key != 'defaultval':
                s += str(self.__dict__[key]) + ' if ' + str(key) + ',\n'
        s += str(self.defaultval) + ' otherwise.\n'
        return s

    def __str__(self):
        return self.__repr__()


class MultipleValue(yaml.YAMLObject):
    yaml_tag = u'!multiplevalue'

    number = None
    value = None

    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number is None:
            return "*" + str(self.value)
        else:
            return str(self.number) + "*" + str(self.value)

__VARS_DATABASE = None


##############
# Public API #
##############


def get_abinit_variables():
    """Returns the database with the description of the ABINIT variables."""
    global __VARS_DATABASE

    if __VARS_DATABASE is None:
        pickle_file = os.path.join(os.getenv("HOME"), ".abinit", "abipy", "abinit_vars.pickle")

        if os.path.exists(pickle_file):
            try:
                with open(pickle_file, "rb") as fh:
                    __VARS_DATABASE = pickle.load(fh)
            except Exception as exc:
                cprint("Error while trying to read variables from pickle file: %s" % pickle_file, "red")
                cprint("Please remove the file and rerun.", "red")
                raise exc

        else:
            # Make directory and file if not present.
            if not os.path.exists(os.path.dirname(pickle_file)):
                os.makedirs(os.path.dirname(pickle_file))

            #print("Reading database from YAML file and generating pickle version. It may take a while...")
            from abipy import data as abidata
            yaml_path = abidata.var_file('abinit_vars.yml')
            __VARS_DATABASE = VariableDatabase.from_file(yaml_path)

            # Save object to pickle file so that can we can reload it from pickle instead of yaml (slower)
            with open(pickle_file, "wb") as fh:
                pickle.dump(__VARS_DATABASE, fh)

    return __VARS_DATABASE


class VariableDatabase(OrderedDict):
    """Stores the mapping varname --> variable object."""

    @classmethod
    def from_file(cls, yaml_path):
        with open(yaml_path, "rt") as fh:
            var_list = yaml.load(fh)

        # Build ordered dict with variables in alphabetical order.
        var_list = sorted(var_list, key=lambda v: v.varname)
        return cls([(v.varname, v) for v in var_list])

    @lazy_property
    def characteristics(self):
        """List of characteristics."""
        from abipy import data as abidata
        with open(abidata.var_file('characteristics.yml'), 'rt') as f:
            return yaml.load(f)

    @lazy_property
    def sections(self):
        """List of sections"""
        from abipy import data as abidata
        with open(abidata.var_file('sections.yml'), 'rt') as f:
            return yaml.load(f)

    @lazy_property
    def name2section(self):
        """
        Dictionary mapping the name of the variable to the section.
        """
        d = {}
        for name, var in self.items():
            d[name] = var.section
        return d

    def group_by_section(self, names):
        """
        Group a list of variable in sections.

        Args:
            names: string or list of strings with ABINIT variable names.

        Return:
            Ordered dict mapping section_name to the list of variable names belonging to the section.
            The dict uses the same ordering as those in `self.sections`
        """
        d = defaultdict(list)

        for name in list_strings(names):
            try:
                sec = self.name2section[name]
                d[sec].append(name)
            except KeyError as exc:
                raise KeyError("{}\n\nIf the key is a valid Abinit variable, try to remove\n"
                               "~/.abinit/abipy/abinit_vars.pickle and rerun.".format(exc))

        return OrderedDict([(sec, d[sec]) for sec in self.sections if d[sec]])

    def apropos(self, varname):
        """Return the list of :class:`Variable` objects that are related` to the given varname"""
        var_list = []
        for v in self.values():
            if (v.text and varname in v.text or
               (v.dimensions is not None and varname in str(v.dimensions)) or
               (v.requires is not None and varname in v.requires) or
               (v.excludes is not None and varname in v.excludes)):
                var_list.append(v)

        return var_list

    def vars_with_section(self, sections):
        """
        List of :class:`Variable` associated to the given sections.
        sections can be a string or a list of strings.
        """
        sections = set(list_strings(sections))
        varlist = []
        for v in self.values():
            if v.section in sections:
                varlist.append(v)
        return varlist

    def vars_with_char(self, chars):
        """
        List of :class:`Variable` with the specified characteristic.
        chars can be a string or a list of strings.
        """
        chars = ["[[" + c + "]]" for c in list_strings(chars)]
        varlist = []
        for v in self.values():
            if v.characteristics is None: continue
            if any(c in v.characteristics for c in chars):
                varlist.append(v)

        return varlist

    def json_dumps_varnames(self):
        """JSON string with the list of variable names extracted from the database."""
        return json.dumps(list(self.keys()))


def docvar(varname):
    """Return the `Variable` object associated to this name."""
    return get_abinit_variables()[varname]


def abinit_help(varname, info=True, stream=sys.stdout):
    """
    Print the abinit documentation on the ABINIT input variable `varname`
    """
    database = get_abinit_variables()
    if isinstance(varname, Variable): varname = varname.varname
    try:
        var = database[varname]
    except KeyError:
        return stream.write("Variable %s not in database" % varname)

    html = "<h2>Default value:</h2> %s <br/><h2>Description</h2> %s" % (
        str(var.defaultval), str(var.text))
    text = html2text.html2text(html)
    if info: text += str(var.info)
    # FIXME: There are unicode chars in abinit doc (Greek symbols)
    text = text.replace("[[", "\033[1m").replace("]]", "\033[0m")

    try:
        stream.write(text)
    except UnicodeEncodeError:
        stream.write(text.encode('ascii', 'ignore'))
    stream.write("\n")


def repr_html_from_abinit_string(text):
    """
    Given a string `text` with an Abinit input file, replace all variables
    with HTML links pointing to the official documentation. Return new string.
    """
    var_database = get_abinit_variables()

    # https://stackoverflow.com/questions/6116978/python-replace-multiple-strings
    # define desired replacements here e.g. rep = {"condition1": "", "condition2": "text"}
    # ordered dict and sort by length is needed because variable names can overlap e.g. kpt, kptopt pair
    import re
    rep = {vname: var.html_link(label=vname) for vname, var in var_database.items()}
    rep = OrderedDict([(re.escape(k), rep[k]) for k in sorted(rep.keys(), key=lambda n: len(n), reverse=True)])

    # Use these three lines to do the replacement
    pattern = re.compile("|".join(rep.keys()))
    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], text)
    return text.replace("\n", "<br>")
