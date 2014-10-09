import os
import pwd
import imp

from os.path import exists, join as pj
from collections import OrderedDict
from warnings import warn

__all__ = [
    "abipy_env",
]


class FrozenDict(dict):
    "A dictionary that does not permit to redefine its keys"

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        if key in self:
          raise Exception("Cannot overwrite existent key: %s" % str(key))
        dict.__setitem__(self, key, val)

    def update(self, *args, **kwargs):
        #print 'in frozen update', args, kwargs
        for (k, v) in dict(*args, **kwargs).items():
            self[k] = v

class AbipyEnvironment(FrozenDict):
   """
   Frozen dictionary with the abipy environment variables read from the configuration file.

   Values can be accessed either with the standard syntax for dictionaries:

   >>> d = AbipyEnvironment({"foo" : "bar"})
   >>> d["foo"] 
   'bar'

   or with the attribute syntax:

   >>> d.foo
   'bar'
   """

   def __init__(self, *args, **kwargs):
      FrozenDict.__init__(self, *args, **kwargs)

   def __getattribute__(self, name):
     try:
        return self[name]
     except KeyError:
        return dict.__getattribute__(self, name)
                                                                   
   def __setattr__(self, name, value):
     if name in self:
        raise Exception("Cannot overwrite existent key: %s" % name)
     else:
        return dict.__settattribute__(self, name, value)

   def get_uservar(self, varname, kwargs):
      "Return the value of the variable varname. Note that kwargs has the precedence over self"
      assert varname in self
      try: 
        return kwargs[varname]
      except KeyError:
        return self[varname]

   def get_default(self, varname, kwargs, default=None):
      """
      Return the value of the variable varname, relying first on kwargs,
      then on self, and finaly on default
      """
      if varname in kwargs:
         return kwargs[varname]
      elif varname in self:
         return self[varname]
      else:
         return default



class EnvVar(object):
    "An Abipy environment variable"

    def __init__(self, name, type, default, allowed_values, help):
        self.name            = name
        self.type            = type 
        self.default         = default
        self.allowed_values  = allowed_values
        self.help            = help
        if self.help[-1] != "\n": self.help += "\n"

    def totemplate(self):
        "Return the string used to write the abipyrc template file"
        docstr = "# " + self.help.replace("\n","\n#")[:-1]

        if self.allowed_values: 
            docstr += "# Allowed values: %s\n" % str(self.allowed_values)
        
        if self.type == "string":
            assignment = "%s = '%s' \n#\n" % (self.name, self.default)
        else:
            assignment = "%s = %s \n#\n" % (self.name, self.default)
                                                                                            
        if self.default is None:    
            assignment = "#" + assignment

        return "#\n".join([docstr , assignment])


class AbipyrcParser(object):
    """
    Responsible for parsing the configuaration file, managing the options and writing the template file
    """

    _type_converter = {
       "string"   : str, 
       "int"      : int,
       "float"    : float,
       "complex"  : complex,
    }

    def __init__(self):
        self._evars = OrderedDict()

    def add_option(self, opt_name, **kwargs): #type, default, allowed_values, help)
        "Register a new configuration option"

        type = kwargs.pop("type","string")
        default = kwargs.pop("default", None)
        allowed_values = kwargs.pop("allowed_values", [])

        convert = self._type_converter[type]
        allowed_values = [convert(v) for v in allowed_values]

        try:
            help  = kwargs.pop("help")
        except KeyError:
            raise KeyError("Help string must be provided")

        if kwargs:
            raise ValueError("Unknown keywords %s" % kwargs.keys())

        if opt_name not in self._evars:
            self._evars[opt_name] = EnvVar(opt_name, type, default, allowed_values, help)
        else:
            raise ValueError("Option %s is already defined" % opt_name)

    def get_converter(self, opt_name):
        """
        Return the converter associated to the option name opt_name. None if opt_name is not a valid option.
        """
        try:
            type = self._evars[opt_name].type
            return self._type_converter[type]
        except:
            return None

    #def get_default(self, opt_name)
    #    """
    #    Return the default value associated to the option name opt_name.
    #    None if opt_name is not a valid option.
    #    """
    #    try:
    #        return self._evars[opt_name].default
    #    except:
    #        return None

    def known_options(self):
       return self._evars.keys()

    def optval_is_allowed(self, opt_name, opt_value):
        "Check whether the value opt_value of option opt_name is valid"
        allowed_values = self._evars[opt_name].allowed_values
        isok = True
        if allowed_values: isok = opt_value in allowed_values
        #if not isok: print "Not allowed: ",opt_name,opt_value,type(opt_value),allowed_values,type(allowed_values)
        return isok

    def parse(self, filename):
        """
        Parse the configuration file
        :return: An instance of :class:`AbipyEnvironment`
        """
        #print 'Reading configuration parameters from file %s' % filename

        env = AbipyEnvironment()

        if not exists(filename):
            #warn("Abipy configuration file %s does not exist!" % filename)
            return env

        module_name = filename
        module = imp.load_source(module_name, filename)

        keys = [ k for k in dir(module) if not k.startswith("_") ]

        unknown, wrong_value, wrong_conversion = [], [], []

        for k in keys:
            convert = self.get_converter(k)
            if convert is None:
                unknown.append(k)
                continue

            try:
                value = convert(module.__dict__[k])
                if self.optval_is_allowed(k, value): 
                    env[k] = value
                else:
                    wrong_value.append(k)
            except:
                wrong_conversion.append((k, module.__dict__[k]))

        if unknown:
            warn("Unknown variables found in configuration file:\n %s " % "\n".join([u for u in unknown]) )

        if wrong_value:
            warn("Variables with wrong value found in configuration file:\n %s " % "\n".join([w for w in wrong_value]) )
        
        if wrong_conversion:
            msg = "Wrong conversion for\n " + "\n".join([str(tup) for tup in wrong_conversion])
            warn(msg)

        return env

    def write_template(self, fileobj):
        """
        Write the template file.
        :arg fileobj: file-like object or string.
        """

        from StringIO import StringIO
        class MyStringIO(StringIO):
            def writen(self, string, nn=1): 
                self.write(string+"\n")
                if nn>1: self.write((nn-1)*"#\n")

        stio = MyStringIO()

        stio.writen("#")
        stio.writen("# Configuration file for abipy. Automatically generated by " + __file__, nn=2)
        stio.writen("# This file sets default values. Uncomment lines to activate them.", nn=2)
        ##
        ## It can be found in
        ##           $HOME/.abitools/abitoolsrc_template
        ##
        ## You should copy it to
        ##           $HOME/.abitools/abitoolsrc
        ##
        ## Otherwise, it will be overwritten
        ## at each installation.

        for (varname, var) in self._evars.items():
            stio.write(var.totemplate())

        stio.seek(0)

        if hasattr(fileobj, "writelines"):
            fileobj.writelines([l for l in stio])
        else:
            with open(fileobj, "w") as fh: 
                fh.writelines([l for l in stio])

##########################################################################################
#
# The Abipy environment variables must be declared and described here.
#
#   default        --> Default value, None (standard value) if an explicit value must be provided
#                      by the user (in this case the line is commented in the template file.)
#   type           --> string, int, float, complex, list. Default is string
#   allowed_values --> list with the allowed values, [] if not constraint is imposed (default)
#                      Used for consistency check when reading the configuration file.
#   help           --> string with a brief description of the meaning of the variable (mandatory).
#
##########################################################################################

parser = AbipyrcParser()

parser.add_option("pseudodir", help="Directory where to look for pseudopotentials.")

parser.add_option("bindir", help= "Directory where to look for executables.")

parser.add_option("jobtype", default="sge", allowed_values=["sge", "pbs", "slurm"], help="Resource manager.")

parser.add_option("shell", default="csh", allowed_values=["bash", "csh"], help="Shell type.")

parser.add_option("mpirun", default="", help="MPI runner. E.g. 'mpiexec -n 64'"),

parser.add_option("modules", help= "List of modules that will be loaded via the command module load <mod>.\n" + 
                                   " e.g ['intel-compilers/12.0.4.191', 'FFTW/3.3' ]" )

parser.add_option("lines_before", help = "List of lines executed before abinit is called.\n" + 
                                         " e.g ['cd ${PBS_O_WORKDIR}', 'limit stacksize unlimited']" )

parser.add_option("lines_after", help = "List of lines executed after abinit is called, e.g. ['echo Job_done!']")

parser.add_option("other_lines", help="List of command lines your job submission system you would like to add.\n")

parser.add_option("visualizer", default="xcrysden", allowed_values=["xcrysden",], 
                   help="Visualization program for crystalline structures, densities ...")

##########################################################################################

# Configuration file

# Avoid using the HOME dir of root when sudo is used.
username= pwd.getpwuid(os.getuid())[0]
homepath = os.path.expanduser("~"+username+"/")
#homepath = os.getenv("HOME"),

abipy_cfg_dir = pj(homepath, ".abinit", "abipy")
abipyrc_path = pj(abipy_cfg_dir, "abipyrc")

# Create the dictionary with the environment
# Other modules should import abipyenv to have access to the abipy variables.
abipy_env = parser.parse(abipyrc_path)

def write_abipy_cfg_data():
    "Write the configuration files used by abipy"

    # Create directories (recursively)
    if not exists(abipy_cfg_dir): 
        os.makedirs(abipy_cfg_dir)

    if not exists(abipyrc_path):
        # Write the configuration file.
        parser.write_template(abipyrc_path)
    else:
        # Do not overwrite the old config file, create a file with the new template.
        parser.write_template(abipyrc_path + ".new")


if __name__ == "__main__":
    write_abipy_cfg_data()
