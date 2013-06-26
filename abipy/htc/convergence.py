from __future__ import print_function, division
import sys
import collections
import warnings
import numpy as np

from os.path import join as pj
from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS

from abipy.profile import abipy_env
from pymatgen.io.abinitio.utils import find_file
from . import Launcher, MassLauncher, InputFile

__all__ = [
    "ConvMassLauncher",
]

# The user will have to put the correct data in input_data before submitting/running the calculations.
# We don't support calculations done with the get* variables when the output file is produced by another run
# since this option causes a dependency among the runs that are difficult to manage
# (the run command now supports the -j option so that
# small calculations can be done with threads without submitting the jobs. Submitting job that are inter-dependent
# is out of the question. The use of get* inside the run is fine, though)
# The masslauncher, indeed, is designed to deal with a large number of independent calculations in which
# the input dataset is either generated automatically (e.g. input files) or is already available
# (e.g. WFK file for GW or Phonons or DEN file for NSCF runs)


##########################################################################################

class Parameter(object):
    "Simple container storing the name of the parameter and its values"
    def __init__(self, name, values):
        self.name = name
        self.values = values

    def __str__(self): return " %s = %s " % (self.name, self.values)
    def __len__(self): return len(self.values)
    def __iter__(self):
        for v in self.values: yield v


class ParameterSpace(object):
    """
    Store the set of parameteres used to generate the input files of the MassLauncher.

    Usage example:

        conv_space = ParameterSpace(ecut = range(10,40,10), tsmear = [0.01, 0.02])

    :note: Does not support more than 2 parameters.
    """
    def __init__(self, *parameters, **kwargs):
        self.parameters = list(parameters)

        for (key, values) in kwargs.items():
            self.parameters.append( Parameter(key, tuple(values)) )

        self.parameters = tuple(self.parameters)
        self.shape = tuple( [len(par) for par in self.parameters] )
        self.size = 1
        for d in self.shape: self.size *= d

    def __len__(self): return self.size

    def __iter__(self):
        for d in self.iterate(): yield d

    @property
    def parnames(self):
        return [p.name for p in self.parameters]

    @property
    def dim(self): return len(self.shape)

    def varname_index(self, varname):
        return self.parnames.index(varname)

    def flat(self):
        return [d for d in self.iterate()]

    def flatvalues(self):
        if self.dim == 1:
            return self.parameters[0].values
        else:
            return [d.values() for d in self.flat()]

    def iterate(self):
        """
        Iterate over the parameter space.

        Returns:
            Dictionary with mapping parameter_name -> parameter_value
        """
        if self.dim == 1:
            par0 = self.parameters[0]
            for v0 in par0.values:
                yield {par0.name : v0}

        elif self.dim == 2:
            par0 = self.parameters[0]
            par1 = self.parameters[1]
            for v0 in par0:
                for v1 in par1:
                    yield {par0.name : v0, par1.name : v1}
        else:
            raise NotImplementedError("dimension %s not supported " % self.dim)

    def slices_1d(self, data, along_varname=None):
        if self.dim == 1:
            return self.parameters[0].values, np.reshape(data,(1,len(data))), {}

        elif self.dim == 2:
            data = np.asarray(data)
            data.shape = self.shape
            axis = 0
            if along_varname is not None: axis = self.varname_index(along_varname)
            if axis == 0:
                slices = [col for col in data.T]
                xx     = self.parameters[0].values
                consts = self.parameters[1]
            else:
                slices = [row for row in data]
                xx     = self.parameters[1].values
                consts = self.parameters[0]

            return xx, slices, consts
        else:
            raise NotImplementedError("dimension %s not supported " % self.dim)

##########################################################################################

class ConvMassLauncherArgParser(MassLauncher.ArgParser):
    """
    Argument parser for ConvMassLauncher. Extend MassLauncher.ArgParser adding the options:

        plot
        extract?
    """

    def __init__(self, *args, **kwargs):
        MassLauncher.ArgParser.__init__(self, *args, **kwargs)

        for cmd in self.commands: self.unregister_subparser(cmd)

        p_plot = self.subparsers.add_parser('plot', help='Plot data')

        p_plot.add_argument('varnames', nargs='*', default = ["etotal"],
                            help='Name of the (scalar) variables to plot (default etotal)')

        p_plot.add_argument('-x', '--xvarname', default=SUPPRESS,
                            help='Name of the (scalar) variable displayed along the x-axis (used when two parameters are used)')

        self.register_subparser("plot", p_plot)

class ConvMassLauncher(object):
# TODO
#class ConvMassLauncher(MassLauncher):
# Here I would like to inherit from MassLauncher to extend the object
# Unfortunately the definition of__getitem__ and __setitem__
# in MassLauncher complicates the inheritance as one should replace these methods
# in the super class. However this might break the methods of the subclass
# For the moment, I'm forced to wrap MassLauncher but this approach is not ideal
# since this will force me to write a lot of boilerplate code just the expose the
# methods and the attributes of the superclass.
    """
    Usage example::

        from abipy.htc import ParameterSpace, ConvMassLauncher
        paramspace = ParameterSpace( ecut = range(10,30,10), tsmear = [0.01, 0.02] )

        # For a converge study wrt one parameter e.g. ecut one can use:
        #paramspace = ParameterSpace( ecut = range(10,30,10))

        convergence = ConvMassLauncher("atom.tmpl","14si.pspnc", paramspace)
        convergence.execute()
    """

    ArgParser = ConvMassLauncherArgParser
    argparser = ConvMassLauncherArgParser()

    Error = Launcher.Error

    def __init__(self, name, template, pseudo, parspace, **kwargs):

        # Get arguments for the MassLaunchers either from kwargs or from the abipy environment.
        pseudodir   = abipy_env.get_uservar("pseudodir", kwargs)
        jobtype     = abipy_env.get_uservar("jobtype", kwargs)
        self.bindir = abipy_env.get_uservar("bindir", kwargs)

        self.parspace = parspace
        ncalcs = len(self.parspace)

        #MassLauncher.__init__(self,
        #                      ncalcs, name,
        #                      jobtype=jobtype,
        #                      pseudodir=pseudodir,
        #                      pseudos=[pseudo,],
        #                      executable=executable)

        self.mlaunchers = mlaunchers = MassLauncher(ncalcs, name,
                                                    jobtype=jobtype,
                                                    pseudodir=pseudodir,
                                                    pseudos=[pseudo,],
                                                    executable=pj(self.bindir, "abinit")
                                                    )

        # TODO: pseudo can be read from the template file.
        if isinstance(template, InputFile):
            raise NotImplementedError("")

        elif isinstance(template, collections.Mapping):
            # Initialize input files from dictionary
            mlaunchers.set_variables(variables=template)

        else:
            # Assume string with filename.
            mlaunchers.read(template)

        for (launcher, params) in zip(mlaunchers, self.parspace.flat()):
            launcher.set_variables(params)

    def execute(self, *args, **kwargs):

        if not args:
            args = sys.argv[1:]

        options, args, kwargs = self.argparser.myparse_args(args)

        if self.argparser.can_handle(options.command):
            # This is my command e.g. plot --> Execute my code.
            getattr(self, options.command)(*args, **kwargs)

        else:
            # Delegate the superclass.
            self.mlaunchers.execute(*args, **kwargs)

    def run(self, *args, **kwargs):
        self.mlaunchers.run(*args, **kwargs)

    def make(self, *args, **kwargs):
        self.mlaunchers.make(*args, **kwargs)

    def destroy(self, *args, **kwargs):
        self.mlaunchers.destroy(*args, **kwargs)

    def clean(self, *args, **kwargs):
        self.mlaunchers.clean(*args, **kwargs)

    def read_gsr_values(self, varname):

        # TODO: we should define a hierarchy of exceptions for htc
        status = self.mlaunchers.report()

        if any([s != 'Completed' for s in status]):
            raise self.Error("Calculation not completed, status = %s" % status)

        #output_files = [ l.last_output() for l in self.mlaunchers]
        output_files = self.mlaunchers.last_output()
        #print output_files

        ml_odatfiles = self.mlaunchers.odat_files()
        #print 'ml_odatfiles', ml_odatfiles

        # Here we read etotal from ab_out, very fragile!
        # moreover only etotal is supported!
        #etotals = []
        #for ofname in output_files:
        #  with open(ofname, "r") as fh:
        #    for line in fh:
        #      if "etotal" in line:
        #        etotals.append( float(line.split()[1]) )
        #        break
        #print "etotals", etotals

        # Here we read varname from the netcdf file GSR.
        # Elegant, flexible and easy to maintain.
        gs_etotals = []
        for launcher in self.mlaunchers:
            # Find the correct output file
            out_files = launcher.odat_files()
            gsfname = find_file(out_files, "GSR")
            if gsfname is None:
                raise RuntimeError("Cannot find GSR file among " % out_files)

            # Open the GSR file and read the total energy
            with GSR_Reader(gsfname) as gsfile:
                gs_etotals.append( gsfile.get_value(varname) )

        #print "gs_etotals",gs_etotals

        return gs_etotals

    def plot(self, *args, **kwargs):
        """
        Plot the value of varname as function of ecut.
        """
        import matplotlib.pyplot as plt

        varnames = kwargs.get("varnames")
        assert len(varnames) == 1
        varname = varnames[0]

        try:
            data = self.read_gsr_values(varname=varname)
        except self.Error as e:
            warnings.warn(str(e))
            return

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        #xx = self.parspace.flatvalues()
        parnames = self.parspace.parnames

        #print parnames
        #print self.parspace.shape
        #print self.parspace.dim
        #print self.parspace.size
        #print "xx",xx
        for par, datum in zip(self.parspace.flat(), data):
            print("par, datum", par, datum)

        xname = kwargs.get("xvarname", self.parspace.parnames[0])

        #TODO recheck the API of this routine
        xx, slices, constpar = self.parspace.slices_1d(data, along_varname=xname)
        #print xx
        #print yy
        print("data",data)
        print("slices",slices)

        lines = []
        for yy in slices:
            line, = ax.plot(xx, yy, "-->", linewidth=3.0, markersize=10)
            lines.append(line)

        if constpar: # 2D case.
            legend_entries = ["%s = %s" % (constpar.name, v) for v in constpar ]
            ax.legend(lines, legend_entries, 'upper right', shadow=True)
        else:
            ax.legend(lines, [varname,], 'upper right', shadow=True)

        # Set xticks and labels.
        ax.grid(True)
        ax.set_xlabel(xname)
        ax.set_ylabel(varname)
        ax.set_xticks(xx)

        plt.show()

##########################################################################################

