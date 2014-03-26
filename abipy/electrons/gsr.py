"""GSR file."""
from __future__ import print_function, division

import collections
import numpy as np

from abipy.iotools import AbinitNcFile, Has_Structure, Has_ElectronBands
from abipy.tools import list_strings
from .ebands import ElectronsReader

__all__ = [
    "GSR_File",
    "GSR_Plotter",
]

class GSR_File(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    File containing the results of a Ground-state calculation.

    Usage example:
                                                                  
    .. code-block:: python
        
        gsr = GSR_File("foo_GSR.nc")
        gsr.ebands.plot()
    """
    def __init__(self, filepath):
        super(GSR_File, self).__init__(filepath)

        with GSR_Reader(filepath) as r:
            # Initialize the electron bands from file
            self._ebands = r.read_ebands()

            #self.structure.set_forces(r.read_forces())
            #self.structure.set_stress(r.read_stress())

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    @property
    def ebands(self):
        """`ElectronBands` object."""
        return self._ebands

    @property
    def structure(self):
        """`Structure` object."""
        return self.ebands.structure

    @property
    def kpoints(self):
        """Iterable with the Kpoints."""
        return self.ebands.kpoints

    @property
    def nsppol(self):
        return self.ebands.nsppol

    @property
    def mband(self):
        return self.ebands.mband


class GSR_Reader(ElectronsReader):
    """
    This object reads the results stored in the _GSR (Ground-State Results)
    file produced by ABINIT. It provides helper function to access the most
    important quantities.
    """
    #def __init__(self, filepath):
    #    """Initialize the object from a filepath."""
    #    super(GSR_Reader, self).__init__(filepath)

    #def read_forces(self):

    #def read_stress(self):

    #def read_energies(self):
    #    return AttrDict()


class GSR_Plotter(collections.Iterable):
    """
    This object receives a list of `GSR_File` objects and provides
    methods to inspect/analyze the results (useful for convergence studies)

    Usage example:
                                                                  
    .. code-block:: python
        
        plotter = GSR_Plotter()
        plotter.add_file("foo_GSR.nc")
        plotter.add_file("bar_GSR.nc")
        plotter.plot_variables("ecut", "etotal")
    """
    def __init__(self):
        self._gsr_files = collections.OrderedDict()

    def __len__(self):
        return len(self._gsr_files)

    def __iter__(self):
        return iter(self._gsr_files.values())

    def __str__(self):
        s = ""
        for gsr in self:
            s += str(gsr) + "\n"
        return s

    def add_files(self, filepaths):
        """Add a list of filenames to the plotter"""
        for filepath in list_strings(filepaths):
            self.add_file(filepath)

    def add_file(self, filepath):
        """Add a filename to the plotter"""
        from abipy.abilab import abiopen
        gsr = abiopen(filepath)
        self._gsr_files[gsr.filepath] = gsr

        # Initialize/check useful quantities.
        #
        # 1) Number of spins
        if not hasattr(self, "nsppol"): 
            self.nsppol = gsr.nsppol
        assert self.nsppol == gsr.nsppol

    @property
    def filepaths(self):
        """Filepaths of the GSR files."""
        return self._gsr_files.keys()

    @property
    def param_name(self):
        """The name of the parameter whose value is checked for convergence."""
        try: 
            return self._param_name

        except AttributeError:
            self.set_param_name(param_name=None)
            return self.param_name

    def _get_param_list(self):
        """Return a dictionary with the values of the parameters extracted from the GSR files."""
        param_list = collections.defaultdict(list)
                                                               
        for gsr in self:
            for pname in gsr.params.keys():
                param_list[pname].append(gsr.params[pname])

        return param_list

    def set_param_name(self, param_name):
        """
        Set the name of the parameter whose value is checked for convergence.
        if param_name is None, we try to find its name by inspecting 
        the values in the gsr.params dictionaries.
        """
        self._param_name = param_name

    def prepare_plot(self):
        """
        This method must be called before plotting data.
        It tries to figure the name of paramenter we are converging 
        by looking at the set of parameters used to compute the different GSR files.
        """
        param_list = self._get_param_list()

        param_name, problem = None, False
        for key, value_list in param_list.items():
            if any(v != value_list[0] for v in value_list):
                if param_name is None:
                    param_name = key
                else:
                    problem = True
                    warnings.warn("Cannot perform automatic detection of convergence parameter.\n" + 
                                  "Found multiple parameters with different values. Will use filepaths as plot labels.")

        self.set_param_name(param_name if not problem else None)

        if self.param_name is None:
            # Could not figure the name of the parameter.
            xvalues = range(len(self))
        else:
            xvalues = param_list[self.param_name]
                                                  
            # Sort xvalues and rearrange the files.
            items = sorted([iv for iv in enumerate(xvalues)], key=lambda item: item[1])
            indices = [item[0] for item in items]
                                                                                             
            files = self._gsr_files.values()
                                                                                             
            newd = collections.OrderedDict()
            for i in indices:
                gsr = files[i]
                newd[gsr.filepath] = gsr
                                                                                             
            self._gsr_files = newd

            # Use sorted xvalues for the plot.
            param_list = self._get_param_list()
            xvalues = param_list[self.param_name]

        self.set_xvalues(xvalues)

    @property
    def xvalues(self):
        """The values used for the X-axis."""
        return self._xvalues 

    def set_xvalues(self, xvalues):
        """xvalues setter."""
        assert len(xvalues) == len(self)
        self._xvalues = xvalues

    def decorate_ax(self, ax, **kwargs):
        ax.grid(True)
        if self.param_name is not None:
            ax.set_xlabel(self.param_name)

        ax.set_ylabel('Energy [eV]')
        ax.legend(loc="best")

        title = kwargs.pop("title", None)
        if title is not None:
            ax.set_title(title)
                                                                                 
        # Set ticks and labels. 
        if self.param_name is None:
            # Could not figure the name of the parameter ==> Use the basename of the files
            ticks, labels = range(len(self)), [f.basename for f in self]
        else:
            ticks, labels = self.xvalues, [f.params[self.param_name] for f in self]

        ax.set_xticks(ticks, minor=False)
        ax.set_xticklabels(labels, fontdict=None, minor=False)

    def plot_variables(self, varname_x, varname_y, hspan=None, **kwargs):
        """
        Ex:
            plot_variables("ecut", "etotal")
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        # Read the value of varname from the files.
        xx, yy = [], []
        for filepath in self.filepaths:
            with GSR_Reader(filepath) as r:
                xx.append(r.read_value(varname_x))
                yy.append(r.read_value(varname_y))

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        ax.plot(xx, yy, "o-", **kwargs)

        if hspan is not None:
            last = yy[-1]
            ax.axhspan(last-hspan, last+hspan, facecolor='0.5', alpha=0.5)

        if title is not None:
            fig.suptitle(title)
                                 
        if show:
            plt.show()
                                 
        if savefig is not None:
            fig.savefig(savefig)
                                 
        return fig

    #def eos_fit(eos_name='murnaghan'):
    #    """Fit E(V)"""
    #    # Read the value of varname from the files.
    #    volumes, energies = [], []
    #    for filepath in self.filepaths:
    #        with GSR_Reader(filepath) as r:
    #            volumes.append(r.read_value(varname_x))
    #            energies.append(r.read_value("etotal"))

    #    # For the list of available models, see EOS.MODELS
    #    from pymatgen.io.abinitio.eos import EOS
    #    eos = EOS(eos_name=eos_name)

    #    # Note that eos.fit expects lengths in Angstrom, energies are in eV.
    #    # To specify different units use len_units and ene_units 
    #    fit = eos.fit(volumes, energies, vol_unit="bohr**3", ene_unit="Ha")
    #    #fit.plot()
    #    return fit
