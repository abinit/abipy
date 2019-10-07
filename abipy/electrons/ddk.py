# coding: utf-8
"""DDK.nc file."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.core.units as units

from monty.string import marquee # is_string, list_strings,
from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from abipy.core.func1d import Function1D
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools import gaussian, duck
from abipy.electrons.ebands import ElectronsReader, ElectronDos
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims


#class DdksAnalyzer(object):
#    """
#    This object received three ddks files with derivatives along the 3 reduced directions.
#    """
#
#    def __init__(self, ddk_paths):
#        # Open files here. The order of input ddk_paths does not matter as long as we can fill all ddks.
#        self.ddks = ddks = 3 * [None]
#        for path in ddk_paths:
#            ddk = DdkFile(path)
#            ddks[ddk.idir - 1] = ddk
#        if any(ddk is None for ddk in ddks):
#            raise ValueError("Cannot find 3 DDK files with different idir.")
#
#        # Consistency check
#        errors = []
#        eapp = errors.append
#        if any(ddk.structure != ddks[0].structure for ddk in ddks[1:]):
#            eapp("Structures in DDK files do not agree with each other.")
#        if any(ddk.kptopt != ddks[0].kptopt for ddk in ddks[1:]):
#            eapp("Found different values of kptopt.")
#        if any(ddk.ebands.kpoints != ddks[0].ebands.kpoints for ddk in ddks[1:]):
#            eapp("Found different list of k-kpoints.")
#        for aname in ("nsppol", "nspden", "nspinor"):
#            if any(getattr(ddk.ebands, aname) != getattr(ddks[0].ebands, aname) for ddk in ddks[1:]):
#                eapp("Found different value of %s" % aname)
#        if errors:
#            raise ValueError("\n".join(errors))
#
#        # Get useful dimensions.
#        eb0 = ddks[0].ebands
#        self.nsppol, self.nspinor, self.nspden = eb0.nsppol, eb0.nspinor, eb0.nspden
#        self.nband = eb0.nband
#        self.kpoints = eb0.kpoints
#        self.nkpt = len(eb0.kpoints)
#        self.weights = eb0.kpoints.weights
#        self.eigens = eb0.eigens
#
#    @lazy_property
#    def v_skb(self):
#        v_skb = np.empty((self.nsppol, self.nkpt, self.nband, 3))
#        for i, ddk in enumerate(self.ddks):
#            v_skb[:, :, :, i] = ddk.reader.read_ddk_diagonal()
#        return v_skb
#
#    @add_fig_kwargs
#    def plot_ebands_with_doses(self, ebands_kpath, doses, ylims=None, **kwargs):
#        """
#        Plot band structure and doses
#
#        Args:
#            ebands_kpath:
#            doses:
#            ylims: Set the data limits for the x-axis in eV. Accept tuple e.g. ``(left, right)``
#                or scalar e.g. ``left``. If left (right) is None, default values are used
#
#        Return: |matplotlib-Figure|
#        """
#        # Build grid plot.
#        import matplotlib.pyplot as plt
#        from matplotlib.gridspec import GridSpec
#        fig = plt.figure()
#        ncols = 3
#        width_ratios = [2, 0.2, 0.2]
#        gspec = GridSpec(1, ncols, width_ratios=width_ratios)
#        gspec.update(wspace=0.05)
#
#        ax_ebands = plt.subplot(gspec[0])
#        ax_doses = []
#        for i in range(2):
#            ax = plt.subplot(gspec[i + 1], sharey=ax_ebands)
#            ax_doses.append(ax)
#            ax.grid(True)
#            set_axlims(ax, ylims, "y")
#
#        # Plot electron bands.
#        ebands_kpath.plot(ax=ax_ebands, ylims=ylims, show=False)
#
#        # Plot DOSes.
#        #doses.edos.plot
#        #vdos.edos.plot
#        #ax.set_ylabel("")
#
#        return fig


class DdkFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File containing the DDK matrix elements for a single `idir` direction

    Usage example:

    .. code-block:: python

        with DdkFile("out_DDK.nc") as ddk:
            ddk.ebands.plot()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(DdkFile, self).__init__(filepath)
        self.reader = DdkReader(filepath)

        # Get info on perturbation and k-point sampling.
        #self.kptopt = self.reader.read_value("kptopt")
        #assert self.kptopt == 2
        #pertcase = self.reader.read_value("pertcase")
        #self.idir = ((pertcase - 1) % 3) + 1
        #self.ipert = (pertcase - self.idir) // 3 + 1

    def close(self):
        """Close the file."""
        self.reader.close()

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def edos(self):
        """
        |ElectronDos| object with the electron DOS computed by Abinit from the input WFK file.
        """
        # NB: This is optional
        mesh = self.reader.read_value("edos_mesh") * abu.Ha_eV
        # nctkarr_t("edos_dos", "dp", "edos_nw, nsppol_plus1"), &
        # dos(nw,0:nsppol) Total DOS, spin up and spin down component.
        spin_dos = self.reader.read_value("edos_dos") / abu.Ha_eV
        nelect = self.ebands.nelect
        fermie = self.ebands.fermie

        return ElectronDos(mesh, spin_dos[1:], nelect, fermie=fermie)

    #@lazy_property
    #def jdos(self)
        # NB: This is optional

    #def make_fermisurfer_dir(self, workdir)
    #    eb3d = self.ebands.get_ebands3d()
    #    eb3d.make_fermisurfer_dir(workdir)

    @lazy_property
    def params(self):
        """OrderedDict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.structure.plot(show=False)
        #yield self.ebands.plot(show=False)

    def write_notebook(self, nbpath=None, title=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        nb.cells.extend([
            nbv.new_code_cell("ddk = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ddk)"),
            nbv.new_code_cell("ddk.ebands.plot();"),
            #nbv.new_code_cell("ddk.ebands.kpoints.plot();"),
            nbv.new_code_cell("""\
if ddk.ebands.kpoints.is_ibz:
    ddk.ebands.get_edos().plot();"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class DdkReader(ElectronsReader):
    """
    This object reads the results stored in the DDK.nc file
    It provides helper function to access the most important quantities.
    """
    def __init__(self, filepath):
        super(DdkReader, self).__init__(filepath)
        nband_sk = self.read_nband_sk()
        if np.any(nband_sk != nband_sk[0, 0]):
            raise NotImplementedError("Found different number of bands per k-point, spin.\nnband_sk: %s\n" %
                    str(nband_sk))
        self.mband = nband_sk[0, 0]

    #def read_ddk_diagonal(self):
    #    """
    #    Read the group velocities i.e the diagonal matrix elements.
    #    Return (nsppol, nkpt) |numpy-array| of real numbers.
    #    """
    #    var = self.read_variable("h1_matrix")
    #    vels = np.diagonal(var[:, :, :, :, :], axis1=2, axis=3)
    #    # Cartesian? Ha --> eV?
    #    return np.real(vels).copy() * (units.Ha_to_eV / units.bohr_to_ang)

    #def read_ddk_skbb(self):
    #    return self.read_value("h1_matrix", cplx_mode="cplx") * (units.Ha_to_eV / units.bohr_to_ang)
