# coding: utf-8
"""
TRANSPORT.nc file.

Warning: This fileformat is deprecated and will be removed when Abinit 9.2 is released
"""

import numpy as np
import abipy.core.abinit_units as abu

from monty.functools import lazy_property
from monty.string import marquee
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.abio.robots import Robot


__all__ = [
    "TransportFile",
]


class TransportFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super().__init__(filepath)
        self.reader = TransportReader(filepath)

        #self.fermi = self.ebands.fermie * abu.eV_Ha
        self.tmesh = self.reader.tmesh
        #self.transport_ngkpt = self.reader.read_value("transport_ngkpt")
        #self.transport_extrael = self.reader.read_value("transport_extrael")
        #self.transport_fermie = self.reader.read_value("transport_fermie")

    @property
    def ntemp(self):
        """Number of temperatures."""
        return len(self.tmesh)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def __str__(self):
        """String representation"""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation"""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="KS Electron Bands"))
        app("")

        # Transport section.
        app(marquee("Transport calculation", mark="="))
        app("Number of temperatures: %d" % self.ntemp)
        app("Mobility:")
        app("Temperature [K]     Electrons [cm^2/Vs]     Holes [cm^2/Vs]")
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            mobility_mu_e = self.get_mobility_mu(0, itemp)
            mobility_mu_h = self.get_mobility_mu(1, itemp)
            app("%14.1lf %18.6lf %18.6lf" % (temp, mobility_mu_e, mobility_mu_h))

        return "\n".join(lines)

    @add_fig_kwargs
    def plot_edos(self, ax=None, **kwargs):
        """
        Plot the electronic density of states

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        wmesh, dos, idos = self.reader.read_dos()
        ax.plot(wmesh, dos[0], **kwargs)
        ax.grid(True)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('States/eV')

        return fig

    @add_fig_kwargs
    def plot_vvtau_dos(self, component='xx', ax=None, colormap='jet', fontsize=8, **kwargs):
        """
        Plot velocity * lifetime density of states.

        Args:
            component: Component to plot: "xx", "yy" "xy" ...
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib colormap.
            fontsize (int): fontsize for titles and legend

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            wmesh, vvdos = self.reader.read_vvdos_tau(itemp, component=component)
            ax.plot(wmesh, vvdos, c=cmap(itemp / self.ntemp), label='T = %dK' % temp)

        ax.grid(True)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('VVDOS')
        ax.set_yscale('log')
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_mobility(self, component='xx', ax=None, colormap='jet', fontsize=8, **kwargs):
        """
        Read the mobility from the netcdf file and plot it

        Args:
            component: Component to plot: "xx", "yy" "xy" ...
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib colormap.
            fontsize (int): fontsize for titles and legend

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            wmesh, mu = self.reader.read_mobility(0, itemp, component,0)
            ax.plot(wmesh, mu, c=cmap(itemp / self.ntemp), label='T = %dK' % temp)

        ax.grid(True)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel(r'mobility $\mu(\epsilon_F)$ [cm$^2$/Vs]')
        ax.set_yscale('log')
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def get_mobility_mu(self, eh, itemp, component='xx', ef=None, spin=0):
        """
        Get the value of the mobility at a chemical potential Ef

        Args:
            eh:
            itemp: Index of the temperature.
            component: Component to plot: "xx", "yy" "xy" ...
            ef: Value of the doping in eV. The default None uses the chemical potential at the temperature
            spin: Spin index.
        """
        from scipy import interpolate
        if ef is None: ef = self.reader.read_value('transport_mu_e')[itemp]
        wmesh, mobility = self.reader.read_mobility(eh, itemp, component, spin)
        f = interpolate.interp1d(wmesh, mobility)

        return f(ef)

    #@add_fig_kwargs
    #def plot_onsanger(self, nn=0, ax=None, **kwargs):
    #    """
    #    Plot Onsanger

    #    Args:
    #        ax: |matplotlib-Axes| or None if a new figure should be created.

    #    Return: |matplotlib-Figure|
    #    """
    #    ax, fig, plt = get_ax_fig_plt(ax=ax)
    #    return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        Return figures plotting the transport data
        """
        yield self.plot_edos(show=False, title="Density of states")
        yield self.plot_vvtau_dos(show=False, title="VVDOS")
        yield self.plot_mobility(show=False, title="Mobility")

    def close(self):
        """Close the file."""
        self.reader.close()

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.plot_edos();"),
            nbv.new_code_cell("ncfile.plot_vvtau_dos();"),
            nbv.new_code_cell("ncfile.plot_mobility();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class TransportReader(ElectronsReader):
    """
    This class reads the results stored in the TRANSPORT.nc file
    It provides helper function to access the most important quantities.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        super().__init__(filepath)

        ktmesh = self.read_value("kTmesh")
        self.tmesh = ktmesh / abu.kb_HaK
        self.nsppol = self.read_dimvalue('nsppol')

    def read_vvdos(self, component='xx', spin=1):
        """
        Read the group velocity density of states
        The vvdos_vals array has 3 dimensions (3, 3, nsppolplus1, nw)

          1. 3x3 components of the tensor
          2. the spin polarization + 1 for the sum
          3. the number of frequencies
        """
        # nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three"), &
        i, j = abu.s2itup(component)
        wmesh = self.read_variable("vvdos_mesh")[:] * abu.Ha_eV
        vals = self.read_variable("vvdos_vals")
        vvdos = vals[i,j,spin,:]

        return wmesh, vvdos

    def read_vvdos_tau(self, itemp, component='xx', spin=1):
        """
        Read the group velocity density of states
        The vvdos_vals array has 3 dimensions (3, 3, nsppolplus1, nw)

          1. 3x3 components of the tensor
          2. the spin polarization + 1 for the sum
          3. the number of frequencies
        """
        # nctkarr_t('vvdos_tau', "dp", "edos_nw, nsppol_plus1, three, three, ntemp"), &
        i, j = abu.s2itup(component)
        wmesh = self.read_variable("vvdos_mesh")[:] * abu.Ha_eV
        vals = self.read_variable("vvdos_tau")
        vvdos_tau = vals[itemp,i,j,spin,:] / (2 * abu.Ha_s)

        return wmesh, vvdos_tau

    def read_dos(self):
        """
        Read the density of states (in eV units)
        """
        # Total DOS, spin up and spin down component.
        # nctkarr_t("edos_dos", "dp", "edos_nw, nsppol_plus1")
        wmesh = self.read_value("edos_mesh") * abu.Ha_to_eV
        dos = self.read_value("edos_dos") / abu.Ha_to_eV
        idos = self.read_value("edos_idos")
        return wmesh, dos, idos

    def read_onsager(self, itemp):
        """
        Read the Onsager coefficients computed in the transport driver in Abinit
        """
        # nctkarr_t('L0', "dp", "edos_nw, nsppol, three, three, ntemp"), &
        L0 = np.moveaxis(self.read_variable("L0")[itemp,:], [0,1,2,3], [3,2,0,1])
        L1 = np.moveaxis(self.read_variable("L1")[itemp,:], [0,1,2,3], [3,2,0,1])
        L2 = np.moveaxis(self.read_variable("L2")[itemp,:], [0,1,2,3], [3,2,0,1])

        return L0, L1, L2

    def read_transport(self, itemp):
        # nctkarr_t('sigma',   "dp", "edos_nw, nsppol, three, three, ntemp"), &
        sigma = np.moveaxis(self.read_variable("sigma")[itemp,:],     [0,1,2,3], [3,2,0,1])
        kappa = np.moveaxis(self.read_variable("kappa")[itemp,:],     [0,1,2,3], [3,2,0,1])
        seebeck = np.moveaxis(self.read_variable("seebeck")[itemp,:], [0,1,2,3], [3,2,0,1])
        pi = np.moveaxis(self.read_variable("pi")[itemp,:],           [0,1,2,3], [3,2,0,1])

        return sigma, kappa, seebeck, pi

    def read_mobility(self, eh, itemp, component, spin):
        """
        Read mobility from the TRANSPORT.nc file
        The mobility is computed separately for electrons and holes.
        """
        # nctkarr_t('mobility',"dp", "edos_nw, nsppol, three, three, ntemp, two"), &
        i, j = abu.s2itup(component)
        wvals = self.read_variable("vvdos_mesh")
        mobility = self.read_variable("mobility")[eh,itemp,i,j,spin,:]

        return wvals, mobility

    #def read_evk_diagonal(self):
    #    """
    #    Read the group velocities i.e the diagonal matrix elements.
    #    Return (nsppol, nkpt) |numpy-array| of real numbers.
    #    """
    #    vels = self.read_variable("vred_diagonal")
    #    # Cartesian? Ha --> eV?
    #    return vels * (abu.Ha_to_eV / abu.Bohr_Ang)


class TransportRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple TRANSPORT.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: TransportRobot
    """

    EXT = "TRANSPORT"

    @add_fig_kwargs
    def plot_mobility_conv(self, eh=0, component='xx', itemp=0, spin=0, fontsize=14, ax=None, **kwargs):
        """
        Plot the convergence of the mobility obtained in a list of files

        Args:
            eh: 0 for electrons, 1 for holes
            component: Component to plot ('xx', 'xy', ...)
            itemp: Index of the temperature.
            spin: Spin index.
            fontsize: fontsize for legends and titles
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        i, j = abu.s2itup(component)

        res = []
        for ncfile in self.abifiles:
            kptrlatt = ncfile.reader.read_value('kptrlatt')
            kptrlattx = kptrlatt[0, 0]
            kptrlatty = kptrlatt[1, 1]
            kptrlattz = kptrlatt[2, 2]
            #nkpt = ncfile.nkpt
            mobility = ncfile.reader.read_value('mobility_mu')[itemp][i,j][spin][eh]
            res.append([kptrlattx, mobility])

        res.sort(key=lambda t: t[0])
        res = np.array(res)

        size = 14
        if eh == 0:
            ax.set_ylabel(r'Electron mobility (cm$^2$/(V$\cdot$s))', size=size)
        elif eh == 1:
            ax.set_ylabel(r'Hole mobility (cm$^2$/(V$\cdot$s))', size=size)
        else:
            raise ValueError("Invalid value for eh argument: %s" % eh)

        from fractions import Fraction
        ratio1 = Fraction(kptrlatty, kptrlattx)
        ratio2 = Fraction(kptrlattz, kptrlattx)
        text1 = '' if ratio1.numerator == ratio1.denominator else \
                r'$\frac{{{0}}}{{{1}}}$'.format(ratio1.numerator, ratio1.denominator)
        text2 = '' if ratio2.numerator == ratio2.denominator else \
                r'$\frac{{{0}}}{{{1}}}$'.format(ratio2.numerator, ratio2.denominator)

        ax.set_xlabel(r'Homogeneous $N_k \times$ ' + text1 + r'$N_k \times$ ' + text2 + r'$N_k$ $\mathbf{k}$-point grid',
                      size=size)

        ax.plot(res[:,0], res[:,1], **kwargs)

        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        #yield self.plot_lattice_convergence(show=False)
        #yield self.plot_gsr_convergence(show=False)
        #for fig in self.get_ebands_plotter().yield_figs(): yield fig
        #self.plot_mobility_conv(eh=0, component='xx', itemp=0, spin=0, fontsize=14, ax=None, **kwargs):

    #def get_panel(self):
    #    """
    #    Build panel with widgets to interact with the |GsrRobot| either in a notebook or in panel app.
    #    """
    #    from abipy.panels.transportfile import TransportRobotPanel
    #    return TransportRobotPanel(self).get_panel()

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GsrRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        # Mixins
        #nb.cells.extend(self.get_baserobot_code_cells())
        #nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


if __name__ == "__main__":
    import sys
    robot = TransportRobot.from_files(sys.argv[1:])
    print(robot)

    #import matplotlib.pyplot as plt
    #plt.figure(0, figsize=(14,9))
    #plt.tick_params(labelsize=14)
    #ax = plt.gca()

    robot.plot_mobility_conv(ax=None, color='k', marker='o', label=r'$N_{{q_{{x,y,z}}}}$ = $N_{{k_{{x,y,z}}}}$')

    #fileslist = ['conv_fine/k27x27x27/q27x27x27/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k30x30x30/q30x30x30/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k108x108x108/q108x108x108/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k120x120x120/q120x120x120/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k132x132x132/q132x132x132/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k144x144x144/q144x144x144/Sio_DS1_TRANSPORT.nc',]

    #plot_mobility_conv(ax, fileslist, color='k', marker='o', label=r'$N_{{q_{{x,y,z}}}}$ = $N_{{k_{{x,y,z}}}}$')

    #fileslist = ['conv_fine/k27x27x27/q54x54x54/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k30x30x30/q60x60x60/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k66x66x66/q132x132x132/Sio_DS1_TRANSPORT.nc',
    #             'conv_fine/k72x72x72/q144x144x144/Sio_DS1_TRANSPORT.nc']

    #plot_mobility_conv(ax, fileslist, color='r', marker='x', label=r'$N_{{q_{{x,y,z}}}}$ = $2 N_{{k_{{x,y,z}}}}$')

    #plt.legend(loc='best',fontsize=14)
    #plt.show()
