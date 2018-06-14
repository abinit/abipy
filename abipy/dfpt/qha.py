# coding: utf-8
from __future__ import print_function, division, absolute_import

import numpy as np
import os
from scipy.interpolate import UnivariateSpline

from monty.collections import dict2namedtuple
from pymatgen.analysis.eos import EOS
from abipy.core.func1d import Function1D
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.phonons import PhdosFile, PhononBandsPlotter
import abipy.core.abinit_units as abu


eVA3_GPa = 160.21766208 # 1 eV/A^3 to GPa


class QHA(object):
    """
    Object to extract results in the quasi-harmonic approximation from electronic and phonon calculations
    at different volumes.
    Provides some basic methods and plotting utils, plus a converter to write input files for phonopy-qha or to
    generate an instance of phonopy.qha.QHA. These can be used to obtain other quantities and plots.
    Does not include electronic entropic contributions for metals.
    """

    def __init__(self, structures, doses, energies, eos_name='vinet', pressure=0):
        """
        Args:
            structures: list of structures at different volumes.
            doses: list of |PhononDos| at volumes corresponding to the structures.
            energies: list of SCF energies for the structures in eV.
            eos_name: string indicating the expression used to fit the energies. See pymatgen.analysis.eos.EOS.
            pressure: value of the pressure in GPa that will be considered in the p*V contribution to the energy.
        """

        self.structures = structures
        self.doses = doses
        self.energies = np.array(energies)
        self.eos = EOS(eos_name)
        self.pressure = pressure

        self.volumes = np.array([s.volume for s in structures])
        self.iv0 = np.argmin(energies)

    @property
    def nvols(self):
        return len(self.volumes)

    @property
    def natoms(self):
        return len(self.structures[0])

    def fit_energies(self, tstart=0, tstop=800, num=100):
        """
        Performs a fit of the energies as a function of the volume at different temperatures.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns:
            `namedtuple` with the following attributes::

                tot_en: numpy array with shape (nvols, num) with the energies used for the fit
                fits: list of subclasses of pymatgen.analysis.eos.EOSBase, depending on the type of
                    eos chosen. Contains the fit for the energies at the different temperatures.
                min_en: numpy array with the minimum energies for the list of temperatures
                min_vol: numpy array with the minimum volumes for the list of temperatures
                temp: numpy array with the temperatures considered

        """

        tmesh = np.linspace(tstart, tstop, num)

        # array with phonon energies and shape (n_vol, n_temp)
        ph_energies = np.array([dos.get_free_energy(tstart, tstop, num).values for dos in self.doses])

        tot_en = self.energies[np.newaxis, :].T + ph_energies + self.volumes[np.newaxis, :].T * self.pressure / eVA3_GPa

        # list of fits objects, one for each temperature
        fits = [self.eos.fit(self.volumes, e) for e in tot_en.T]

        # list of minimum volumes and energies, one for each temperature
        min_volumes = np.array([fit.v0 for fit in fits])
        min_energies = np.array([fit.e0 for fit in fits])

        return dict2namedtuple(tot_en=tot_en, fits=fits, min_en=min_energies, min_vol=min_volumes, temp=tmesh)

    def set_eos(self, eos_name):
        """
        Updates the EOS used for the fit.

        Args:
            eos_name: string indicating the expression used to fit the energies. See pymatgen.analysis.eos.EOS.
        """

        self.eos = EOS(eos_name)

    @classmethod
    def from_files(cls, gsr_files_paths, phdos_files_paths):
        """
        Creates an instance of QHA from a list og GSR files and a list PHDOS.nc files.
        The list should have the same size and the volumes should match.

        Args:
            gsr_files_paths: list of paths to GSR files
            phdos_files_paths: list of paths to PHDOS.nc files

        Returns:
            A new instance of QHA
        """

        energies = []
        structures = []
        for gp in gsr_files_paths:
            with GsrFile.from_file(gp) as g:
                energies.append(g.energy)
                structures.append(g.structure)

        doses = [PhdosFile(dp).phdos for dp in phdos_files_paths]

        return cls(structures, doses, energies)

    @add_fig_kwargs
    def plot_energies(self, tstart=0, tstop=800, num=10, ax=None, **kwargs):
        """
        Plots the energies as a function of volume at different temperatures.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 10.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """

        f = self.fit_energies(tstart, tstop, num)

        ax, fig, plt = get_ax_fig_plt(ax)
        xmin, xmax = np.floor(self.volumes.min() * 0.97), np.ceil(self.volumes.max() * 1.03)
        x = np.linspace(xmin, xmax, 100)

        for fit, e, t in zip(f.fits, f.tot_en.T - self.energies[self.iv0], f.temp):
            ax.scatter(self.volumes, e, label=t, color='b', marker='x', s=5)
            ax.plot(x, fit.func(x) - self.energies[self.iv0], color='b', lw=1)

        ax.plot(f.min_vol, f.min_en - self.energies[self.iv0] , color='r', lw=1, marker='x', ms=5)

        ax.set_xlabel('V (${\AA}^3$)')
        ax.set_ylabel('E (eV)')

        return fig

    def get_thermal_expansion_coeff(self, tstart=0, tstop=800, num=100):
        """
        Calculates the thermal expansion coefficient as a function of temperature, using
        finite difference on the fitted values of the volume as a function of temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns: |Function1D|
        """

        f = self.fit_energies(tstart, tstop, num)

        dt = f.temp[1] - f.temp[0]
        alpha = (f.min_vol[2:] - f.min_vol[:-2]) / (2 * dt) / f.min_vol[1:-1]

        return Function1D(f.temp[1:-1], alpha)

    @add_fig_kwargs
    def plot_thermal_expansion_coeff(self, tstart=0, tstop=800, num=100, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """

        ax, fig, plt = get_ax_fig_plt(ax)

        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 2

        if "color" not in kwargs:
            kwargs["color"] = "b"

        alpha = self.get_thermal_expansion_coeff(tstart, tstop, num)

        ax.plot(alpha.mesh, alpha.values, **kwargs)
        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')

        ax.set_xlim(tstart, tstop)

        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig

    @add_fig_kwargs
    def plot_vol_vs_t(self, tstart=0, tstop=800, num=100, ax=None, **kwargs):
        """
        Plots the volume as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        f = self.fit_energies(tstart, tstop, num)

        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 2

        if "color" not in kwargs:
            kwargs["color"] = "b"

        ax.plot(f.temp, f.min_vol, **kwargs)
        ax.set_xlabel('T (K)')
        ax.set_ylabel(r'V (${\AA}^3$)')

        ax.set_xlim(tstart, tstop)

        return fig

    @add_fig_kwargs
    def plot_phbs(self, phbands, temperatures=None, t_max=1000, colormap="plasma", **kwargs):
        """
        Given a list of |PhononBands| plots the band structures with a color depending on
        the temperature using a |PhononBandsPlotter|.
        If temperatures are not given they will be deduced inverting the dependence of volume
        with respect to the temperature. If a unique volume could not be identified an error will
        be raised.

        Args:
            phbands: list of |PhononBands| objects.
            temperatures: list of temperatures.
            t_max: maximum temperature considered for plotting.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """

        if temperatures is None:
            tv = self.get_t_for_vols([b.structure.volume for b in phbands], t_max=t_max)
            temperatures = []
            for b, t in zip(phbands, tv):
                if len(t) != 1:
                    raise ValueError("Couldn't find a single temperature for structure with "
                                     "volume {}. Found {}: {}".format(b.structure.volume, len(t), list(t)))
                temperatures.append(t[0])

        temperatures_str = ["{:.0f} K".format(t) for t in temperatures]

        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)
        colors = [cmap(t / max(temperatures)) for t in temperatures]
        labels_phbs = zip(temperatures_str, phbands)

        pbp = PhononBandsPlotter(labels_phbs)
        pbp._LINE_COLORS = colors
        pbp._LINE_STYLES = ['-']

        fig = pbp.combiplot(show=False, **kwargs)

        # ax, fig, plt = get_ax_fig_plt(ax)
        # cmap = plt.get_cmap(colormap)
        # for b, t in zip(phbands, temperatures):
        #     color = cmap(t / max(temperatures))
        #     b.plot(ax=ax, units=units, match_bands=True, show=False, color=color, label="")
        #     ax.lines[-1].set_label("{%d} K".format(t))
        #
        # ax.legend(loc="best")

        return fig

    def get_vol_at_t(self, t):
        """
        Calculates the volume corresponding to a specific temperature.

        Args:
            t: a temperature in K

        Returns:
            The volume
        """

        f = self.fit_energies(t, t, 1)

        return f.min_vol[0]

    def get_t_for_vols(self, vols, t_max=1000):
        """
        Find the temperatures corresponding to a specific volume.
        The search is performed interpolating the V(T) dependence with a spline and
        finding the roots with of V(t) - v.
        It may return more than one temperature for a volume in case of non monotonic behavior.

        Args:
            vols: list of volumes
            t_max: maximum temperature considered for the fit

        Returns:
            A list of lists of temperatures. For each volume more than one temperature can
            be identified.
        """

        if not isinstance(vols, (list, tuple, np.ndarray)):
            vols = [vols]


        f = self.fit_energies(0, t_max, t_max+1)

        temps = []
        for v in vols:
            spline = UnivariateSpline(f.temp, f.min_vol - v, s=0)
            temps.append(spline.roots())

        return temps

    def write_phonopy_qha_inputs(self, tstart=0, tstop=2100, num=211, path=None):
        """
        Writes nvols thermal_properties-i.yaml files that can be used as inputs for phonopy-qha.
        Notice that phonopy apparently requires the value of the 300 K temperature to be present
        in the list. Choose the values of tstart, tstop and num to satisfy this condition.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 211.
            path: a path to a folder where the files will be stored
        """

        if path is None:
            path = os.getcwd()

        np.savetxt(os.path.join(path, 'e-v.dat'), np.array([self.volumes, self.energies]).T, fmt='%.10f')

        # generator for thermal_properties.yaml extracted from phonopy:
        # phonopy.phonon.thermal_properties.ThermalProperties._get_tp_yaml_lines
        for j, d in enumerate(self.doses):
            lines = []
            lines.append("# Thermal properties / unit cell (natom)")
            lines.append("")
            lines.append("unit:")
            lines.append("  temperature:   K")
            lines.append("  free_energy:   kJ/mol")
            lines.append("  entropy:       J/K/mol")
            lines.append("  heat_capacity: J/K/mol")
            lines.append("")
            lines.append("natom: %5d" % (self.natoms))

            lines.append("num_modes: %d" % (3 * self.natoms))
            lines.append("num_integrated_modes: %d" % (3 * self.natoms))

            lines.append("")
            lines.append("zero_point_energy: %15.7f" % (d.zero_point_energy * abu.e_Cb * abu.Avogadro ))
            lines.append("high_T_entropy:    %15.7f" % 0) # high_T_entropy is not used in QHA
            lines.append("")
            lines.append("thermal_properties:")
            fe = d.get_free_energy(tstart, tstop, num) * abu.e_Cb * abu.Avogadro
            entropy = d.get_entropy(tstart, tstop, num) * abu.e_Cb * abu.Avogadro
            cv = d.get_cv(tstart, tstop, num) * abu.e_Cb * abu.Avogadro
            temperatures = cv.mesh
            for i, t in enumerate(temperatures):
                lines.append("- temperature:   %15.7f" % t)
                lines.append("  free_energy:   %15.7f" % (fe.values[i] / 1000))
                lines.append("  entropy:       %15.7f" % entropy.values[i])
                # Sometimes 'nan' of C_V is returned at low temperature.
                if np.isnan(cv.values[i]):
                    lines.append("  heat_capacity: %15.7f" % 0)
                else:
                    lines.append("  heat_capacity: %15.7f" % cv.values[i])
                lines.append("  energy:        %15.7f" %
                             (fe.values[i] / 1000 + entropy.values[i] * t / 1000))

                lines.append("")

            with open(os.path.join(path, "thermal_properties-{}.yaml".format(j)), 'wt') as f:
                f.write("\n".join(lines))

    def get_phonopy_qha(self, tstart=0, tstop=2100, num=211, eos='vinet', t_max=None, energy_plot_factor=None):
        """
        Creates an instance of phonopy.qha.QHA that can be used generate further plots and output data.
        The object is returned right after the construction. The "run()" method should be executed
        before getting results and plots.
        Notice that phonopy apparently requires the value of the 300 K temperature to be present
        in the list. Choose the values of tstart, tstop and num to satisfy this condition.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 211.
            eos: the expression used to fit the energies in phonopy. Possible values are "vinet",
                "murnaghan" and "birch_murnaghan". Passed to phonopy's QHA.
            t_max: maximum temperature. Passed to phonopy's QHA.
            energy_plot_factor: factor multiplying the energies. Passed to phonopy's QHA.

        Returns:
            An instance of phonopy.qha.QHA
        """

        try:
            from phonopy.qha import QHA as QHA_phonopy
        except ImportError as exc:
            print("Phonopy is required to generate the QHA phonopy object")
            raise exc

        fe = np.zeros((num, self.nvols))
        entropy = np.zeros((num, self.nvols))
        cv = np.zeros((num, self.nvols))
        for j, d in enumerate(self.doses):
            fe[:,j] = (d.get_free_energy(tstart, tstop, num) * abu.e_Cb * abu.Avogadro / 1000).values
            entropy[:,j] = (d.get_entropy(tstart, tstop, num) * abu.e_Cb * abu.Avogadro).values
            cv[:,j] = (d.get_cv(tstart, tstop, num) * abu.e_Cb * abu.Avogadro).values
        temperatures = self.doses[0].get_cv(tstart, tstop, num).mesh

        en = self.energies + self.volumes * self.pressure / eVA3_GPa

        qha_p = QHA_phonopy(self.volumes, en, temperatures, cv, entropy, fe, eos, t_max, energy_plot_factor)

        return qha_p