# coding: utf-8
"""Objects to run and analyze the frozen phonons generated from displacements of atoms"""
from __future__ import print_function, division, absolute_import

import numpy as np
import scipy.optimize as optimize

from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from abipy.core.abinit_units import phfactor_ev2units, amu_emass, Bohr_Ang, eV_Ha
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt


def quadratic_fit_function(xx, aa, bb):
    """
    Simple quadratic function used as a default for the fitting.

    Args:
        xx: the variable
        aa: the coefficient of the quadratic term
        bb: the constant term
    """

    return aa * xx ** 2 + bb

class FrozenPhonon(object):
    """
    Class defining a set of structures with displaced atoms.
    Provides methods to generate, interpolate and plot the data.
    """

    def __init__(self, original_structure, original_displ_cart, structures, normalized_displ_cart, etas,
                 qpt_frac_coords, scale_matrix, energies=None):
        """

        Args:
            original_structure: the original |Structure| object without any displacement.
            original_displ_cart: the original displacement used for generating the structures in cartesian
                coordinates. Should be linked to original_structure and should contain the phdispl_cart obtained
                from the phonon calculation.
            structures: list of structures with the displaced atoms
            normalized_displ_cart: displacement applied to the supercell to generate the list of structures.
                The normalization defines the largest displacement of one atom at 1 Angstrom.
            etas: list of (signed) factors that will multiply the normalized displacement to obtain the list of
                structure.
            qpt_frac_coords: fractional coordinates of the qpoint corresponding to the displacement.
            scale_matrix: the scaling matrix between the original cell and the supercell.
            energies: energies in eV corresponding to the structures.
        """

        self.original_structure = original_structure
        self.original_displ_cart = original_displ_cart
        self.structures = structures
        self.normalized_displ_cart = normalized_displ_cart
        self.etas = np.array(etas)
        self.qpt_frac_coords = qpt_frac_coords
        self.scale_matrix = scale_matrix
        self._energies = energies

    @property
    def energies(self):
        """
        Energies in eV corresponding to the the structures.
        """
        return self._energies

    @energies.setter
    def energies(self, energies):
        """
        Set the number of energies corresponding to the structures, should have the same length as structures
        """
        if len(energies) != self.n_displ:
            raise ValueError("The size of the enegies list does not match the number of internal structures.")
        self._energies = energies

    @property
    def qpt_cart_coords(self):
        """
        The cartesian coordinates of the qpoint.
        """
        return self.original_structure.lattice.reciprocal_lattice.get_cartesian_coords(self.qpt_frac_coords)

    @property
    def n_displ(self):
        """
        Number of displacements.
        """
        return len(self.structures)

    @lazy_property
    def ieta0(self):
        """
        The index corresponding to the structure with no displacements.
        """
        try:
            return self.etas.tolist().index(0)
        except ValueError:
            raise ValueError("The structure with no displacement is not present in the list.")

    @classmethod
    def from_phbands(cls, phbands, qpt_frac_coords, imode, etas, scale_matrix=None, max_supercell=None):
        """
        Create an instace of FrozenPhonon using the eigendisplacements from a |PhononBands|

        Args:
            phbands: a |PhononBands| instance.
            qpt_frac_coords: q vector in reduced coordinate in reciprocal space or index of the qpoint.
            imode: index of the mode.
            etas: list of amplitudes of the displacement to be applied to the system. Will correspond to the
                largest displacement of one atom in Angstrom.
            scale_matrix: the scaling matrix of the supercell. If None a scaling matrix suitable for
                the qpoint will be determined.
            max_supercell: mandatory if scale_matrix is None, ignored otherwise. Defines the largest
                supercell in the search for a scaling matrix suitable for the q point.

        Returns:
            A FrozenPhonon.
        """

        qind = phbands.qindex(qpt_frac_coords)
        original_displ_cart = phbands.phdispl_cart[qind, imode].reshape((-1, 3))

        # first extract data with the normalized displacement (eta=1)
        normalized_fp = phbands.get_frozen_phonons(qind, imode, 1, scale_matrix, max_supercell)

        structures = []

        for n in etas:
            structures.append(phbands.get_frozen_phonons(qind, imode, n, normalized_fp.scale_matrix).structure)

        return cls(phbands.structure, original_displ_cart, structures, normalized_fp.displ, etas,
                   phbands.qpoints[qind].frac_coords, normalized_fp.scale_matrix)

    @lazy_property
    def mass_factor(self):
        """
        The factor accounting for the different masses and displacement of each atom
        """
        # the norms of each atomic displacement of the normalized displacement
        displ_norms = np.linalg.norm(self.normalized_displ_cart, axis=1)

        mass_factor = sum(site.specie.atomic_mass * (df**2)
                          for site, df in zip(self.structures[0], displ_norms / displ_norms.max())) / 2.
        mass_factor *= amu_emass
        return mass_factor

    def _freq_to_quad_coeff(self, freq):
        """
        Helper function to convert from a frequency in eV to the quadratic coefficient corresponding to the
        energies in eV and the etas in Angstrom. Uses the conversions to a.u. to get the correct value.
        """

        return freq ** 2 * self.mass_factor * eV_Ha / Bohr_Ang ** 2

    def _quad_coeff_to_freq(self, coeff):
        """
        Helper function to convert from the quadratic coefficient corresponding to the fit of energies in eV and
        the etas in Angstrom to a frequency in eV. Uses the conversions to a.u. to get the correct value.
        """

        return np.sqrt(coeff / self.mass_factor / eV_Ha) * Bohr_Ang

    def fit_to_frequency(self, fit_function=None, units="eV", min_fit_eta=None, max_fit_eta=None):
        """
        Uses the energies and the displacements to calculate the phonon frequency corresponding to the quadratic
        term of the fit.
        The fit is performed with scipy.optimize.curve_fit based on the function given in input and can also be
        limited number to a subset of the values of the displacements.


        Args:
            fit_function: a function that will be used to fit the data. The first parameter should be the coefficient
                of the quadratic term. If None a simple quadratic fit will be used.
            units: units of the output frequency. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            min_fit_eta: if not None represents minimum value allowed for the (signed) eta to be used in the fit.
            max_fit_eta: if not None represents maximum value allowed for the (signed) eta to be used in the fit.

        Returns:
            A namedtuple with 'freq': the values of the frequency extracted from the fit,
            'fit_params': the parameters obtained from the fit, 'cov': the estimated covariance of fit_params
            (see scipy.optimize.curve_fit documentation for more details).
        """

        if self.energies is None:
            raise ValueError("The energies are required to calculate the fit")

        if min_fit_eta is None:
            min_fit_eta = self.etas.min()
        if max_fit_eta is None:
            max_fit_eta = self.etas.max()

        indices = np.where((min_fit_eta <= self.etas) & (self.etas <= max_fit_eta))

        if fit_function is None:
            fit_function = quadratic_fit_function

        etas = self.etas[indices]
        energies = np.array(self.energies)[indices]

        params, cov = optimize.curve_fit(fit_function, etas, energies)

        # frequency in eV.
        freq = self._quad_coeff_to_freq(params[0])

        return dict2namedtuple(freq=freq*phfactor_ev2units(units), fit_params=params, cov=cov)

    @add_fig_kwargs
    def plot_fit_energies(self, fit_function=None, min_fit_eta=None, max_fit_eta=None, freq=None,
                          ax=None, **kwargs):
        """
        Fits the displacements etas to the energies. See fit_to_frequency() for more details.

        Args:
            fit_function: a function that will be used to fit the data. The first parameter should be the coefficient
                of the quadratic term. If None a simple quadratic fit will be used.
            min_fit_eta: if not None represents minimum value allowed for the (signed) eta to be used in the fit.
            max_fit_eta: if not None represents maximum value allowed for the (signed) eta to be used in the fit.
            freq: if not None the quadratic function with this frequency as a coefficient will be added to the plot.
                freq in eV. Requires the 0 displacement to be present in the list of etas.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns:
            |matplotlib-Figure|
        """

        if fit_function is None:
            fit_function = quadratic_fit_function

        fit_data = self.fit_to_frequency(fit_function, min_fit_eta=min_fit_eta, max_fit_eta=max_fit_eta)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        if "color" not in kwargs and "c" not in kwargs:
            kwargs["color"] = "blue"

        kwargs_points = dict(kwargs)
        kwargs_points["linestyle"] = "None"
        if "markersize" not in kwargs_points and "ms" not in kwargs_points:
            kwargs_points["markersize"] = 10
        if "marker" not in kwargs_points and "m" not in kwargs_points:
            kwargs_points["marker"] = "o"

        ax.plot(self.etas, self.energies, **kwargs_points)

        if "linewidth" not in kwargs and "lw" not in kwargs:
            kwargs["linewidth"] = 2
        kwargs["linestyle"] = "-"

        etas_plot = np.linspace(self.etas.min(), self.etas.max(), 100)
        ax.plot(etas_plot, fit_function(etas_plot, *fit_data.fit_params), **kwargs)

        if freq is not None:
            kwargs["linestyle"] = "--"
            e0 = self.energies[self.ieta0]
            en_freq = quadratic_fit_function(etas_plot, self._freq_to_quad_coeff(freq), e0)
            ax.plot(etas_plot, en_freq, **kwargs)

        ax.set_xlabel("Displacements (Ang)")
        ax.set_ylabel("Energy (eV)")

        return fig

    @add_fig_kwargs
    def plot_anharmonic_contribution(self, freq, relative=False, ax=None, **kwargs):
        """
        Plots the the absolute relative difference between the energies extracted from the frequency as
        quadratic coefficient and the calculated energies, giving an estimate of the anharmonic contribution
        Requires the 0 displacement to be present in the list of etas.

        Args:
            freq: phonon frequncy in eV
            relative: if True the plot will represent the relative difference with respect to the expected value
                obtained from the frequency, rather than the absolute difference.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns:
            |matplotlib-Figure|
        """

        if self.energies is None:
            raise ValueError("The energies are required to calculate the fit")

        self._freq_to_quad_coeff(freq)

        if "color" not in kwargs and "c" not in kwargs:
            kwargs["color"] = "blue"
        if "linewidth" not in kwargs and "lw" not in kwargs:
            kwargs["linewidth"] = 2

        e0 = self.energies[self.ieta0]
        en_freq = quadratic_fit_function(self.etas,self._freq_to_quad_coeff(freq), e0)

        diff = np.abs(en_freq - np.array(self.energies))

        if relative:
            diff = [d/(e-e0) * 100 if e != 0 else 0 for d, e in zip(diff, en_freq)]

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ax.plot(self.etas, diff, **kwargs)

        ax.set_xlabel("Displacements (Ang)")
        if relative:
            ax.set_ylabel("Relative energy difference (%)")
        else:
            ax.set_ylabel("Energy difference (eV)")

        return fig
