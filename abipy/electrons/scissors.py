# coding: utf-8
"""Scissors operator."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from six.moves import cPickle as pickle
from collections import OrderedDict
from monty.collections import AttrDict
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt


__all__ = [
    "Scissors",
    "ScissorsBuilder",
]


class ScissorsError(Exception):
    """Base class for the exceptions raised by :class:`Scissors`"""


class Scissors(object):
    """
    This object represents an energy-dependent scissors operator.
    The operator is defined by a list of domains (energy intervals)
    and a list of functions defined in these domains.
    The domains should fulfill the constraints documented in the main constructor.

    .. note::

        eV units are assumed.

    The standard way to create this object is via the methods provided by the factory class :class:`ScissorBuilder`.
    Once the instance has been created, one can correct the band structure by calling the `apply` method.
    """
    Error = ScissorsError

    def __init__(self, func_list, domains, residues, bounds=None):
        """
        Args:
            func_list: List of callable objects. Each function takes an eigenvalue and returns
                the corrected value.
            domains: Domains of each function. List of tuples [(emin1, emax1), (emin2, emax2), ...]
            bounds: Specify how to handle energies that do not fall inside one of the domains.
                At present, only constant boundaries are implemented.
            residues: A list of the residues of the fitting per domain

        .. note::

            #. Domains should not overlap, cover e0mesh, and given in increasing order.

            #. Holes are permitted but the interpolation will raise an exception if the
               eigenvalue falls inside the hole.

            #. Errors contains a list of the fitting errors per domain

        """
        # TODO Add consistency check.
        self.func_list = func_list
        self.domains = np.atleast_2d(domains)
        self.residues = residues
        assert len(self.func_list) == len(self.domains)

        # Treat the out-of-boundary conditions. func_low and func_high are used to handle energies
        # that are below or above the min/max energy given in domains.
        blow, bhigh = "c", "c"
        if bounds is not None:
            blow, bhigh = bounds[0][0], bounds[0][1]

        if blow.lower() == "c":
            try:
                self.func_low = lambda x: float(bounds[0][1])
            except:
                x_low = self.domains[0,0]
                fx_low = func_list[0](x_low)
                self.func_low = lambda x: fx_low
        else:
            raise NotImplementedError("Only constant boundaries are implemented")

        if bhigh.lower() == "c":
            try:
                self.func_high = lambda x: float(bounds[1][1])
            except:
                x_high = self.domains[1, -1]
                fx_high = func_list[-1](x_high)
                self.func_high = lambda x: fx_high
        else:
            raise NotImplementedError("Only constant boundaries are implemented")

        # This counter stores the number of points that are out of bounds.
        self.out_bounds = np.zeros(3, np.int)

    def apply(self, eig):
        """Correct the eigenvalue eig (eV units)."""
        # Get the list of domains.
        domains = self.domains

        if eig < domains[0,0]:
            # Eig is below the first point of the first domain.
            # Call func_low
            print("left ", eig, " < ", domains[0,0])
            self.out_bounds[0] += 1
            return self.func_low(eig)

        if eig > domains[-1,1]:
            # Eig is above the last point of the last domain.
            # Call func_high
            print("right ", eig, " > ", domains[-1,1])
            self.out_bounds[1] += 1
            return self.func_high(eig)

        # eig is inside the domains: find the domain
        # and call the corresponding function.
        for idx, dms in enumerate(domains):
            if dms[1] >= eig >= dms[0]:
                return self.func_list[idx](eig)

        self.out_bounds[2] += 1
        raise self.Error("Cannot find location of eigenvalue %s in domains:\n%s" % (eig, domains))


class ScissorsBuilder(object):
    """
    This object facilitates the creation of :class:`Scissors` instances.

    Usage:

        builder = ScissorsBuilder.from_file("out_SIGRES.nc")

        # To plot the QP results as function of the KS energy:
        builder.plot_qpe_vs_e0()

        # To select the domains esplicitly (optional but highly recommended)
        builder.build(domains_spin=[[-10, 6.02], [6.1, 20]])

        # To compare the fitted results with the ab-initio data:
        builder.plot_fit()

        # To plot the corrected bands:
        builder.plot_qpbands(abidata.ref_file("si_nscf_WFK.nc"))
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Generate object from (SIGRES.nc) file. Main entry point for client code.
        """
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            return cls(qps_spin=ncfile.qplist_spin, sigres_ebands=ncfile.ebands)

    @classmethod
    def pickle_load(cls, filepath):
        """Load the object from a pickle file."""
        with open(filepath, "rb") as fh:
            d = AttrDict(pickle.load(fh))
            # Costruct the object and compute the scissors.
            new = cls(d.qps_spin, d.sigres_ebands)
            new.build(d.domains_spin, d.bounds_spin)
            return new

    def pickle_dump(self, filepath, protocol=-1):
        """Save the object in Pickle format"""
        assert all(s1 == s2 for s1, s2 in zip(self.domains_spin.keys(), self.bounds_spin.keys()))
        assert all(s1 == s2 for s1, s2 in zip(self.domains_spin.keys(), range(self.nsppol)))

        bounds_spin = None
        if any(v is not None for v in self.bounds_spin.values()):
            bounds_spin = [a.tolist() for a in self.bounds_spin.values()]

        # This trick is needed because we cannot pickle bound methods of the scissors operator.
        d = dict(qps_spin=self._qps_spin,
                 sigres_ebands=self.sigres_ebands,
                 domains_spin=[a for a in self.domains_spin.values()],
                 bounds_spin=bounds_spin)

        with open(filepath, "wb") as fh:
            pickle.dump(d, fh, protocol=protocol)

    def __init__(self, qps_spin, sigres_ebands):
        """
        Args:
            qps_spin: List of :class:`QPlist`, for each spin.
            sigres_ebands: |ElectronBands| obtained from the SIGRES file
        """
        # Sort quasiparticle data by e0.
        self._qps_spin = tuple([qps.sort_by_e0() for qps in qps_spin])

        # Compute the boundaries of the E0 mesh.
        e0min, e0max = np.inf, -np.inf
        for qps in self._qps_spin:
            e0mesh = qps.get_e0mesh()
            e0min = min(e0min, e0mesh[0])
            e0max = max(e0max, e0mesh[-1])

        self._e0min, self._e0max = e0min, e0max

        # The KS bands stored in the sigres file (used to compute automatically the boundaries)
        self.sigres_ebands = sigres_ebands

        # Start with default values for domains.
        self.build()

    @property
    def nsppol(self):
        """Number of spins."""
        return len(self._qps_spin)

    @property
    def e0min(self):
        """Minimum KS energy in eV (takes into account spin)"""
        return self._e0min

    @property
    def e0max(self):
        """Maximum KS energy in eV (takes into account spin)"""
        return self._e0max

    @property
    def scissors_spin(self):
        """Returns a tuple of :class:`Scissors` indexed by the spin value."""
        try:
            return self._scissors_spin
        except AttributeError:
            raise AttributeError("Call self.build to create the scissors operator")

    def build(self, domains_spin=None, bounds_spin=None, k=3):
        """
        Build the scissors operator.

        Args:
            domains_spin: list of domains in eV for each spin. If domains is None,
                domains are computed automatically from the sigres bands
                (two domains separated by the middle of the gap).
            bounds_spin: Options specifying the boundary conditions (not used at present)
            k: Parameter defining the order of the fit.
        """
        nsppol = self.nsppol

        # The parameters defining the scissors operator
        self.domains_spin = OrderedDict()
        self.bounds_spin = OrderedDict()

        if domains_spin is None:
            # Use sigres_ebands and the position of the homo, lumo to compute the domains.
            domains_spin = nsppol * [None]
            e_bands = self.sigres_ebands
            for spin in e_bands.spins:
                gap_mid = (e_bands.homos[spin].eig + e_bands.lumos[spin].eig) / 2
                domains_spin[spin] = [[self.e0min - 0.2 * abs(self.e0min), gap_mid],
                                      [gap_mid, self.e0max + 0.2 * abs(self.e0max)]]
                #print("domains", domains_spin[spin])
        else:
            if nsppol == 1:
                domains_spin = np.reshape(domains_spin, (1, -1, 2))
            elif nsppol == 2:
                assert len(domains_spin) == nsppol
                if bounds_spin is not None: assert len(bounds_spin) == nsppol
            else:
                raise ValueError("Wrong number of spins %d" % nsppol)
            #if len(domains_spin) != nsppol:
            #    raise ValueError("len(domains_spin) == %s != nsppol %s" % (len(domains_spin), nsppol))

        # Construct the scissors operator for each spin.
        scissors_spin = nsppol * [None]
        for spin, qps in enumerate(self._qps_spin):
            bounds = None if not bounds_spin else bounds_spin[spin]
            scissors_spin[spin] = qps.build_scissors(domains_spin[spin], bounds=bounds, k=k, plot=False)

            # Save input so that we can reconstruct Scissors.
            self.domains_spin[spin] = domains_spin[spin]
            self.bounds_spin[spin] = bounds

        self._scissors_spin = scissors_spin
        return domains_spin

    @add_fig_kwargs
    def plot_qpe_vs_e0(self, with_fields="all", **kwargs):
        """Plot the quasiparticle corrections as function of the KS energy."""
        ax_list = None
        for spin, qps in enumerate(self._qps_spin):
            kwargs["title"] = "spin %s" % spin
            fig = qps.plot_qps_vs_e0(with_fields=with_fields, ax_list=ax_list, show=False, **kwargs)
            ax_list = fig.axes

        return fig

    @add_fig_kwargs
    def plot_fit(self, ax=None, fontsize=8, **kwargs):
        """
        Compare fit functions with input quasi-particle corrections.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for titles and legend.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for spin in range(self.nsppol):
            qps = self._qps_spin[spin]
            e0mesh, qpcorrs = qps.get_e0mesh(), qps.get_qpeme0().real

            ax.scatter(e0mesh, qpcorrs, label="Input QP corrections, spin %s" % spin)
            scissors = self._scissors_spin[spin]
            intp_qpc = [scissors.apply(e0) for e0 in e0mesh]
            ax.plot(e0mesh, intp_qpc, label="Scissors operator, spin %s" % spin)

        ax.grid(True)
        ax.set_xlabel('KS energy (eV)')
        ax.set_ylabel('QP-KS (eV)')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def plot_qpbands(self, bands_filepath, bands_label=None, dos_filepath=None, dos_args=None, **kwargs):
        """
        Correct the energies found in the netcdf file bands_filepath and plot the band energies (both the initial
        and the corrected ones) with matplotlib. The plot contains the KS and the QP DOS if dos_filepath is not None.

        Args:
            bands_filepath: Path to the netcdf file containing the initial KS energies to be corrected.
            bands_label String used to label the KS bands in the plot.
            dos_filepath: Optional path to a netcdf file with the initial KS energies on a homogeneous k-mesh
                (used to compute the KS and the QP dos)
            dos_args: Dictionary with the arguments passed to get_dos to compute the DOS
                Used if dos_filepath is not None.

            kwargs: Options passed to the plotter.

        Return: |matplotlib-Figure|
        """
        from abipy.abilab import abiopen, ElectronBandsPlotter

        # Read the KS band energies from bands_filepath and apply the scissors operator.
        with abiopen(bands_filepath) as ncfile:
            ks_bands = ncfile.ebands
            #structure = ncfile.structure

        qp_bands = ks_bands.apply_scissors(self._scissors_spin)

        # Read the band energies computed on the Monkhorst-Pack (MP) mesh and compute the DOS.
        ks_dos, qp_dos = None, None
        if dos_filepath is not None:
            with abiopen(dos_filepath) as ncfile:
                ks_mpbands = ncfile.ebands

            dos_args = {} if not dos_args else dos_args
            ks_dos = ks_mpbands.get_edos(**dos_args)
            # Compute the DOS with the modified QPState energies.
            qp_mpbands = ks_mpbands.apply_scissors(self._scissors_spin)
            qp_dos = qp_mpbands.get_edos(**dos_args)

        # Plot the LDA and the QPState band structure with matplotlib.
        plotter = ElectronBandsPlotter()

        bands_label = bands_label if bands_label is not None else os.path.basename(bands_filepath)
        plotter.add_ebands(bands_label, ks_bands, dos=ks_dos)
        plotter.add_ebands(bands_label + " + scissors", qp_bands, dos=qp_dos)

        #qp_marker: if int > 0, markers for the ab-initio QP energies are displayed. e.g qp_marker=50
        #qp_marker = 50
        #if qp_marker is not None:
        #    # Compute correspondence between the k-points in qp_list and the k-path in qp_bands.
        #    # TODO
        #    # WARNING: strictly speaking one should check if qp_kpoint is in the star of k-point.
        #    # but compute_star is too slow if written in pure python.
        #    x, y, s = [], [], []
        #    for ik_path, kpoint in enumerate(qp_bands.kpoints):
        #        #kstar = kpoint.compute_star(structure.fm_symmops)
        #        for spin in range(self.nsppol):
        #            for ik_qp, qp in enumerate(self._qps_spin[spin]):
        #                #if qp.kpoint in kstar:
        #                if qp.kpoint == kpoint:
        #                    x.append(ik_path)
        #                    y.append(np.real(qp.qpe))
        #                    s.append(qp_marker)
        #    plotter.set_marker("ab-initio QP", [x, y, s])

        return plotter.combiplot(**kwargs)
