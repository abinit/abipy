# coding: utf-8
"""Tools to analyze the output files produced by Lobster."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import glob
import numpy as np

from collections import defaultdict, OrderedDict
from monty.string import marquee
from monty.collections import tree
from monty.io import zopen
from monty.termcolor import cprint
from monty.functools import lazy_property
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.abinit.pseudos import Pseudo
from abipy.core.func1d import Function1D
from abipy.core.mixins import BaseFile, NotebookWriter
from abipy.electrons.gsr import GsrFile
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_visible, set_ax_xylabels
from abipy.tools import duck


class _LobsterFile(BaseFile, NotebookWriter):
    """
    Base class for output files produced by lobster.
    """
    # These class attributes can be redefined in order to customize the plots.

    # Mapping L --> color used in plots.
    color_l = {"s": "red", "p": "blue", "d": "green", "f": "orange"}

    # \U starts an eight-character Unicode escape. raw strings do not work in python2.7
    # and we need a latex symbol to avoid errors in matplotlib --> replace myuparrow --> uparrow

    # Mapping spin --> title used in subplots that depend on (collinear) spin.
    spin2tex = {k: v.replace("myuparrow", "uparrow") for k, v in
            {0: r"$\sigma=\myuparrow$", 1: r"$\sigma=\downarrow$"}.items()}

    # TODOL ls
    params_orbs = {
      "s":  {"style": dict(color="red", ls="-", lw=1), "latex": "s"},
      #
      "p_x": {"style": dict(color="blue", ls="-", lw=1), "latex": "p_x"},
      "p_y": {"style": dict(color="blue", ls="-", lw=1), "latex": "p_y"},
      "p_z": {"style": dict(color="blue", ls="-", lw=1), "latex": "p_z"},
      #
      "d_xy":  {"style": dict(color="green", ls="-", lw=1), "latex": "d_{xy}"},
      "d_yz":  {"style": dict(color="green", ls="-", lw=1), "latex": "d_{yz}"},
      "d_z^2": {"style": dict(color="green", ls="-", lw=1), "latex": "d_{z^2}"},
      "d_xz":  {"style": dict(color="green", ls="-", lw=1), "latex": "d_{xz}"},
      "d_x^2-y^2": {"style": dict(color="green", ls="-", lw=1), "latex": "d_{x^2-y^2}"},
      #
      "4f_y(3x^2-y^2)": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{y(3x^2-y^2)}"},
      "4f_xyz": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{xyz}"},
      "4f_yz^2": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{yz^2}"},
      "4f_z^3": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{z^3}"},
      "4f_xz^2": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{xz^2}"},
      "4f_z(x^2-y^2)": {"style": dict(color="orange", ls="-", lw=1), "latex": "4f_{z(x^2-y^2)}"},
    }

    def __str__(self):
        return self.to_string()

    def close(self):
        """Needed by ABC."""

    #@add_fig_kwargs
    #def plot_with_ebands(self, ebands, fontsize=12, **kwargs):
    #    """
    #    Plot bands + (COHP|COOP|DOSCAR) depending on the content of the file.

    #    Args:
    #        ebands: Path to ncfile with ebands or |ElectronBands| object.
    #        fontsize: fontsize for legends and titles

    #    Returns: |matplotlib-Figure|
    #    """
    #    ebands = ElectronBands.as_ebands(ebands)

    #    import matplotlib.pyplot as plt
    #    from matplotlib.gridspec import GridSpec
    #    fig = plt.figure()

    #    # Build grid.
    #    nrows, ncols = 1, 2
    #    gspec = GridSpec(nrows=nrows, ncols=ncols, width_ratios=(2, 1), wspace=0.05)

    #    # Bands and DOS will share the y-axis
    #    axmat = np.array((nrows, ncols), dtype=object)
    #    for icol in range(ncols):
    #        for irow in range(nrows):
    #            axmat[irow, icol] = plt.subplot(gspec[irow, icol], sharey=None if irow == 0 else axmat[irow, icol-1])

    #    #axmat, fig, plt = get_axarray_fig_plt(None, nrows=self.nsppol, ncols=len(entries) + 1,
    #    #                                        sharex=False, sharey=True, squeeze=False)

    #    for ix, ax in enumerate(axmat.ravel()):
    #        if ix == 0:
    #            ebands.plot(e0="fermie", ax=ax, show=False)
    #        else:
    #            self.plot(ax=ax, exchange_xy=True, show=False)

    #        #if ix != 0:
    #        #    set_visible(ax, False, "ylabel")
    #        #if self.nsppol == 2 and spin == 0:
    #        #    set_visible(ax, False, "xlabel")

    #    return fig


class CoxpFile(_LobsterFile):
    """
    Wrapper class for the crystal orbital projections produced by Lobster.
    Wraps both a COOP and a COHP.
    Can contain both the total and orbitalwise projections.

    .. attribute:: cop_type

        String. Either "coop" or "cohp".

    .. attribute:: type_of_index

        Dictionary mappping site index to element string.

    .. attribute:: energies

        List of energies. Shifted such that the Fermi level lies at 0 eV.

    .. attribute:: total

        A dictionary with the values of the total overlap, not projected over orbitals.
        The dictionary should have the following nested keys: a tuple with the index of the sites
        considered as a pair (0-based, e.g. (0, 1)), the spin (i.e. 0 or 1), a "single"
        or "integrated" string indicating the value corresponding to the value of
        the energy or to the integrated value up to that energy.

    .. attribute:: partial

        A dictionary with the values of the partial crystal orbital projections.
        The dictionary should have the following nested keys: a tuple with the index of the sites
        considered as a pair (0-based, e.g. (0, 1)), a tuple with the string representing the
        projected orbitals for that pair (e.g. ("4s", "4p_x")), the spin (i.e. 0 or 1),
        a "single" or "integrated" string indicating the value corresponding to the value of
        the energy or to the integrated value up to that energy. Each dictionary should contain a
        numpy array with a list of COP values with the same size as energies.

    .. attribute:: averaged

        A dictionary with the values of the partial crystal orbital projections
        averaged over all atom pairs specified. The main key should indicate the spin (0 or 1)
        and the nested dictionary should have "single" and "integrated" as keys.

    .. attribute:: fermie

        value of the fermi energy in eV.
    """

    @property
    def site_pairs_total(self):
        """
        List of site pairs available for the total COP
        """
        return list(self.total.keys())

    @property
    def site_pairs_partial(self):
        """
        List of site pairs available for the partial COP
        """
        return list(self.partial.keys())

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance of CoxpFile from the files produce by Lobster.
        Accepts gzipped files.

        Args:
            filepath: path to the COHPCAR.lobster or COOPCAR.lobster.

        Returns:
            A CoxpFile.
        """
        # From lobster documentation
        # COHPCAR.lobster:
        # File that contains the pCOHPs as requested in the lobsterin file.
        # It resembles the format of TB-LMTO-ASA's COPL file, which is organized as follows:
        # - Starting in line 3, the labels for the interactions are presented, followed by the
        # actual data.
        # - Column 1: energy axis, shifted such that the Fermi level lies at zero eV.
        # - Column 2: pCOHP averaged over all atom pairs specified
        # - Column 3: integrated pCOHP (IpCOHP) averaged over all atom pairs
        # - Column 4: pCOHP of the first interaction
        # - Column 5: IpCOHP of the first interaction
        # - and so on...
        # Note that in a spin-polarized calculation, the first set of the columns (2, 3, ..., 2N+3)
        # belongs to the first (up) spin and the other set (2N+4, 2N+5, ..., 4N+5) belongs to the
        # second (down) spin. Here N is the number of interactions.

        float_patt = r'-?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
        header_patt = re.compile(r'\s+(\d+)\s+(\d+)\s+(\d+)\s+('+float_patt+
                                 r')\s+('+float_patt+r')\s+('+float_patt+r')')
        pair_patt = re.compile(r'No\.\d+:([a-zA-Z]+)(\d+)(?:\[([a-z0-9_\-^]+)\])?->([a-zA-Z]+)(\d+)(?:\[([a-z0-9_\-^]+)\])?')

        new = cls(filepath)

        with zopen(filepath, "rt") as f:
            # Find the header
            for line in f:
                match = header_patt.match(line.rstrip())
                if match:
                    n_column_groups = int(match.group(1))
                    n_spin = int(match.group(2))
                    n_en_steps = int(match.group(3))
                    min_en = float(match.group(4))
                    max_en = float(match.group(5))
                    new.fermie = float(match.group(6))
                    break
            else:
                raise ValueError("Can't find the header in file {}".format(filepath))

            n_pairs = n_column_groups-1

            count_pairs = 0
            pairs_data = []
            new.type_of_index = {}
            # Parse the pairs considered
            for line in f:
                match = pair_patt.match(line.rstrip())
                if match:
                    # Adds a tuple: [type1, index1, orbital1, type2, index2, orbital2]
                    # with orbital1, orbital2 = None if the pair is not orbitalwise
                    type1, index1, orbital1, type2, index2, orbital2 = match.groups()
                    # 0-based indexing
                    index1 = int(index1) - 1
                    index2 = int(index2) - 1
                    pairs_data.append([type1, index1, orbital1, type2, index2, orbital2])
                    if index1 in new.type_of_index: assert new.type_of_index[index1] == type1
                    new.type_of_index[index1] = type1
                    if index2 in new.type_of_index: assert new.type_of_index[index2] == type2
                    new.type_of_index[index2] = type2
                    count_pairs += 1
                    if count_pairs == n_pairs:
                        break

            spins = [0, 1][:n_spin]

            data = np.fromstring(f.read(), dtype=np.float, sep=' ').reshape([n_en_steps, 1+n_spin*n_column_groups*2])

            # Initialize and fill results
            new.energies = data[:, 0]
            new.averaged = defaultdict(dict)
            new.total = tree()
            new.partial = tree()

            for i, s in enumerate(spins):
                base_index = 1+i*n_column_groups*2
                new.averaged[s]['single'] = data[:, base_index].copy()
                new.averaged[s]['integrated'] = data[:, base_index+1].copy()
                # NB (i, j) --> (j, i) symmetry is enforced to make API easier.
                for j, p in enumerate(pairs_data):
                    index1 = p[1]
                    index2 = p[4]
                    if p[2] is not None:
                        # Partial
                        single = data[:, base_index+2*(j+1)].copy()
                        integrated = data[:, base_index+2*(j+1)+1].copy()
                        new.partial[(index1, index2)][(p[2], p[5])][s]['single'] = single
                        new.partial[(index2, index1)][(p[5], p[2])][s]['single'] = single
                        new.partial[(index1, index2)][(p[2], p[5])][s]['integrated'] = integrated
                        new.partial[(index2, index1)][(p[5], p[2])][s]['integrated'] = integrated
                    else:
                        # Total
                        single = data[:, base_index+2*(j+1)].copy()
                        integrated = data[:, base_index+2*(j+1)+1].copy()
                        new.total[(index1, index2)][s]['single'] = single
                        new.total[(index2, index1)][s]['single'] = single
                        new.total[(index1, index2)][s]['integrated'] = integrated
                        new.total[(index2, index1)][s]['integrated'] = integrated

        new.cop_type = "unknown"
        if "COOPCAR.lobster" in filepath: new.cop_type = "coop"
        if "COHPCAR.lobster" in filepath: new.cop_type = "cohp"
        new.nsppol = len(new.averaged)

        return new

    @lazy_property
    def functions_pair_lorbitals(self):
        """
        Extracts a dictionary with keys pair, orbital, spin and containing a |Function1D| object resolved
        for l orbitals.
        """
        if not self.partial:
            raise RuntimeError("Partial orbitals not calculated.")

        results = tree()
        for pair, pair_data in self.partial.items():
            # Check if the symmetric has already been calculated
            if (pair[1], pair[0]) in results:
                for orbs, orbs_data in results[(pair[1], pair[0])].items():
                    results[pair][(orbs[1], orbs[0])] = orbs_data
                continue

            # For each look at all orbital possibilities
            for orbs, orbs_data in pair_data.items():
                k = (orbs[0].split("_")[0], orbs[1].split("_")[0])
                if k in results[pair]:
                    for spin in orbs_data.keys():
                        results[pair][k][spin] = results[pair][k][spin] + Function1D(self.energies,orbs_data[spin]['single'])
                else:
                    for spin in orbs_data.keys():
                        results[pair][k][spin] = Function1D(self.energies, orbs_data[spin]['single'])

        return results

    @lazy_property
    def functions_pair_morbitals(self):
        """
        Extracts a dictionary with keys pair, orbital, spin and containing a |Function1D| object resolved
        for l and m orbitals.
        """
        if not self.partial:
            raise RuntimeError("Partial orbitals not calculated.")

        results = tree()
        for pair, pair_data in self.partial.items():
            for orbs, orbs_data in pair_data.items():
                for spin in orbs_data.keys():
                    results[pair][orbs][spin] = Function1D(self.energies, orbs_data[spin]['single'])
        return results

    @lazy_property
    def functions_pair(self):
        """
        Extracts a dictionary with keys pair, spin and containing a |Function1D| object for the total COP.
        """
        results = tree()
        for pair, pair_data in self.total.items():
            for spin in pair_data.keys():
                results[pair][spin] = Function1D(self.energies, pair_data[spin]['single'])
        return results

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []; app = lines.append
        if verbose:
            #app(marquee("File Info", mark="="))
            app(self.filestat(as_string=True))
        app("")
        app("%s: Number of energies: %d, from %.3f to %.3f (eV) with E_fermi set 0 (was %.3f)" % (
            self.cop_type.upper(), len(self.energies), self.energies[0], self.energies[-1], self.fermie))
        app("has_partial_projections: %s, nsppol: %d" % (bool(self.partial), self.nsppol))
        app("Number of pairs: %d" % len(self.total))
        for i, pair in enumerate(self.total):
            type0, type1 = self.type_of_index[pair[0]], self.type_of_index[pair[1]]
            app("[%d] %s@%s --> %s@%s" % (i, type0, pair[0], type1, pair[1]))

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(what="d", show=False)
        yield self.plot(what="i", show=False)

    @add_fig_kwargs
    def plot(self, what="d", spin=None, ax=None, exchange_xy=False, fontsize=12, **kwargs):
        """
        Plot COXP averaged values (DOS or IDOS depending on what).

        Args:
            what: string selecting what will be plotted. "d" for DOS, "i" for IDOS
            spin: Select spin index if not None.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            exchange_xy: True to exchange x-y axis.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        if not self.averaged: return None
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        xlabel = r"$E - E_f\;(eV)$"
        if self.cop_type == "coop":
            ysign = +1
            ylabel = {"d": "COOP", "i": "ICOOP"}[what]
        elif self.cop_type == "cohp":
            ysign = -1
            ylabel = {"d": "-COHP", "i": "-ICOHP"}[what]
        else:
            raise ValueError("Wrong cop_type: `%s`" % str(self.cop_type))

        spins = range(self.nsppol) if spin is None else [spin]
        for spin in spins:
            opts = {"color": "black", "linewidth": 2.0, "ls": "-"} if spin == 0 else \
                   {"color": "red", "linewidth": 2.0, "ls": "--"}
            opts.update(kwargs)
            key = {"d": "single", "i": "integrated"}[what]
            xs, ys = self.energies, ysign * self.averaged[spin][key]
            if exchange_xy: xs, ys = ys, xs
            ax.plot(xs, ys, label=self.spin2tex[spin] if self.nsppol == 2 else None,
                    **opts)

        if exchange_xy:
            # Add vertical line to signal the zero.
            ax.axvline(c="k", ls=":", lw=1)
        else:
            ax.axhline(c="k", ls=":", lw=1)

        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)

        if self.nsppol == 2:
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_average_pairs(self, with_site_index, what="single", exchange_xy=False,
                            fontsize=8, **kwargs):
        """
        Plot COXP total overlap for all sites containg `with_site_index` and average sum
        (multiplied by the number of pairs)

        Args:
            with_site_index: int of list of integers selecting the site index to be included.
            what: "single" for COXP DOS, "integrated" for IDOS.
            exchange_xy: True to exchange x-y axis.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        if not duck.is_listlike(with_site_index): with_site_index = [with_site_index]

        # Create list of pairs.
        all_pairs = []
        for index in with_site_index:
            pairs = [p for p in self.total if index in p]
            all_pairs.append(pairs)
            if not pairs:
                cprint("Cannot find pairs containing  site index %s" % str(index), "yellow")
            else:
                print("Pairs containing site index:", index)
                for pair in pairs:
                    type0, type1 = self.type_of_index[pair[0]], self.type_of_index[pair[1]]
                    print("\t%s@%s --> %s@%s" % (type0, pair[0], type1, pair[1]))

        # Build (1, nsppol) grid
        #ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=1, ncols=self.nsppol,
                                                sharex=True, sharey=True, squeeze=False)

        xlabel, ylabel, ysign = r"$E - E_f\;(eV)$", "COOP", +1
        if self.cop_type == "cohp":
            ylabel, ysign = "-COHP", -1

        # self.total[pair][spin][what]
        for spin, ax in enumerate(ax_list.ravel()):
            ax.grid(True)
            sum_all = np.zeros_like(self.energies)
            for index, pairs in zip(with_site_index, all_pairs):
                if not pairs: continue
                ys = np.zeros_like(self.energies)
                for pair in pairs:
                    ys += self.total[pair][spin][what]
                #ys /= len(pairs)
                ys /= len(self.total)
                sum_all += ys

                #label, style = self.get_labelstyle_from_spin_pair(spin, pair)
                label = "Average over pairs with %s@%s" % (self.type_of_index[index], index)
                style = {}
                xs, ys = self.energies, ysign * ys
                if exchange_xy: xs, ys = ys, xs
                ax.plot(xs, ys, label=label, **style)

            # Plot average * num_pairs for this spin.
            xs, ys = self.energies, ysign * self.averaged[spin][what] * 2 * len(self.total)
            if exchange_xy: xs, ys = ys, xs
            style = {}
            ax.plot(xs, ys, label="Average", **style)

            # Compare with sum
            xs, ys = self.energies, ysign * sum_all
            if exchange_xy: xs, ys = ys, xs
            ax.plot(xs, ys, label="Sumall", **style)

            if exchange_xy:
                # Add vertical line to signal the zero.
                ax.axvline(c="k", ls=":", lw=1)
            else:
                ax.axhline(c="k", ls=":", lw=1)

            if self.nsppol == 2: ax.set_title(self.spin2tex[spin])
            set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_site_pairs_total(self, from_site_index, what="single", exchange_xy=False, ax=None,
                              fontsize=8, **kwargs):
        """
        Plot COXP total overlap (DOS or IDOS) for all sites listed in `from_site_index`

        Args:
            from_site_index: int of list of integers selecting the first site of the pairs to be included.
            what: "single" for COXP DOS, "integrated" for IDOS.
            exchange_xy: True to exchange x-y axis.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        # Handle list of sites.
        if duck.is_listlike(from_site_index):
            nrows = len(from_site_index)
            sharex, sharey = True, False
            if exchange_xy: sharex, sharey = sharey, sharex
            ax_list, fig, plt = get_axarray_fig_plt(ax, nrows=nrows, ncols=1,
                                                    sharex=sharex, sharey=sharey, squeeze=False)
            # Recursive call.
            for ax, index in zip(ax_list.ravel(), from_site_index):
                self.plot_site_pairs_total(index, what=what, exchange_xy=exchange_xy, ax=ax,
                    fontsize=fontsize, show=False)

            return fig

        # Single site
        pairs = [p for p in self.total if from_site_index == p[0]]
        if not pairs:
            cprint("Cannot find pairs starting from site index %s" % str(from_site_index), "yellow")
            return None

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        xlabel, ylabel, ysign = r"$E - E_f\;(eV)$", "COOP", +1
        if self.cop_type == "cohp":
            ylabel, ysign = "-COHP", -1

        # self.total[pair][spin][what]
        for pair in pairs:
            for spin in range(self.nsppol):
                xs, ys = self.energies, ysign * self.total[pair][spin][what]
                if exchange_xy: xs, ys = ys, xs
                label, style = self.get_labelstyle_from_spin_pair(spin, pair)
                ax.plot(xs, ys, label=label, **style)

        if exchange_xy:
            # Add vertical line to signal the zero.
            ax.axvline(c="k", ls=":", lw=1)
        else:
            ax.axhline(c="k", ls=":", lw=1)

        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def get_labelstyle_from_spin_pair(self, spin, pair):
        """Return label and linestyle for spin and pair indices"""
        type0, type1 = self.type_of_index[pair[0]], self.type_of_index[pair[1]]
        label = r"${%s}@{%s} \rightarrow {%s}@{%s}$" % (type0, pair[0], type1, pair[1])
        style = {"color": "black", "linewidth": 1.0} if spin == 0 else \
                {"color": "red", "linewidth": 1.0}
        # TODO: Improve style
        style = dict(lw=2, ls="-" if spin == 0 else "--") #, color=color)
        return label, style

    @add_fig_kwargs
    def plot_site_pairs_partial(self, from_site_index, what="single", exchange_xy=True, ax=None,
                                fontsize=8, **kwargs):
        """
        Plot partial crystal orbital projections (DOS or IDOS) for all sites listed in `from_site_index`

        Args:
            from_site_index: int of list of integers selecting the first site of the pairs to be included.
            what: "single" for COXP DOS, "integrated" for IDOS.
            exchange_xy: True to exchange x-y axis.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        if not self.partial:
            cprint("File does not contain partial projections", "yellow")
            return None

        # Handle list of sites.
        if duck.is_listlike(from_site_index):
            nrows, ncols = 1, len(from_site_index)
            sharex, sharey = True, False
            if exchange_xy: sharex, sharey = sharey, sharex
            ax_list, fig, plt = get_axarray_fig_plt(ax, nrows=nrows, ncols=ncols,
                                                    sharex=sharex, sharey=sharey, squeeze=False)
            # Recursive call.
            for ax, index in zip(ax_list.ravel(), from_site_index):
                self.plot_site_pairs_partial(index, what=what, exchange_xy=exchange_xy, ax=ax,
                    fontsize=fontsize, show=False)

            return fig

        # Single site.
        pairs = [p for p in self.partial if from_site_index == p[0]]
        if not pairs:
            cprint("Cannot find pairs starting from site index %s" % str(from_site_index), "yellow")
            return None

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        xlabel, ylabel, ysign = r"$E - E_f\;(eV)$", "pCOOP", +1
        if self.cop_type == "cohp":
            ylabel, ysign = "-COHP", -1

        # [(0, 1)]["4s", "4p_x"][spin]["single"]
        for pair in pairs:
            for orbs, d in self.partial[pair].items():
                for spin in range(self.nsppol):
                    xs, ys = self.energies, ysign * d[spin][what]
                    if exchange_xy: xs, ys = ys, xs
                    label, style = self.get_labelstyle_from_spin_pair_orbs(spin, pair, orbs)
                    ax.plot(xs, ys, label=label, **style)

        if exchange_xy:
            # Add vertical line to signal the zero.
            ax.axvline(c="k", ls=":", lw=2)
        else:
            ax.axhline(c="k", ls=":", lw=2)

        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def get_labelstyle_from_spin_pair_orbs(self, spin, pair, orbs):
        """Return label and linestyle for spin, pair indices and orbs tuple."""
        type0, type1 = self.type_of_index[pair[0]], self.type_of_index[pair[1]]
        label = r"$%s_{%s}@%s \rightarrow %s_{%s}@%s$" % (type0, orbs[0], pair[0], type1, orbs[1], pair[1])
        # TODO: Improve style
        #style = self.params_orbs[orbs[0]]["style"]
        style = dict(lw=2, ls="-" if spin == 0 else "--") #, color=color)
        return label, style

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("coxpfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(coxpfile)"),
            nbv.new_code_cell("coxpfile.plot(what='d');"),
            nbv.new_code_cell("coxpfile.plot(what='i');"),
            nbv.new_code_cell("coxpfile.plot_site_pairs_total(from_site_index=[0,]);"),
            nbv.new_code_cell("coxpfile.plot_site_pairs_partial(from_site_index=[0,]);"),
            nbv.new_code_cell("#coxpfile.plot_with_ebands(ebands='filepath');"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class ICoxpFile(_LobsterFile):
    """
    Wrapper class for the integrated crystal orbital projections up to the Fermi energy
    produced from Lobster.
    May contain the output stored in ICOHPLIST.lobster and ICOOPLIST.lobster

    .. attribute::  cop_type

        String. Either "coop" or "cohp".

    .. attribute:: values

        A dictionary with the following keys: a tuple with the index of the sites
        considered as a pair (0-based, e.g. (0,1)), the spin (i.e. 0 or 1)

    .. attribute:: type_of_index

        Dictionary mappping site index to element string.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance of ICoxpFile from the files produce by Lobster.
        Accepts gzipped files.

        Args:
            filepath: path to the ICOHPLIST.lobster or ICOOPLIST.lobster.

        Returns:
            A ICoxpFile.
        """
        float_patt = r'-?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
        header_patt = re.compile(r'.*?(over+\s#\s+bonds)?\s+for\s+spin\s+(\d).*')
        data_patt = re.compile(r'\s+\d+\s+([a-zA-Z]+)(\d+)\s+([a-zA-Z]+)(\d+)\s+('+
                               float_patt+r')\s+('+float_patt+r')(\d+)?')

        new = cls(filepath)
        new.values = tree()

        spin = None
        avg_num_bonds = False
        new.type_of_index = {}
        with zopen(filepath, "rt") as f:
            for line in f:
                match = header_patt.match(line.rstrip())
                if match:
                    spin = [0, 1][int(match.group(2))-1]
                    avg_num_bonds = match.group(1) is not None
                match = data_patt.match(line.rstrip())
                if match:
                    type1, index1, type2, index2, dist, avg, n_bonds = match.groups()
                    # 0-based indexing
                    index1 = int(index1) - 1
                    index2 = int(index2) - 1
                    new.type_of_index[index1] = type1
                    new.type_of_index[index2] = type2
                    avg_data = {'average': float(avg), 'distance': dist, 'n_bonds': int(n_bonds) if n_bonds else None}
                    new.values[(index1, index2)][spin] = avg_data
                    new.values[(index2, index1)][spin] = avg_data

        new.cop_type = "unknown"
        if "ICOOPLIST.lobster" in filepath: new.cop_type = "coop"
        if "ICOHPLIST.lobster" in filepath: new.cop_type = "cohp"

        return new

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []; app = lines.append
        app("Number of pairs: %d" % len(self.values))
        app(self.dataframe.to_string(index=False))

        return "\n".join(lines)

    @lazy_property
    def dataframe(self):
        """|pandas-DataFrame| with results."""
        # self.values[pair][spin]
        import pandas as pd
        rows = []
        for pair, d in self.values.items():
            for spin in sorted(d.keys()):
                rows.append(OrderedDict([
                    ("index0", pair[0]),
                    ("index1", pair[1]),
                    ("type0", self.type_of_index[pair[0]]),
                    ("type1", self.type_of_index[pair[1]]),
                    ("spin", spin),
                    ("average", d[spin]["average"]),
                    ("distance", d[spin]["distance"]),
                    ("n_bonds", d[spin]["n_bonds"]),
                    ("pair", (pair[0], pair[1])),
                ]))

        return pd.DataFrame(rows, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """Barplot with average values."""
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        import seaborn as sns
        sns.barplot(x="average", y="pair", hue="spin", data=self.dataframe, ax=ax)
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("icoxp_file = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(icoxp_file)"),
            nbv.new_code_cell("icoxp_file.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class LobsterDoscarFile(_LobsterFile):
    """
    Total and partial dos extracted from lobster DOSCAR.
    The fermi energy is always at the zero value.

    .. attribute:: energies

        List of energies. Shifted such that the Fermi level lies at 0 eV.

    .. attribute:: total_dos

        A dictionary with spin as a key (i.e. 0 or 1) and containing the values of the total DOS.
        Should have the same size as energies.

    .. attribute:: pdos

        A dictionary with the values of the projected DOS.
        The dictionary should have the following nested keys: the index of the site (0-based),
        the string representing the projected orbital (e.g. "4p_x"), the spin (i.e. 0 or 1).
        Each dictionary should contain a numpy array with a list of DOS values with the
        same size as energies.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance from the DOSCAR.lobster file.
        Accepts gzipped files.

        Args:
            filepath: path to the DOSCAR.lobster.

        Returns:
            A LobsterDoscarFile.
        """
        with zopen(filepath, "rt") as f:
            dos_data = f.readlines()

        new = cls(filepath)

        new.nsites = int(dos_data[0].split()[0])
        n_energies = int(dos_data[5].split()[2])
        new.fermie = float(dos_data[5].split()[3])

        n_spin = 1 if len(dos_data[6].split()) == 3 else 2
        spins = [0, 1][:n_spin]

        # extract np array for total dos
        tdos_data = np.fromiter((d for l in dos_data[6:6+n_energies] for d in l.split()),
                                dtype=np.float).reshape((n_energies, 1+2*n_spin))

        new.energies = tdos_data[:,0].copy()
        new.total_dos = {}
        for i_spin, spin in enumerate(spins):
            new.total_dos[spin] = tdos_data[:,1+2*i_spin].copy()

        new.pdos = tree()
        new.type_of_index = {}
        # read partial doses
        for i_site in range(new.nsites):
            i_first_line = 5+(n_energies+1)*(i_site+1)

            # read orbitals
            # 6.01503759     -14.03508772   401       2.29842595       1.00000000; Z= 31; 4s 4p_y 4p_z 4p_x
            tokens = dos_data[i_first_line].split(';')
            orbitals = tokens[-1].split()
            Z = int(tokens[-2].split()[-1])
            el = Element.from_Z(Z)
            new.type_of_index[i_site] = el.symbol

            # extract np array for partial dos
            pdos_data = np.fromiter((d for l in dos_data[i_first_line+1:i_first_line+1+n_energies] for d in l.split()),
                dtype=np.float).reshape((n_energies, 1+n_spin*len(orbitals)))

            for i_orb, orb in enumerate(orbitals):
                for i_spin, spin in enumerate(spins):
                    new.pdos[i_site][orb][spin] = pdos_data[:, i_spin+n_spin*i_orb+1]

        new.nsppol = len(new.total_dos)
        return new

    def to_string(self, verbose=0):
        """String representation with Verbosity level `verbose`."""
        lines = []; app = lines.append
        app("Number of energies: %d, from %.3f to %.3f (eV) with E_fermi set to 0 (was %.3f)" % (
            len(self.energies), self.energies[0], self.energies[-1], self.fermie))
        app("nsppol: %d" % (self.nsppol))
        app("Number of sites in projected DOS: %d" % len(self.pdos))
        for i_site, dsite in self.pdos.items():
            app("%d --> {%s}" % (i_site, ", ".join(dsite.keys())))

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        for site_index in range(self.nsites):
            yield self.plot_pdos_site(site_index, show=False)

    @add_fig_kwargs
    def plot(self, spin=None, ax=None, exchange_xy=False, fontsize=12, **kwargs):
        """
        Plot DOS.

        Args:
            spin:
            exchange_xy: True to exchange x-y axis.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        spins = range(self.nsppol) if spin is None else [spin]
        for spin in spins:
            opts = {"color": "black", "linewidth": 2.0, "ls": "-"} if spin == 0 else \
                   {"color": "red", "linewidth": 2.0, "ls": "--"}
            opts.update(kwargs)
            spin_sign = +1 if spin == 0 else -1
            xs, ys = self.energies, spin_sign * self.total_dos[spin]
            if exchange_xy: xs, ys = ys, xs
            ax.plot(xs, ys, label=self.spin2tex[spin] if self.nsppol == 2 else None,
                    **opts)

        xlabel, ylabel = r"$E - E_f\;(eV)$", "DOS"
        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)

        if self.nsppol == 2:
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_pdos_site(self, site_index, ax=None, exchange_xy=False, fontsize=8, **kwargs):
        """
        Plot projected DOS

        Args:
            site_index
            exchange_xy: True to exchange x-y axis.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        # Handle list of sites.
        if duck.is_listlike(site_index):
            nrows = len(site_index)
            sharex, sharey = True, False
            if exchange_xy: sharex, sharey = sharey, sharex
            ax_list, fig, plt = get_axarray_fig_plt(ax, nrows=nrows, ncols=1,
                                                    sharex=sharex, sharey=sharey, squeeze=False)
            # Recursive call.
            for ax, index in zip(ax_list.ravel(), site_index):
                self.plot_pdos_site(index, exchange_xy=exchange_xy, ax=ax, fontsize=fontsize, show=False)

            return fig

        # Single site.
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        # self.pdos[site_index]["4p_x"][spin]
        for orb, d in self.pdos[site_index].items():
            for spin in range(self.nsppol):
                spin_sign = +1 if spin == 0 else -1
                xs, ys = self.energies, spin_sign * d[spin]
                if exchange_xy: xs, ys = ys, xs
                n, lm = orb[0], orb[1:]
                params = self.params_orbs[lm]
                label = "$%s@%s_{%s}$" % (self.type_of_index[site_index], site_index, params["latex"])
                ax.plot(xs, ys, label=label, **params["style"])

        xlabel, ylabel = r"$E - E_f\;(eV)$", "PDOS"
        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)

        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("lobdos = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(lobdos)"),
            nbv.new_code_cell("lobdos.plot();"),
            nbv.new_code_cell("pdos.plot_pdos_site(site_index=[0,]);"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class LobsterInput(object):
    """
    This object stores the basic variables for a Lobster input and generates the lobsterin file.
    """

    accepted_basis_sets = {"bunge", "koga", "pbevaspfit2015"}

    available_advanced_options = {"basisRotation", "writeBasisFunctions", "onlyReadVasprun.xml", "noMemoryMappedFiles",
                                  "skipPAWOrthonormalityTest", "doNotIgnoreExcessiveBands", "doNotUseAbsoluteSpilling",
                                  "skipReOrthonormalization", "doNotOrthogonalizeBasis", "forceV1HMatrix",
                                  "noSymmetryCorrection", "symmetryDetectionPrecision", "useOriginalTetrahedronMethod",
                                  "useDecimalPlaces", "forceEnergyRange"}

    def __init__(self, basis_set=None, basis_functions=None, include_orbitals=None, atom_pairs=None, dist_range=None,
                 orbitalwise=True, start_en=None, end_en=None, en_steps=None, gaussian_smearing=None,
                 bwdf=None, advanced_options=None):
        """
        Args
            basis_set: String containing one of the possible basis sets available: bunge, koga, pbevaspfit2015
            basis_functions: list of strings giving the symbol of each atom and the basis functions: "Ga 4s 4p"
            include_orbitals: string containing which types of valence orbitals to use. E.g. "s p d"
            atom_pairs: list of tuples containing the couple of elements for which the COHP analysis will be
             performed. Index is 1-based.
            dist_range: list of tuples, each containing the minimum and maximum distance (in Angstrom) used to
             automatically generate atom pairs. Each tuple can also contain two atomic symbol to restric the match to
             the specified elements. examples: (0.5, 1.5) or (0.5, 1.5, 'Zn', 'O')
            start_en: starting energy with respect to the Fermi level (in eV)
            end_en: ending energy with respect to the Fermi level (in eV)
            en_steps: number of energy increments
            gaussian_smearing: smearing in eV when using gaussian broadening
            bwdf: enables the bond-weighted distribution function (BWDF). Value is the binning interval
            advanced_options: dict with additional advanced options. See lobster user guide for further details
        """
        if basis_set and basis_set.lower() not in self.accepted_basis_sets:
            raise ValueError("Wrong basis set {}".format(basis_set))
        self.basis_set = basis_set
        self.basis_functions = basis_functions or []
        self.include_orbitals = include_orbitals
        self.atom_pairs = atom_pairs or []
        self.dist_range = dist_range or []
        self.orbitalwise = orbitalwise
        self.start_en = start_en
        self.end_en = end_en
        self.en_steps = en_steps
        self.gaussian_smearing = gaussian_smearing
        self.bwdf = bwdf
        self.advanced_options = advanced_options or {}

        if not all(opt in self.available_advanced_options for opt in self.advanced_options.keys()):
            raise ValueError("Unknown adavanced options")

    @classmethod
    def _get_basis_functions_from_abinit_pseudos(cls, pseudos):
        """
        Extracts the basis function used from the PAW abinit pseudopotentials

        Args:
            pseudos: a list of Pseudos objects.
        """
        basis_functions = []
        for p in pseudos:
            if not hasattr(p, "valence_states"):
                raise RuntimeError("Only PAW pseudos in PAWXML format are supported by Lobster interface.")
            el = p.symbol + " ".join(str(vs['n'] + OrbitalType(int(vs['l'])).name)
                                     for vs in p.valence_states.values() if 'n' in vs)
            basis_functions.append(el)
        return basis_functions

    def set_basis_functions_from_abinit_pseudos(self, pseudos):
        """
        Sets the basis function used from the PAW abinit pseudopotentials

        Args:
            pseudos: a list of Pseudos objects.
        """
        basis_functions = self._get_basis_functions_from_abinit_pseudos(pseudos)

        self.basis_functions = basis_functions

    @classmethod
    def _get_basis_functions_from_potcar(cls, potcar):
        """
        Extracts the basis function used from a POTCAR.

        Args:
            potcar: a pymatgen.io.vasp.inputs.Potcar object
        """
        basis_functions = []
        for p in potcar:
            basis_functions.append(p.element +" "+ " ".join(str(vs[0]) + vs[1] for vs in p.electron_configuration))
        return basis_functions

    def set_basis_functions_from_potcar(self, potcar):
        """
        Sets the basis function used from a POTCAR.

        Args:
            potcar: a pymatgen.io.vasp.inputs.Potcar object
        """
        basis_functions = self._get_basis_functions_from_potcar(potcar)

        self.basis_functions = basis_functions

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []

        if self.basis_set:
            lines.append("basisSet {}".format(self.basis_set))

        for bf in self.basis_functions:
            lines.append("basisFunctions {}".format(bf))

        for ap in self.atom_pairs:
            line = "cohpBetween atom {} atom {}".format(*ap)
            if self.orbitalwise:
                line += " orbitalwise"
            lines.append(line)

        for dr in self.dist_range:
            line = "cohpGenerator from {} to {}".format(dr[0], dr[1])
            if len(dr) > 2:
                line += " type {} type {}".format(dr[2], dr[3])
            if self.orbitalwise:
                line += " orbitalwise"
            lines.append(line)

        if self.start_en:
            lines.append("COHPStartEnergy {}".format(self.start_en))

        if self.end_en:
            lines.append("COHPEndEnergy {}".format(self.end_en))

        if self.en_steps:
            lines.append("COHPSteps {}".format(self.en_steps))

        if self.gaussian_smearing:
            lines.append("gaussianSmearingWidth {}".format(self.gaussian_smearing))

        if self.bwdf:
            lines.append("BWDF {}".format(self.bwdf))

        for k, v in self.advanced_options.items():
            lines.append("{} {}".format(k, v))

        return "\n".join(lines)

    @classmethod
    def from_dir(cls, dirpath, dE=0.01, **kwargs):
        """
        Generates an instance of the class based on the output folder of a DFT calculation.
        Reads the information from the pseudopotentials in order to determine the
        basis functions.

        Args:
            dirpath: the path to the calculation directory. For abinit it should contain the
                "files" file and GSR file, for vasp it should contain the vasprun.xml and the POTCAR.
            dE: The spacing of the energy sampling in eV.
            kwargs: the inputs for the init method, except for basis_functions, start_en,
                end_en and en_steps.

        Returns:
            A LobsterInput.
        """
        # Try to determine the code used for the calculation
        dft_code = None
        if os.path.isfile(os.path.join(dirpath, 'vasprun.xml')):
            dft_code = "vasp"
            vr = Vasprun(os.path.join(dirpath, 'vasprun.xml'))

            en_min = np.min([bands_spin for bands_spin in vr.eigenvalues.values()])
            en_max = np.max([bands_spin for bands_spin in vr.eigenvalues.values()])
            fermie = vr.efermi

            potcar = Potcar.from_file(os.path.join(dirpath, 'POTCAR'))
            basis_functions = cls._get_basis_functions_from_potcar(potcar)

        elif glob.glob(os.path.join(dirpath, '*.files')):
            dft_code = "abinit"
            ff = glob.glob(os.path.join(dirpath, '*.files'))[0]
            with open(ff, "rt") as files_file:
                ff_lines = files_file.readlines()
            out_path = ff_lines[3].strip()
            if not os.path.isabs(out_path):
                out_path = os.path.join(dirpath, out_path)

            with GsrFile.from_file(out_path + '_GSR.nc') as gsr:
                en_min = gsr.ebands.eigens.min()
                en_max = gsr.ebands.eigens.max()
                fermie = gsr.ebands.fermie

            pseudo_paths = []
            for l in ff_lines[5:]:
                l = l.strip()
                if l:
                    if not os.path.isabs(l):
                        l = os.path.join(dirpath, l)
                    pseudo_paths.append(l)

            pseudos = [Pseudo.from_file(p) for p in pseudo_paths]

            basis_functions = cls._get_basis_functions_from_abinit_pseudos(pseudos)
        else:
            raise ValueError('Unable to determine the code used in dir {}'.format(dirpath))

        start_en = en_min + fermie
        end_en = en_max - fermie

        # shift the energies so that are divisible by dE and the value for the fermi level (0 eV) is included
        start_en = np.floor(start_en/dE)*dE
        end_en= np.ceil(end_en/dE)*dE

        en_steps = int((end_en-start_en)/dE)

        return cls(basis_functions=basis_functions, start_en=start_en, end_en=end_en, en_steps=en_steps, **kwargs)

    def write(self, dirpath='.'):
        """
        Write the input file 'lobsterin' in dirpath.
        """
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        # Write the input file.
        with open(os.path.join(dirpath, 'lobsterin'), "wt") as f:
            f.write(str(self))


class LobsterAnalyzer(NotebookWriter):

    @classmethod
    def from_dir(cls, dirpath, prefix=""):
        """
        Generates an instance of the class based on the output folder of a DFT calculation.

        Args:
            dirpath: the path to the calculation directory.
        """
        dirpath = os.path.abspath(dirpath)
        k2ext = {
            "coop_path": "COOPCAR.lobster",
            "cohp_path": "COHPCAR.lobster",
            "icohp_path": "ICOHPLIST.lobster",
            "lobdos_path": "DOSCAR.lobster",
        }
        kwargs = {k: None for k in k2ext}
        for k, ext in k2ext.items():
            # Handle [prefix]COHPCAR.lobster
            p = "%s*%s" % (prefix, ext)
            paths = glob.glob(os.path.join(dirpath, p))
            # Treat [prefix]COHPCAR.lobster.gz
            if not paths:
                p = "%s*%s.gz" % (prefix, ext)
                paths = glob.glob(os.path.join(dirpath, p))

            if paths:
                #print(k, "-->", paths)
                if len(paths) > 1:
                    raise RuntimeError("Found multiple files matching glob pattern: %s" % str(paths))
                kwargs[k] = paths[0]

        return cls(dirpath, prefix, **kwargs)

    def __init__(self, dirpath, prefix, coop_path=None, cohp_path=None, icohp_path=None, lobdos_path=None):
        self.coop = CoxpFile.from_file(coop_path) if coop_path else None
        self.cohp = CoxpFile.from_file(cohp_path) if cohp_path else None
        self.icohp = ICoxpFile.from_file(icohp_path) if icohp_path else None
        self.doscar = LobsterDoscarFile.from_file(lobdos_path) if lobdos_path else None
        self.dirpath = dirpath
        self.prefix = prefix

        for a in ("coop", "cohp", "icohp", "doscar"):
            obj = getattr(self, a)
            if obj is not None:
                self.nsppol = obj.nsppol
                break

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []; app = lines.append
        for aname, header in [("coop", "COOP File"), ("cohp", "COHP File"),
                              ("icohp", "ICHOHPLIST File"), ("doscar", "Lobster DOSCAR"),]:
            obj = getattr(self, aname, None)
            if obj is None: continue
            app(marquee(header, mark="="))
            app(obj.to_string(verbose=verbose))
            app("")

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)

    @add_fig_kwargs
    def plot(self, entries=("coop", "cohp", "doscar"), spin=None, **kwargs):
        """
        Plot COOP + COHP + DOSCAR.

        Args:
            entries: "cohp" to plot COHP, "coop" for COOP
            spin:

        Returns: |matplotlib-Figure|
        """
        entries = [e for e in entries if getattr(self, e, None)]
        if not entries: return None

        axmat, fig, plt = get_axarray_fig_plt(None, nrows=self.nsppol, ncols=len(entries),
                                              sharex=False, sharey=True, squeeze=False)

        spins = range(self.nsppol) if spin is None else [spin]
        for spin in spins:
            for ix, (label, ax) in enumerate(zip(entries, axmat[spin])):
                obj = getattr(self, label)
                obj.plot(ax=ax, spin=spin, exchange_xy=True, show=False)
                if ix != 0:
                    set_visible(ax, False, "ylabel")
                if self.nsppol == 2 and spin == 0:
                    set_visible(ax, False, "xlabel")

        return fig

    @add_fig_kwargs
    def plot_coxp_with_dos(self, from_site_index, what="cohp", with_orbitals=False, exchange_xy=True,
                           fontsize=8, **kwargs):
        """
        Plot COHP (COOP) for all sites in from_site_index and Lobster DOS on len(from_site_index) + 1 subfigures.

        Args:
            from_site_index: List of site indices for COHP (COOP)
            what: "cohp" to plot COHP, "coop" for COOP
            with_orbitals: True if orbital projections are wanted.
            exchange_xy: True to exchange x-y axis. By default, use x-axis for DOS-like quantities.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        from_site_index = from_site_index if duck.is_listlike(from_site_index) else [from_site_index]
        ncols = len(from_site_index) + 1
        sharex, sharey = True, False
        if exchange_xy: sharex, sharey = sharey, sharex
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=1, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()

        # Plot (COHP|COOP) (total|projections) on the first ncols - 1 axes.
        coxp_file = getattr(self, what)
        for ix, (ax, index) in enumerate(zip(ax_list[:-1], from_site_index)):
            if with_orbitals:
                coxp_file.plot_site_pairs_partial(index, exchange_xy=exchange_xy, ax=ax, fontsize=fontsize, show=False)
            else:
                coxp_file.plot_site_pairs_total(index, exchange_xy=exchange_xy, ax=ax, fontsize=fontsize, show=False)

            if ix != 0:
                set_visible(ax, False, "ylabel")
                #set_visible(ax, False, "xlabel")

        # Plot DOS on the last ax.
        ax = ax_list[-1]
        self.doscar.plot(ax=ax, exchange_xy=exchange_xy, fontsize=fontsize, show=False)
        set_visible(ax, False, "ylabel" if exchange_xy else "xlabel")

        return fig

    #@add_fig_kwargs
    #def plot_with_ebands(self, ebands, entries=("coop", "cohp", "doscar"), **kwargs):
    #    """
    #    Plot bands + COHP, COOP, DOSCAR.

    #    Args:
    #        ebands:
    #        entries: "cohp" to plot COHP, "coop" for COOP
    #        fontsize: fontsize for legends and titles

    #    Returns: |matplotlib-Figure|
    #    """
    #    ebands = ElectronBands.as_ebands(ebands)

    #    entries = [e for e in entries if getattr(self, e, None)]
    #    if not entries: return None

    #    import matplotlib.pyplot as plt
    #    from matplotlib.gridspec import GridSpec
    #    fig = plt.figure()

    #    # Build grid.
    #    nrows, ncols = self.nsppol, len(entries) + 1
    #    width_ratios = [2] + [1] * len(entries)
    #    gspec = GridSpec(nrows=nrows, ncols=ncols, width_ratios=width_ratios, wspace=0.05)

    #    # Bands and DOS will share the y-axis
    #    axmat = np.array((nrows, ncols), dtype=object)
    #    for icol in range(ncols):
    #        for irow in range(nrows):
    #            axmat[irow, icol] = plt.subplot(gspec[irow, icol], sharey=None if irow == 0 else axmat[irow, icol-1])

    #    #axmat, fig, plt = get_axarray_fig_plt(None, nrows=self.nsppol, ncols=len(entries) + 1,
    #    #                                        sharex=False, sharey=True, squeeze=False)

    #    for spin in range(self.nsppol):
    #        for ix, (label, ax) in enumerate(zip(entries, axmat[spin, :])):
    #            if ix == 0:
    #                ebands.plot(spin=spin, e0="fermie", ax=ax, show=False)
    #            else:
    #                obj = getattr(self, label)
    #                obj.plot(ax=ax, spin=spin, exchange_xy=True, show=False)

    #            if ix != 0:
    #                set_visible(ax, False, "ylabel")
    #            if self.nsppol == 2 and spin == 0:
    #                set_visible(ax, False, "xlabel")

    #    return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("lobana = abilab.LobsterAnalyzer.from_dir(dirpath='%s', prefix='%s')" % (
                self.dirpath, self.prefix)),
            nbv.new_code_cell("print(lobana)"),
            nbv.new_code_cell("lobana.plot();"),
            nbv.new_code_cell("lobana.plot_coxp_with_dos(from_site_index=[0,], what='cohp', with_orbitals=False);"),
            #nbv.new_code_cell("lobana.plot_with_ebands(ebands='path_to_file_with_ebands');"),
        ])

        return self._write_nb_nbpath(nb, nbpath)
