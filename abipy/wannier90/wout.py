# coding: utf-8
"""Interface to the wout output file produced by Wannier90."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd

from collections import OrderedDict
from monty.string import marquee
from abipy.core.mixins import BaseFile, Has_Structure, NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt


class WoutFile(BaseFile, Has_Structure, NotebookWriter):
    """
    Main output file produced by Wannier90

    Usage example:

    .. code-block:: python

        with abilab.abiopen("foo.wout") as wout:
            print(wout)
            wout.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: WoutFile
    """
    def __init__(self, filepath):
        super(WoutFile, self).__init__(filepath)
        self.warnings = []
        self.use_disentangle = False
        self.conv_df, self.dis_df = None, None

        with open(self.filepath, "rt") as fh:
            self.lines = fh.readlines()

        self._parse_dims()
        try:
            self._parse_iterations()
        except Exception as exc:
            print("Exception in _parse_iterations:\n", exc)

    def close(self):
        """Close file. Required by abc protocol."""

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app("Wannier90 version: %s" % self.version)
        app("Number of Wannier functions: %d" % self.nwan)
        app("K-grid: %s" % self.grid_size)
        if self.use_disentangle:
            app("Using DISENTANGLE algorithm")
            #for k, v in self.params_section["DISENTANGLE"].items():
            #    app("%s: %s" % (k, v))
        app("")

        if self.dis_df is not None:
            # Print first and last n cycles.
            app(marquee("DISENTANGLE", mark="="))
            n = 5 if not verbose else 20
            if len(self.dis_df) > 2 * n:
                app(pd.concat([self.dis_df.head(n), self.dis_df.tail(n)]).to_string(index=False))
            else:
                app(self.dis_df.to_string(index=False))
            app("")

        if self.conv_df is not None:
            # Print first and last n cycles.
            app(marquee("WANNIERISE", mark="="))
            n = 5 if not verbose else 20
            if len(self.conv_df) > 2 * n:
                app(pd.concat([self.conv_df.head(n), self.conv_df.tail(n)]).to_string(index=False))
            else:
                app(self.conv_df.to_string(index=False))
            app("")

        #if verbose:

        return "\n".join(lines)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def _parse_dims(self):
        """
        Parse basic dimensions and get structure from the header of the file.
        """
        self.version, self._structure, self.grid_size = None, None, None
        # Init dictionary with parameters.
        self.params_section = OrderedDict([(s, OrderedDict()) for s in
            ("MAIN", "WANNIERISE", "PLOTTING", "DISENTANGLE")])
        params_done = False

        for iln, line in enumerate(self.lines):
            # Check for any warnings
            if 'Warning' in line:
                self.warnings.append(line)
                continue

            if "Time to read parameters" in line:
                params_done = True
                continue

            # Get release string.
            if "Release:" in line:
                i = line.find("Release:")
                self.version = line[i:].split()[1]
                continue

            # Parse lattice.
            if "Lattice Vectors" in line and self._structure is None:
                #              Lattice Vectors (Ang)
                #    a_1     0.000000   2.715473   2.715473
                #    a_2     2.715473   0.000000   2.715473
                #    a_3     2.715473   2.715473   0.000000
                lattice = np.array([list(map(float, self.lines[iln+j].split()[1:])) for j in range(1, 4)])
                continue

            # Parse atoms.
            if "|   Site   " in line and self._structure is None:
                # *----------------------------------------------------------------------------*
                # |   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |
                # +----------------------------------------------------------------------------+
                # | Si   1   0.00000   0.00000   0.00000   |    0.00000   0.00000   0.00000    |
                # | Si   2   0.25000   0.25000   0.25000   |    1.35774   1.35774   1.35774    |
                # *----------------------------------------------------------------------------*
                frac_coords, species = [], []
                i = iln + 2
                while True:
                    l = self.lines[i].strip()
                    if l.startswith("*"): break
                    i += 1
                    tokens = l.replace("|", " ").split()
                    species.append(tokens[0])
                    frac_coords.append(np.array(list(map(float, tokens[2:5]))))

                self._structure = Structure(lattice, species, frac_coords)
                continue

            # Parse kmesh.
            if "Grid size" in line:
                # Grid size =  2 x  2 x  2      Total points =    8
                tokens = line.split("=")[1].split("Total")[0].split("x")
                self.grid_size = np.array(list(map(int, tokens)))
                continue

            if not params_done and any(sname in line for sname in self.params_section):
                #*---------------------------------- MAIN ------------------------------------*
                #|  Number of Wannier Functions               :                 4             |
                #|  Wavefunction spin channel                 :                up             |
                #*----------------------------------------------------------------------------*
                # Use params_done to avoid parsing the second section with WANNIERISE
                key = line.replace("*", "").replace("-", "").strip()
                i = iln + 1
                l = self.lines[i].strip()
                while not l.startswith("*-"):
                    tokens = [s.strip() for s in l.replace("|", "").split(":")]
                    self.params_section[key][tokens[0]] = tokens[1]
                    i += 1
                    l = self.lines[i].strip()
                continue

        # Extract important metadata from sections and convert from string.
        self.nwan = int(self.params_section["MAIN"]["Number of Wannier Functions"])
        if self.params_section["DISENTANGLE"].get("Using band disentanglement", "F") == "T":
            self.use_disentangle = True

    def _parse_iterations(self):
        """
        Parse iteration steps if not already done and store results in self.

        Return: 0 if success.
        """
        # Don't parse it again if already done.
        if self.conv_df is not None: return 0

        if self.use_disentangle:
            # Parse Disentanglement cycles
            # +---------------------------------------------------------------------+<-- DIS
            # |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
            # +---------------------------------------------------------------------+<-- DIS
            #       1       3.91743302       3.66269149       6.955E-02      0.28    <-- DIS
            #       2       3.66269149       3.66269149       2.021E-14      0.29    <-- DIS
            # <<<      Delta < 1.000E-10  over  3 iterations     >>>
            # <<< Disentanglement convergence criteria satisfied >>>
            in_diis = 0
            data = OrderedDict([(s, []) for s in
                ("iter", "omegaI_im1", "omegaI_i", "delta_frac", "time")])

            for line in self.lines:
                line = line.strip()
                if not line.endswith("<-- DIS"):
                    if in_diis: break
                    continue
                in_diis += 1
                if in_diis >= 4:
                    toks = line.split()
                    for i, key in enumerate(data.keys()):
                        data[key].append(int(toks[i]) if i == 0 else float(toks[i]))

            self.dis_df = pd.DataFrame.from_dict(data)

        # Parse Wannierization cycles.
        for start, line in enumerate(self.lines):
            if "Initial State" in line: break
        else:
            return 1

        self.wf_centers = [[] for _ in range(self.nwan)]
        self.wf_spreads = [[] for _ in range(self.nwan)]
        data = OrderedDict([(s, []) for s in
            ("iter",  "delta_spread", "rms_gradient", "spread", "time",
             "O_D", "O_OD", "O_TOT",
            )])

        lines = self.lines[start:][:]
        while lines:
            line = lines.pop(0).strip()
            if not line.startswith("Initial State") and not line.startswith("Cycle:"): continue
            step = int(line.split()[-1]) if line.startswith("Cycle:") else 0
            for iw in range(self.nwan + 1):
                # WF centre and spread    1  (  0.042127,  0.071712, -0.424794 )    10.42287858
                # Sum of centres and spreads (  0.933074, -0.071343, -0.800933 )    42.11245002
                tokens = lines.pop(0).replace("(", " ").replace(")", " ").replace(",", "").split()
                spread = float(tokens[-1])
                center = list(map(float, tokens[-4:-1]))
                if iw != self.nwan:
                    self.wf_centers[iw].append(center)
                    self.wf_spreads[iw].append(spread)

            line = lines.pop(0)
            assert not line.strip()

            #      0     0.421E+02     0.0000000000       42.1124500153       0.57  <-- CONV
            #        O_D=     32.4608805 O_OD=      5.9016636 O_TOT=     42.1124500 <-- SPRD
            # ------------------------------------------------------------------------------
            # Cycle:      1
            #  WF centre and spread    1  ( -0.141379, -0.258009, -0.488755 )     8.73529609
            #  Sum of centres and spreads (  0.794432, -0.178304, -0.823458 )    32.08361376
            #      1    -0.100E+02    10.2774612971       32.0836137609       0.57  <-- CONV
            #        O_D=     22.6617600 O_OD=      5.6719478 O_TOT=     32.0836138 <-- SPRD
            # Delta: O_D= -0.9799120E+01 O_OD= -0.2297158E+00 O_TOT= -0.1002884E+02 <-- DLTA
            # ------------------------------------------------------------------------------
            # Parse CONV and add values to data
            conv_toks = lines.pop(0).split()
            for i, k in enumerate(("iter", "delta_spread", "rms_gradient", "spread", "time")):
                data[k].append(int(conv_toks[i]) if i == 0 else float(conv_toks[i]))

            sprd_toks = lines.pop(0).split()
            data["O_D"].append(float(sprd_toks[1]))
            data["O_OD"].append(float(sprd_toks[3]))
            data["O_TOT"].append(float(sprd_toks[5]))
            if step > 0:
                dlta_tokens = lines.pop(0).split()

        self.conv_df = pd.DataFrame.from_dict(data)

        # Convert to numpy array (nwan, nstep, 3) and (nwan, nstep)
        self.wf_centers = np.array(self.wf_centers)
        self.wf_spreads = np.array(self.wf_spreads)

        return 0

    @add_fig_kwargs
    def plot(self, fontsize=12, **kwargs):
        """
        Plot the convergence of the Wannierise cycle.

        Args:
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        if self._parse_iterations() != 0:
            print("Wout files does not contain Wannierization cycles. Returning None")
            return None

        items = ["delta_spread", "rms_gradient", "spread"]
        if self.use_disentangle:
            items += ["omegaI_i"]

        # Build grid of plots.
        num_plots, ncols, nrows = len(items), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        marker = "."
        for ax, item in zip(ax_list, items):
            ax.grid(True)
            ax.set_xlabel("Iteration Step")
            ax.set_ylabel(item)
            s = 1
            if item == "omegaI_i":
                # Plot Disentanglement cycles
                ax.plot(self.dis_df.iter[s:], self.dis_df[item][s:], marker=marker)
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                ax2 = inset_axes(ax, width="60%", height="40%", loc="upper right")
                ax2.grid(True)
                ax2.set_title("delta_frac", fontsize=8)
                ax2.plot(self.dis_df.iter[s:], self.dis_df["delta_frac"][s:], marker=marker)

            else:
                ax.plot(self.conv_df.iter[s:], self.conv_df[item][s:], marker=marker)

        return fig

    @add_fig_kwargs
    def plot_centers_spread(self, fontsize=8, **kwargs):
        """
        Plot the convergence of the Wannier centers and spread
        as function of iteration number

        Args:
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        if self._parse_iterations() != 0:
            print("Wout files does not contain Wannierization cycles. Returning None")
            return None

        # Build grid of plots.
        # nwan subplot with evolution of the WF center + last subplot with all spreads
        num_plots, ncols, nrows = self.nwan + 1, 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        marker = "."
        for iax in range(num_plots):
            ax = ax_list[iax]
            ax.grid(True)
            ax.set_xlabel("Iteration Step")
            s = 1
            if iax < self.nwan:
                ax.set_ylabel("Center of WF #%s" % (iax + 1))
                for idir in range(3):
                    ax.plot(self.conv_df.iter[s:], self.wf_centers[iax, s:, idir], marker=marker,
                            label={0: "x", 1: "y", 2: "z"}[idir] if iax == 0 else None)
            else:
                ax.set_ylabel("WF Spread")
                for iw in range(self.nwan):
                    ax.plot(self.conv_df.iter[s:], self.wf_spreads[iw, s:], marker=marker,
                            label="WF#%d" % (iw + 1))

            if iax in (0, self.nwan):
                ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.plot_centers_spread(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("wout = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(wout)"),
            nbv.new_code_cell("wout.structure.plot();"),
            nbv.new_code_cell("wout.plot();"),
            nbv.new_code_cell("wout.plot_centers_spread();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)
