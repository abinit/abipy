# coding: utf-8
"""
Objects related to the computation of Debye-Waller tensors from the generalized phonon DOS.
"""
import numpy as np
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import list_strings, marquee
from monty.collections import dict2namedtuple
from monty.termcolor import cprint
from abipy.core.mixins import Has_Structure
from abipy.tools.plotting import add_fig_kwargs, set_axlims, get_axarray_fig_plt, set_visible
from abipy.tools.printing import print_dataframe


class _Component(object):
    """
    Object used to select/plot the components of the DW tensor.
    """

    def __init__(self, name, ij, **plot_kwargs):
        self.name = name
        self.ij = ij
        self.i, self.j = None, None
        if self.ij is not None:
            self.i, self.j = self.ij[0], self.ij[1]
        self.plot_kwargs = plot_kwargs

    def get_tavg_label(self, what, with_units=False):
        n = dict(displ="U", vel="V")[what]
        unit = ""
        if with_units:
            unit = r"\;%s" % (dict(displ=r"\AA^2", vel="v")[what])

        if self.name == "trace":
            return r"$\langle {%s}^2 \rangle%s$" % (n, unit)
        else:
            return r"$\langle {%s}_{%s} \rangle%s$" % (n, self.name, unit)

    def eval33w(self, mat33w):
        #assert mat33w.shape[:2] == (3, 3)
        if self.ij is not None:
            return mat33w[self.i, self.j]
        if self.name == "trace":
            return mat33w.trace()

        raise TypeError("Don't know how to extract data for `%s`" % self.name)


class MsqDos(Has_Structure):
    """
    This object stores the generalized phonon DOS in CARTESIAN coords.
    This DOS-like quantity allows one to calculate Debye Waller factors as a function of Temperature
    by integration with 1/omega and the Bose-Einstein factor.
    This object is usually instanciated via the read_msq_dos method of the PhdosReader
    and stored as attribute of the PhdosFile instance.

    See :cite:`Lee1995` for the further details about the internal implementation and
    :cite:`Trueblood1996` for the different conventions used by crystallographers.

    See also http://atztogo.github.io/phonopy/formulation.html#thermal-displacement
    """
    C = _Component
    ALL_COMPS = OrderedDict([
        ("trace", C(name="trace", ij=None, color="k")),
        ("xx", C(name="xx", ij=(0, 0), color="r", ls="-")),
        ("yy", C(name="yy", ij=(1, 1), color="g", ls="-")),
        ("zz", C(name="zz", ij=(2, 2), color="b", ls="-")),
        ("xy", C(name="xy", ij=(0, 1), color="c", ls="--")),
        ("xz", C(name="xz", ij=(0, 2), color="m", ls="--")),
        ("yz", C(name="yz", ij=(1, 2), color="y", ls="--")),
        # Symmetric components.
        #("yx", (1, 0)),
        #("zx", (2, 0)),
        #("zy", (2, 1)),
    ])
    del C

    def __init__(self, structure, wmesh, values, amu_symbol):
        """
        Arg:
            structure: |Structure| object.
            wmesh: Frequency mesh in eV
            values: (natom, 3, 3, nomega) array with generalized DOS in Cartesian coordinates.
            amu_symbol: Dictionary element.symbol -> mass in atomic units.
        """
        self._structure = structure
        self.wmesh = wmesh * abu.eV_Ha ####
        self.nw = len(self.wmesh)
        self.values = values * abu.Ha_eV ###
        self.amu_symbol = amu_symbol
        assert len(self.values) == len(self.structure)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def __str__(self):
        """Invoked by str"""
        return self.to_string()

    def to_string(self, verbose=0):
        """
        Human-readable string with useful information on the object.

        Args:
            verbose: Verbosity level.
        """
        lines = []; app = lines.append

        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(marquee(r"Fullfilment of \int dw g_ij(w) = \delta_ij", mark="="))
        app("")
        from scipy.integrate import simps
        for iatom, site in enumerate(self.structure):
            d = simps(self.values[iatom], x=self.wmesh)
            app("For site: %s" % site)
            app(str(d))
            app("Trace: %.4f, determinant: %.4f" % (d.trace(), np.linalg.det(d)))
            app("")

        for fmt in ("cartesian", "cif"):
            df = self.get_dataframe(temp=300, view="inequivalent", fmt=fmt, verbose=verbose)
            s = print_dataframe(df, title="Format: %s" % fmt, file="string")
            lines.extend(s.splitlines())

        title = marquee("Constraints on tensor components in reduced coords induced by site symmetries", mark="=")
        s = print_dataframe(self.structure.site_symmetries.get_tensor_rank2_dataframe(), file="string", title=title)
        lines.extend(s.splitlines())

        #if verbose > 1:
            #max_err = self.check_site_symmetries(verbose=verbose)

        return "\n".join(lines)

    def get_json_doc(self, tstart=0, tstop=600, num=11):
        """
        Return dictionary with results. Used by emmet builder.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate.
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        # (natom, 3, 3, nt)
        ucart_t = self.get_msq_tmesh(tmesh, what_list=("displ")).displ

        # Convert tensor to U_CIF format
        ucif_t = np.empty_like(ucart_t)
        for it in range(num):
            ucif_t[:,:,:,it] = self.convert_ucart(ucart_t[:,:,:,it], fmt="cif")

        jdoc = {
            "natom": len(self.structure),
            "nomega": self.nw,             # Number of frequencies
            "ntemp": len(tmesh),            # Number of temperatures
            "tmesh": tmesh,                 # Temperature mesh in K
            "wmesh": self.wmesh,            # Frequency mesh in ??
            "gdos_aijw": self.values,       # Generalized DOS in Cartesian coords. (natom, 3, 3, nomega) array.
            "amu_symbol": self.amu_symbol,  # Dict symbol --> Atomic mass in a.u.
            "structure": self.structure,    # Structure object
            "ucif_t": ucif_t,               # U tensors (natom, 3, 3, ntemp)  as a function of T for T in tmesh in CIF format.
            "ucif_string_t300k": self.get_cif_string(temp=300, symprec=None),  # String with U tensor at T=300 in Cif format.
        }

        return jdoc

    def get_msq_tmesh(self, tmesh, iatom_list=None, what_list=("displ", "vel")):
        """
        Compute mean square displacement for each atom in `iatom_list` as a function of T.
        in Cartesian coordinates and atomic-units.

        Args:
            tmesh: array-like with temperatures in Kelvin.
            iatom_list: List of atom sites to compute. None if all aomts are wanted.
            what_list: "displ" for displacement, "vel" for velocity tensor.

        Return:
            namedtuple with (tmesh=tmesh, displ=msq_d, vel=msq_v)

            msq_d = np.empty((natom, 3, 3, nt))
        """
        tmesh = np.array(tmesh)
        nt = len(tmesh)

        # Frequency mesh starts at iomin to avoid 1/0 and ignore eventual negative frequencies.
        for iomin, w in enumerate(self.wmesh):
            if w > 1e-12: break
        else:
            raise ValueError("Cannot find index such that wmesh[i] > 1e-12 !!!")
        wvals = self.wmesh[iomin:]
        nw = len(wvals)

        # We will compute: Ucart(T, k, ij) = 1/M_k \int dw (n(w) + 1/2) g_ij(w) / w for the k-atom in a.u.
        # Calculate Bose-Einstein occupation factors only once for each T (instead of for each atom).
        npht = np.zeros((nt, nw))
        for it, temp in enumerate(tmesh):
            npht[it] = abu.occ_be(wvals, temp * abu.kb_HaK) + 0.5

        natom = len(self.structure)
        msq_d = np.empty((natom, 3, 3, nt))
        msq_v = np.empty((natom, 3, 3, nt))
        what_list = list_strings(what_list)

        # Perform frequency integration to get tensor(T)
        from scipy.integrate import simps
        if iatom_list is not None: iatom_list = set(iatom_list)
        for iatom in range(natom):
            if iatom_list is not None and iatom not in iatom_list: continue
            symbol = self.structure[iatom].specie.symbol
            for it in range(nt):
                fn = self.values[iatom, :, :, iomin:] * npht[it]
                if "displ" in what_list:
                    # Mean square displacement for each atom as a function of T (bohr^2).
                    ys = fn / wvals
                    fact = 1.0 / (self.amu_symbol[symbol] * abu.amu_emass)
                    msq_d[iatom, :, :, it] = simps(ys, x=wvals) * fact * abu.Bohr_Ang ** 2
                if "vel" in what_list:
                    # Mean square velocity for each atom as a function of T (bohr^2/atomic time unit^2)"
                    ys = fn * wvals
                    fact = 1.0 / (self.amu_symbol[symbol] * abu.amu_emass)
                    msq_v[iatom, :, :, it] = simps(ys, x=wvals) * fact # * abu.velocity_at_to_si ** 2

        return dict2namedtuple(tmesh=tmesh, displ=msq_d, vel=msq_v)

    def convert_ucart(self, ucart_mat, fmt):
        """
        Convert the U tensor from Cartesian coordinates to format `fmt`
        Return new tensor. See also :cite:`Grosse-Kunstleve2002`.

        Args:
            ucart_mat: (natom,3,3) array with tensor in Cartesian coords.
            fmt: Output format. Available options: "cif", "ustar", "beta", "cartesian"

        Return: (natom, 3, 3) tensor.
        """
        natom = len(self.structure)
        if fmt == "cartesian":
            return ucart_mat.copy()

        #elif fmt == "B":
        #    # Eq 7
        #    return ucart_mat * 8 * np.pi**2

        elif fmt in ("cif", "ustar", "beta"):
            # Build A matrix
            amat = self.structure.lattice.matrix.T
            ainv = np.linalg.inv(amat)
            new_mat = np.zeros_like(ucart_mat)
            # Eq 3b: A^-1 U_Cart A^-T
            for iatom in range(natom):
                new_mat[iatom] = np.matmul(ainv, np.matmul(ucart_mat[iatom], ainv.T))

            # Now we have Ustar
            if fmt == "ustar": return new_mat
            if fmt == "beta": return new_mat * 8 * np.pi**2  # Eq 6

            # CIF Format Eq 4a
            # Build N matrix (no 2 pi factor)
            lengths = self.structure.lattice.reciprocal_lattice_crystallographic.lengths
            ninv = np.diag(1.0 / np.array(lengths, dtype=float))
            # N^-1 U_star N^-T
            for iatom in range(natom):
                new_mat[iatom] = np.matmul(ninv, np.matmul(new_mat[iatom], ninv.T))

            return new_mat

        raise ValueError("Invalid format: `%s`" % str(fmt))

    def get_dataframe(self, temp=300, fmt="cartesian", view="inequivalent", what="displ", decimals=4,
                      select_symbols=None, verbose=0):
        """
        Return |pandas-DataFrame| with cartesian tensor components as columns and (inequivalent) sites along the rows.

        Args:
            temp: Temperature in Kelvin.
            fmt: "cartesian" for elements in Cartesian coordinates, "cif" for results in reduced coordinates
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            what: "displ" for displament, "vel" for velocity.
            decimals: Number of decimal places to round to.
                If decimals is negative, it specifies the number of positions to the left of the decimal point.
            select_symbols: String or list of strings with chemical symbols. Used to select only atoms of this type.
            verbose: Verbosity level.

        Return: |pandas-DataFrame|
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols=select_symbols, verbose=verbose)

        # [natom, 3, 3, nt=1]
        msq = self.get_msq_tmesh([float(temp)], iatom_list=aview.iatom_list, what_list=what)
        ucart = getattr(msq, what)
        natom = len(self.structure)
        ucart = np.reshape(ucart, (natom, 3, 3))
        values = ucart
        if what == "displ":
            values = self.convert_ucart(ucart, fmt)

        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        inds = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
        rows = []
        for (iatom, wlabel) in zip(aview.iatom_list, aview.wyck_labels):
            site = self.structure[iatom]
            d = OrderedDict()
            d["element"] = site.specie.symbol
            d["site_index"] = iatom
            d["frac_coords"] = np.round(site.frac_coords, decimals=decimals)
            d["cart_coords"] = np.round(site.coords, decimals=decimals)
            d["wyckoff"] = wlabel
            if fmt == "cartesian":
                d["isotropic"] = ucart[iatom].trace() / 3.0
                d["determinant"] = np.linalg.det(ucart[iatom])
            for col, ind in zip(columns, inds):
                d[col] = values[iatom, ind[0], ind[1]]
            rows.append(d)

        import pandas as pd
        return pd.DataFrame(rows, columns=list(rows[0].keys()) if rows else None)

    def write_cif_file(self, filepath, temp=300, symprec=None):
        """
        Write CIF file with structure and anisotropic U tensor in CIF format.

        Args:
            filepath: Name of CIF file. If None, a temporary filepath is created.
            temp: Temperature in Kelvin (used to compute U).
            symprec (float): If not none, finds the symmetry of the structure
                and writes the cif with symmetry information. Passes symprec to the SpacegroupAnalyzer

        Return: Filepath
        """
        if filepath is None:
            import tempfile
            _, filepath = tempfile.mkstemp(suffix=".cif", text=True)

        with open(filepath, "wt") as fh:
            fh.write(self.get_cif_string(temp=temp, symprec=symprec))

        return filepath

    #def jsmol(self, temp=300, symprec=None, verbose=0): # pragma: no cover
    #    """
    #    Args:
    #        symprec (float): If not none, finds the symmetry of the structure
    #            and writes the cif with symmetry information. Passes symprec to the SpacegroupAnalyzer
    #        verbose: Verbosity level.
    #    """
    #    cif_string = self.get_cif_string(temp=temp, symprec=symprec)

    #    try:
    #      from jupyter_jsmol import JsmolView
    #    except ImportError:
    #        raise ImportError("jupyter_jsmol is not installed. See https://github.com/fekad/jupyter-jsmol")

    #    jsmol = JsmolView(color='white')
    #    from IPython.display import display, HTML
    #    FIXME TEMPORARY HACK
    #    display(HTML('<script type="text/javascript" src="/nbextensions/jupyter-jsmol/jsmol/JSmol.min.js"></script>'))
    #    display(jsmol)

    #    cmd = 'load inline "%s" {1 1 1}; ellipsoid;' % cif_string
    #    if verbose: print("executing cmd:", cmd)
    #    jsmol.script(cmd)

    #    return jsmol

    def vesta_open(self, temp=300): # pragma: no cover
        """
        Visualize termal displacement ellipsoids at temperature `temp` (Kelvin) with Vesta_ application.
        """
        filepath = self.write_cif_file(filepath=None, temp=temp)
        cprint("Writing structure + Debye-Waller tensor in CIF format for T = %s (K) to file: %s" % (temp, filepath), "green")
        cprint("In the Vesta GUI, select: Properties -> Atoms -> Show as displament ellipsoids.", "green")
        from abipy.iotools import Visualizer
        visu = Visualizer.from_name("vesta")

        return visu(filepath)()

    def get_cif_string(self, temp=300, symprec=None):
        """
        Return string with structure and anisotropic U tensor in CIF format at temperature `temp` in Kelvin

        Args:
            symprec (float): If not none, finds the symmetry of the structure
                and writes the cif with symmetry information. Passes symprec
                to the SpacegroupAnalyzer
        """
        # Get string with structure in CIF format.
        # Don't use symprec because it changes the order of the sites
        # and we must be consistent with site_labels when writing aniso_U terms!
        from pymatgen.io.cif import CifWriter
        cif = CifWriter(self.structure, symprec=symprec)

        aniso_u = """loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_23
_atom_site_aniso_U_13
_atom_site_aniso_U_12""".splitlines()

        # Compute U tensor in CIF format (reduced coords)
        natom = len(self.structure)
        msq = self.get_msq_tmesh([float(temp)], what_list="displ")
        ucart = getattr(msq, "displ")
        ucart = np.reshape(ucart, (natom, 3, 3))
        ucif = self.convert_ucart(ucart, fmt="cif")

        # Add matrix elements. Use 0 based index
        for iatom, site in enumerate(self.structure):
            site_label = "%s%d" % (site.specie.symbol, iatom)
            m = ucif[iatom]
            aniso_u.append("%s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" %
                    (site_label, m[0, 0], m[1, 1], m[2, 2], m[1, 2], m[0, 2], m[0, 1]))

        return str(cif) + "\n".join(aniso_u)

    def check_site_symmetries(self, temp=300, verbose=0):
        """
        Check site symmetries of the displacement tensor at temperature `temp` in Kelvin.
        Return maximum error.
        """
        msq = self.get_msq_tmesh([float(temp)], what_list="displ")
        ucart = getattr(msq, "displ")
        natom = len(self.structure)
        ucart = np.reshape(ucart, (natom, 3, 3))

        return self.structure.site_symmetries.check_site_symmetries(ucart, verbose=verbose)

    def _get_components(self, components):
        """
        Return list of components to analyze from user input.
        """
        if components == "all":
            return list(self.ALL_COMPS.values())
        elif components == "upper":
            return [self.ALL_COMPS[c] for c in ("xx", "yy", "zz", "yz", "xz", "xy")]
        elif components == "diag":
            return [self.ALL_COMPS[c] for c in ("xx", "yy", "zz")]
        elif components == "offdiag":
            return [self.ALL_COMPS[c] for c in ("xy", "xz", "yz")]
        else:
            return [self.ALL_COMPS[c] for c in list_strings(components)]

    @add_fig_kwargs
    def plot(self, components="upper", view="inequivalent", units="eV", select_symbols=None,
             xlims=None, ylims=None, sharey=False, fontsize=8, verbose=0, **kwargs):
        """
        Plot the generalized phonon DOS g_ij(omega, atom) for each atom in the unit cell.
        One subplot per atom. Each subplot shows the 9 independent components of the symmetric tensor.
        as a function of frequency. By default, only "inequivalent" atoms are shown.

        Args:
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            components: List of cartesian tensor components to plot e.g. ["xx", "xy"].
                "all" for all components. "upper" for the upper triangle, "diag" for diagonal elements.
            units: Units energy axis. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            select_symbols: String or list of strings with chemical symbols. Used to select only atoms of this type.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Set the data limits for the y-axis.
            sharey: True if y-axis should be shared.
            fontsize: Legend and title fontsize.
            verbose: Verbosity level.

        Returns: |matplotlib-Figure|
        """
        # TODO Decide units for internal arrays.
        factor = abu.phfactor_ev2units(units)

        # Select atoms.
        aview = self._get_atomview(view, select_symbols, verbose=verbose)

        num_plots = len(aview.iatom_list)
        nrows, ncols = 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=sharey, squeeze=True)
        ax_list = np.reshape(ax_list, (nrows, ncols)).ravel()
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        xx = self.wmesh * factor
        components = self._get_components(components)

        # For each atom in the view.
        for ix, (ax, iatom, site_label) in enumerate(zip(ax_list, aview.iatom_list, aview.site_labels)):
            irow, icol = divmod(ix, ncols)
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            ax.set_title(site_label, fontsize=fontsize)
            #site = self.structure[iatom]
            #color = cmap(float(iatom) / max((len(iatom_list) - 1), 1))

            # Plot components for this atom on the same ax.
            for c in components:
                yw = c.eval33w(self.values[iatom])
                label = r"$G_{%s}$" % c.name if ix == 0 else None
                ax.plot(xx, yw / factor, label=label, **c.plot_kwargs)

            # Handle labels.
            if irow == nrows - 1:
                ax.set_xlabel('Frequency %s' % abu.phunit_tag(units))
            else:
                set_visible(ax, False, "xlabel", "xticklabels")

            if ix == 0:
                ax.set_ylabel(r"$g_{ij}(\omega)$ 1/%s (Cart coords)" % abu.phunit_tag(units))
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_tensor(self, tstart=0, tstop=600, num=50, components="all", what="displ", view="inequivalent",
                    select_symbols=None, colormap="jet", xlims=None, ylims=None, fontsize=10, verbose=0, **kwargs):
        """
        Plot tensor(T) for each atom in the unit cell.
        One subplot for each component, each subplot show all inequivalent sites.
        By default, only "inequivalent" atoms are shown.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate.
            components: "all" for all components. "diag" for diagonal elements, "offdiag" for off-diagonal terms only.
            what: "displ" for displament, "vel" for velocity.
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            select_symbols: String or list of strings with chemical symbols. Used to select only atoms of this type.
            colormap: matplotlib colormap.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize.
            verbose: Verbosity level.

        Returns: |matplotlib-Figure|
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols=select_symbols, verbose=verbose)

        # One subplot for each component
        diag = ["xx", "yy", "zz"]
        offdiag = ["xy", "xz", "yz"]
        components = {
            "all": diag + offdiag, "diag": diag, "offdiag": offdiag,
        }[components]

        components = self._get_components(components)
        shape = np.reshape(components, (-1, 3)).shape
        nrows, ncols = shape[0], shape[1]

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=True)
        ax_list = np.reshape(ax_list, (nrows, ncols)).ravel()
        cmap = plt.get_cmap(colormap)

        # Compute U(T)
        tmesh = np.linspace(tstart, tstop, num=num)
        msq = self.get_msq_tmesh(tmesh, iatom_list=aview.iatom_list, what_list=what)
        # [natom,3,3,nt] array
        values = getattr(msq, what)

        for ix, (ax, comp) in enumerate(zip(ax_list, components)):
            irow, icol = divmod(ix, ncols)
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            ylabel = comp.get_tavg_label(what, with_units=True)
            ax.set_ylabel(ylabel, fontsize=fontsize)

            # Plot this component for all inequivalent atoms on the same subplot.
            for ii, (iatom, site_label) in enumerate(zip(aview.iatom_list, aview.site_labels)):
                color = cmap(float(ii) / max((len(aview.iatom_list) - 1), 1))
                ys = comp.eval33w(values[iatom])
                ax.plot(msq.tmesh, ys, label=site_label if ix == 0 else None,
                        color=color) #, marker="o")
                if ix == 0:
                    ax.legend(loc="best", fontsize=fontsize, shadow=True)

            if irow == 1:
                ax.set_xlabel('Temperature (K)')
            else:
                set_visible(ax, False, "xlabel", "xticklabels")

        return fig

    @add_fig_kwargs
    def plot_uiso(self, tstart=0, tstop=600, num=50, what="displ", view="inequivalent",
                  select_symbols=None, colormap="jet", xlims=None, ylims=None, sharey=False,
                  fontsize=10, verbose=0, **kwargs):
        """
        Plot phonon PJDOS for each atom in the unit cell.
        One subplot for each component, each subplot show all inequivalent sites.
        By default, only "inequivalent" atoms are shown.

        comparison of Ueq values, which
        are calculated as the mean of the diagonal elements of the harmonic ADP tensor, (d)
        comparison of the ADP anisotropy factor, which is defined as the ratio between maximum Uii
        and minimum Uii values. A ratio of 1 would correspond to an isotropic displacement.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate.
            components: "all" for all components. "diag" for diagonal elements, "offdiag" for off-diagonal terms only.
            what: "displ" for displament, "vel" for velocity.
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            select_symbols: String or list of strings with chemical symbols. Used to select only atoms of this type.
            colormap: matplotlib colormap.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            sharey: True if y-axis should be shared.
            fontsize: Legend and title fontsize.
            verbose: Verbosity level.

        Returns: |matplotlib-Figure|
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols=select_symbols, verbose=verbose)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=2, ncols=1,
                                                sharex=True, sharey=sharey, squeeze=True)
        cmap = plt.get_cmap(colormap)

        # Compute U(T)
        tmesh = np.linspace(tstart, tstop, num=num)
        msq = self.get_msq_tmesh(tmesh, iatom_list=aview.iatom_list, what_list=what)
        # [natom, 3, 3, nt]
        values = getattr(msq, what)
        ntemp = len(msq.tmesh)

        for ix, ax in enumerate(ax_list):
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            if what == "displ":
                ylabel = r"$U_{iso}\;(\AA^2)$" if ix == 0 else \
                         r"Anisotropy factor ($\dfrac{\epsilon_{max}}{\epsilon_{min}}}$)"
            elif what == "vel":
                ylabel = r"$V_{iso}\;(m/s)^2$" if ix == 0 else \
                         r"Anisotropy factor ($\dfrac{\epsilon_{max}}{\epsilon_{min}}}$)"
            else:
                raise ValueError("Unknown value for what: `%s`" % str(what))
            ax.set_ylabel(ylabel, fontsize=fontsize)

            # Plot this component for all inequivalent atoms on the same subplot.
            for ii, (iatom, site_label) in enumerate(zip(aview.iatom_list, aview.site_labels)):
                color = cmap(float(ii) / max((len(aview.iatom_list) - 1), 1))
                #msq.displ[iatom, 3, 3, nt]
                if ix == 0:
                    # ISO calculated as the mean of the diagonal elements of the harmonic ADP tensor
                    ys = np.trace(values[iatom]) / 3.0
                elif ix == 1:
                    # Ratio between maximum Uii and minimum Uii values.
                    # A ratio of 1 would correspond to an isotropic displacement.
                    ys = np.empty(ntemp)
                    for itemp in range(ntemp):
                        eigs = np.linalg.eigvalsh(values[iatom, :, :, itemp], UPLO='U')
                        ys[itemp] = eigs.max() / eigs.min()
                else:
                    raise ValueError("Invalid ix index: `%s" % ix)

                ax.plot(msq.tmesh, ys, label=site_label if ix == 0 else None,
                        color=color) #, marker="o")
                if ix == 0:
                    ax.legend(loc="best", fontsize=fontsize, shadow=True)

            if ix == len(ax_list) - 1:
                ax.set_xlabel("Temperature (K)")
            else:
                set_visible(ax, False, "xlabel", "xticklabels")

        return fig
