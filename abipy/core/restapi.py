"""
This module provides interfaces with the Materials Project REST API v2 to enable
the creation of data structures and pymatgen objects using Materials Project data.
"""
from __future__ import unicode_literals, print_function, division

import sys

from collections import OrderedDict
from pprint import pprint
from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from monty.string import marquee
from pymatgen import SETTINGS
try:
    from pymatgen.ext.matproj import MPRester, MPRestError
except ImportError:
    from pymatgen.matproj.rest import MPRester, MPRestError
from abipy.tools.printing import print_dataframe
from abipy.core.mixins import NotebookWriter


MP_DEFAULT_ENDPOINT = "https://www.materialsproject.org/rest/v2"

MP_KEYS_FOR_DATAFRAME = ("pretty_formula", "e_above_hull", "energy_per_atom",
                         "formation_energy_per_atom", "nsites", "volume",
                         "spacegroup.symbol", "spacegroup.number",
                         "band_gap", "total_magnetization", "material_id")
                         # "unit_cell_formula", "icsd_id", "icsd_ids", "cif", , "tags", "elasticity")


def get_mprester(api_key=None, endpoint=None):
    """
    Args:
        api_key (str): A String API key for accessing the MaterialsProject
            REST interface. Please apply on the Materials Project website for one.
            If this is None, the code will check if there is a `PMG_MAPI_KEY` in
            your .pmgrc.yaml. If so, it will use that environment
            This makes easier for heavy users to simply add
            this environment variable to their setups and MPRester can
            then be called without any arguments.
        endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
            Defaults to the standard Materials Project REST address, but
            can be changed to other urls implementing a similar interface.
    """
    if api_key is None:
        api_key = SETTINGS.get("PMG_MAPI_KEY")
        if api_key is None:
            raise RuntimeError("Cannot find PMG_MAPI_KEY in pymatgen settings. Add it to $HOME/.pmgrc.yaml")

    if endpoint is None: endpoint = MP_DEFAULT_ENDPOINT
    return MyMPRester(api_key=api_key, endpoint=endpoint)


class MyMPRester(MPRester):
    """
    Subclass Materials project Rester.
    See :cite:`Jain2013,Ong2015`.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: MyMPRester
    """
    Error = MPRestError

    def get_phasediagram_results(self, elements):
        """
        Contact the materials project database, fetch entries and build :class:``PhaseDiagramResults`` instance.

        Args:
            elements: List of chemical elements.
        """
        entries = self.get_entries_in_chemsys(elements, inc_structure="final")
        return PhaseDiagramResults(entries)


class PhaseDiagramResults(object):
    """
    Simplified interface to phase-diagram pymatgen API.

    Inspired to:

        https://anaconda.org/matsci/plotting-and-analyzing-a-phase-diagram-using-the-materials-api/notebook

    See also: :cite:`Ong2008,Ong2010`
    """
    def __init__(self, entries):
        self.entries = entries
        from abipy.core.structure import Structure
        for e in entries:
            e.structure = Structure.as_structure(e.structure)

        self.structures = [e.structure for e in entries]
        self.mpids = [e.entry_id for e in entries]

        # Create phase diagram.
        from pymatgen.analysis.phase_diagram import PhaseDiagram
        self.phasediagram = PhaseDiagram(self.entries)

    def plot(self, show_unstable=True, show=True):
        """
        Plot phase diagram.

        Args:
            show_unstable (float): Whether unstable phases will be plotted as
                well as red crosses. If a number > 0 is entered, all phases with
                ehull < show_unstable will be shown.
            show: True to show plot.

        Return:
            plotter object.
        """
        from pymatgen.analysis.phase_diagram import PDPlotter
        plotter = PDPlotter(self.phasediagram, show_unstable=show_unstable)
        if show:
            plotter.show()
        return plotter

    @lazy_property
    def dataframe(self):
        """Pandas dataframe with the most important results."""
        rows = []
        for e in self.entries:
            d = e.structure.get_dict4pandas(with_spglib=True)
            decomp, ehull = self.phasediagram.get_decomp_and_e_above_hull(e)

            rows.append(OrderedDict([
                ("Materials ID", e.entry_id),
                ("spglib_symb", d["spglib_symb"]), ("spglib_num", d["spglib_num"]),
                ("Composition", e.composition.reduced_formula),
                ("Ehull", ehull), # ("Equilibrium_reaction_energy", pda.get_equilibrium_reaction_energy(e)),
                ("Decomposition", " + ".join(["%.2f %s" % (v, k.composition.formula) for k, v in decomp.items()])),
            ]))

        import pandas as pd
        return pd.DataFrame(rows, columns=list(rows[0].keys()) if rows else None)

    def print_dataframes(self, with_spglib=False, file=sys.stdout, verbose=0):
        """
        Print pandas dataframe to file `file`.

        Args:
            with_spglib: True to compute spacegroup with spglib.
            file: Output stream.
            verbose: Verbosity level.
        """
        print_dataframe(self.dataframe, file=file)
        if verbose:
            from abipy.core.structure import dataframes_from_structures
            dfs = dataframes_from_structures(self.structures, index=self.mpids, with_spglib=with_spglib)
            print_dataframe(dfs.lattice, title="Lattice parameters:", file=file)
            if verbose > 1:
                print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):", file=file)


class DatabaseStructures(NotebookWriter):
    """
    Store the results of a query to the MP database.
    This object is immutable, use add_entry to create a new instance.
    """

    def __init__(self, structures, ids, data=None):
        """
        Args:
            structures: List of structure objects
            ids: List of database ids.
            data: List of dictionaries with data associated to the structures (optional).
	"""
        from abipy.core.structure import Structure
        self.structures = list(map(Structure.as_structure, structures))
        self.ids, self.data = ids, data
        assert len(self.structures) == len(ids)
        if data is not None:
            assert len(self.structures) == len(data)

    def __bool__(self):
        """bool(self)"""
        return bool(self.structures)
    __nonzero__ = __bool__  # py2

    def filter_by_spgnum(self, spgnum):
        """Filter structures by space group number. Return new MpStructures object."""
        inds = [i for i, s in enumerate(self.structures) if s.get_space_group_info()[1] == int(spgnum)]
        new_data = None if self.data is None else [self.data[i] for i in inds]
        return self.__class__([self.structures[i] for i in inds], [self.ids[i] for i in inds], data=new_data)

    def add_entry(self, structure, entry_id, data_dict=None):
        """
        Add new entry, return new object.

        Args:
           structure: New structure object.
           entry_id: ID associated to new structure.
           data_dict: Option dictionary with metadata.
        """
        if data_dict is None:
           new_data = None if self.data is None else self.data + [{}]
        else:
            assert self.data is not None
            new_data = self.data + [data_dict]

        return self.__class__(self.structures + [structure], self.ids + [entry_id], data=new_data)

    @property
    def lattice_dataframe(self):
        """pandas DataFrame with lattice parameters."""
        return self.structure_dataframes.lattice

    @property
    def coords_dataframe(self):
        """pandas DataFrame with atomic positions."""
        return self.structure_dataframes.coords

    @lazy_property
    def structure_dataframes(self):
        """Pandas dataframes constructed from self.structures."""
        from abipy.core.structure import dataframes_from_structures
        return dataframes_from_structures(self.structures, index=self.ids, with_spglib=True)

    def print_results(self, fmt="abivars", verbose=0, file=sys.stdout):
        """
        Print pandas dataframe, structures using format `fmt`, and data to file `file`.
        `fmt` is automaticall set to `cif` if structure is disordered.
        Set fmt to None or empty string to disable structure output.
        """
        print("\n# Found %s structures in %s database (use `verbose` to get further info)\n"
                % (len(self.structures), self.dbname), file=file)

        if self.dataframe is not None: print_dataframe(self.dataframe, file=file)
        if verbose and self.data is not None: pprint(self.data, stream=file)

        # Print structures
        print_structures = not (fmt is None or str(fmt) == "None")
        if print_structures:
            for i, structure in enumerate(self.structures):
                print(" ", file=file)
                print(marquee("%s input for %s" % (fmt, self.ids[i]), mark="#"), file=file)
                print("# " + structure.spget_summary(verbose=verbose).replace("\n", "\n# ") + "\n", file=file)
                if not structure.is_ordered:
                    print(structure.convert(fmt="cif"), file=file)
                else:
                    print(structure.convert(fmt=fmt), file=file)
                print(2 * "\n", file=file)

        if len(self.structures) > 10:
            # Print info again
            print("\n# Found %s structures in %s database (use `verbose` to get further info)\n"
                    % (len(self.structures), self.dbname), file=file)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """NOP required by NotebookWriter protocol."""
        yield None

    def write_notebook(self, nbpath=None, title=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("dbs = abilab.restapi.DatabaseStructures.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("import qgrid"),
            nbv.new_code_cell("# dbs.print_results(fmt='cif', verbose=0)"),
            nbv.new_code_cell("# qgrid.show_grid(dbs.lattice_dataframe)"),
            nbv.new_code_cell("# qgrid.show_grid(dbs.coords_dataframe)"),
            nbv.new_code_cell("qgrid.show_grid(dbs.dataframe)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class MpStructures(DatabaseStructures):
    """
    Store the results of a query to the Materials Project database.

    .. inheritance-diagram:: MpStructures
    """
    dbname = "Materials Project"

    @lazy_property
    def dataframe(self):
        """Pandas dataframe constructed from self.data. None if data is not available."""
        if not self.data: return None
        import pandas as pd
        rows = []
        for d, structure in zip(self.data, self.structures):
            d = Dotdict(d)
            d = OrderedDict([(k, d.dotget(k, default=None)) for k in MP_KEYS_FOR_DATAFRAME])
            # Add lattice parameters.
            for k in ("a", "b", "c", "alpha", "beta", "gamma"):
                d[k] = getattr(structure.lattice, k)
            rows.append(d)

        return pd.DataFrame(rows, index=self.ids, columns=list(rows[0].keys()))

    def open_browser(self, browser=None, limit=10):
        """
        Args:
            browser: Open webpage in ``browser``. Use default if $BROWSER if None.
            limit: Max number of tabs opened in browser. None for no limit.
        """
        import webbrowser
        import cgi
        for i, mpid in enumerate(self.ids):
            if limit is not None and i >= limit:
                print("Found %d structures found. Won't open more than %d tabs" % (len(self.ids), limit))
                break
            # https://materialsproject.org/materials/mp-2172/
            url = "https://materialsproject.org/materials/%s/" % mpid
            webbrowser.get(browser).open_new_tab(cgi.escape(url))


class CodStructures(DatabaseStructures):
    """
    Store the results of a query to the COD_ database :cite:`Grazulis2011`.

    .. inheritance-diagram:: CodStructures
    """
    dbname = "COD"

    @lazy_property
    def dataframe(self):
        """
        |pandas-Dataframe| constructed. Essentially geometrical info and space groups found by spglib_
        as COD API is rather limited.
        """
        df = self.lattice_dataframe.copy()
        # Add space group from COD
        df["cod_sg"] = [d.get("sg", "").replace(" ", "") for d in self.data]
        return df


class Dotdict(dict):

    def dotget(self, key, default=None):
        """
        d.dotget["foo.bar"] --> d["foo"]["bar"] if "foo.bar" not in self
        """
        # if key is in dict access as normal
        if key in self:
            return super(Dotdict,self).__getitem__(key)

        # Assume string
        i = -1
        try:
            i = key.find(".")
            if i == -1: return default
        except AttributeError:
            return default

        try:
            root, key = key[:i], key[i+1:]
            if key == ".": return None
            return Dotdict(**self[root])[key]
        except Exception:
            return None