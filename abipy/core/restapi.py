"""
This module provides interfaces with the Materials Project REST API v2 to enable
the creation of data structures and pymatgen objects using Materials Project data.
"""
from __future__ import division, unicode_literals, print_function, division

from collections import OrderedDict
from monty.functools import lazy_property
from pymatgen import SETTINGS
from pymatgen.matproj.rest import MPRester, MPRestError


MP_DEFAULT_ENDPOINT = "https://www.materialsproject.org/rest/v2"

MP_KEYS_FOR_DATAFRAME = ("pretty_formula", "e_above_hull", "energy_per_atom",
                         "formation_energy_per_atom", "nsites", "volume",
                         "spacegroup.symbol", "spacegroup.number",
                         "band_gap", "total_magnetization", "material_id")
                         # "unit_cell_formula", "icsd_id", "icsd_ids", "cif", , "tags", "elasticity")


class MyMPRester(MPRester):

    Error = MPRestError


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


class MpStructures(object):
    """Store the results of a query to the MP database."""

    def __init__(self, structures, mpids, data=None):
	"""
	Args:
	    structures: List of structure objects
	    mpids: List of MaterialsProject ids.
	    data: List of dictionaries with data associated to the structures (optional).
	"""
        self.structures, self.mpids = structures, mpids
        self.data = data

    @lazy_property
    def table(self):
	"""Pandas dataframe constructed from self.data. None if data is not available."""
        if self.data is None: return None
        import pandas as pd
        rows = []
        for d in self.data:
            d = Dotdict(d)
            rows.append(OrderedDict([(k, d.dotget(k, default=None)) for k in MP_KEYS_FOR_DATAFRAME]))
            table = pd.DataFrame(rows, index=[r["material_id"] for r in rows],
                                 columns=list(rows[0].keys()))

    #def __str__(self):
    #    lines = []
    #    return "\n".join(lines)


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
