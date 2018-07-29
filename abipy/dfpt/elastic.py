# coding: utf-8
"""
Objects to analyze elastic and piezoelectric tensors computed by anaddb.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import pandas as pd

from collections import OrderedDict
from monty.string import is_string, list_strings, marquee
from monty.functools import lazy_property
from pymatgen.analysis.elasticity.tensors import Tensor
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.piezo import PiezoTensor
from abipy.core.mixins import Has_Structure
from abipy.flowtk.netcdf import ETSF_Reader


class ElasticData(Has_Structure):
    """
    Container with the different elastic and piezoelectric tensors
    computed by anaddb. Data is tored in pymatgen tensor objects.

    Provides methods to analyze/tabule data
    http://progs.coudert.name/elate/mp?query=mp-2172
    """

    ALL_ELASTIC_TENSOR_NAMES = (
        "elastic_clamped", "elastic_relaxed", "elastic_stress_corr", "elastic_relaxed_fixed_D",
    )

    ALL_PIEZOELECTRIC_TENSOR_NAMES = (
        "piezo_clamped", "piezo_relaxed",
        "d_piezo_relaxed", "g_piezo_relaxed", "h_piezo_relaxed"
    )

    ALL_TENSOR_NAMES = ALL_ELASTIC_TENSOR_NAMES + ALL_PIEZOELECTRIC_TENSOR_NAMES

    TYPE2NAMES = {
        "all": ALL_TENSOR_NAMES,
        "elastic": ALL_ELASTIC_TENSOR_NAMES,
        "piezoelectric": ALL_PIEZOELECTRIC_TENSOR_NAMES,
    }

    def __init__(self, structure, elastic_clamped=None, elastic_relaxed=None, elastic_stress_corr=None,
                 elastic_relaxed_fixed_D=None, piezo_clamped=None, piezo_relaxed=None, d_piezo_relaxed=None,
                 g_piezo_relaxed=None, h_piezo_relaxed=None):
        """
        Args:
            structure: |Structure| object.
            elastic_clamped: clamped-ion elastic tensor in Voigt notation (shape (6,6)) in GPa.
            elastic_relaxed: relaxed-ion elastic tensor in Voigt notation (shape (6,6)) in GPa.
            elastic_stress_corr: relaxed-ion elastic tensor considering the stress left inside cell
                in Voigt notation (shape (6,6)) in GPa.
            elastic_relaxed_fixed_D: relaxed-ion elastic tensor at fixed displacement field
                in Voigt notation (shape (6,6)) in GPa.
            piezo_clamped: clamped-ion piezoelectric tensor in Voigt notation (shape (3,6)) in c/m^2.
            piezo_relaxed: relaxed-ion piezoelectric tensor in Voigt notation (shape (3,6)) in c/m^2.
            d_piezo_relaxed: relaxed-ion piezoelectric d tensor in Voigt notation (shape (3,6)) in pc/m^2.
            g_piezo_relaxed: relaxed-ion piezoelectric g tensor in Voigt notation (shape (3,6)) in m^2/c.
            h_piezo_relaxed: relaxed-ion piezoelectric h tensor in Voigt notation (shape (3,6)) in GN/c.

        .. note::

            Arguments can be either arrays or Tensor objects.
        """
        self._structure = structure
        self.elastic_clamped = self._define_variable(elastic_clamped, ElasticTensor)
        self.elastic_relaxed = self._define_variable(elastic_relaxed, ElasticTensor)
        self.elastic_stress_corr = self._define_variable(elastic_stress_corr, ElasticTensor)
        self.elastic_relaxed_fixed_D = self._define_variable(elastic_relaxed_fixed_D, ElasticTensor)
        self.piezo_clamped = self._define_variable(piezo_clamped, PiezoTensor)
        self.piezo_relaxed = self._define_variable(piezo_relaxed, PiezoTensor)
        self.d_piezo_relaxed = self._define_variable(d_piezo_relaxed, PiezoTensor)
        self.g_piezo_relaxed = self._define_variable(g_piezo_relaxed, Tensor)
        self.h_piezo_relaxed = self._define_variable(h_piezo_relaxed, Tensor)

    def _define_variable(self, tensor_voigt, tensor_class):
        """
        Helper function to set values of a variable
        """
        if isinstance(tensor_voigt, Tensor):
            if not isinstance(tensor_voigt, tensor_class):
                raise TypeError("Expecting tensor class `%s`, received class `%s`" % (
                    tensor_class.__name__, tensor_voigt.__class__.__name__))
            return tensor_voigt
        else:
            return tensor_class.from_voigt(tensor_voigt) if tensor_voigt is not None else None

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    @classmethod
    def from_file(cls, path):
        """
        Builds the object from an anaddb.nc file
        """
        with ETSF_Reader(path) as reader:
            return cls.from_ncreader(reader)

    @classmethod
    def from_ncreader(cls, reader):
        """
        Builds the object from a ETSF_Reader
        """

        structure = reader.read_structure()
        # [6, 6] symmetric tensors (written by Fortran)
        # Produced in ddb_elast
        elastic_clamped = reader.read_value("elastic_constants_clamped_ion", default=None)
        elastic_relaxed = reader.read_value("elastic_constants_relaxed_ion", default=None)
        elastic_stress_corr = reader.read_value("elastic_constants_relaxed_ion_stress_corrected",default=None)
        # ddb_piezo
        elastic_relaxed_fixed_D = reader.read_value("elastic_tensor_relaxed_ion_fixed_D", default=None)

        # [3, 6] tensors
        # Produced in ddb_piezo
        piezo_clamped = reader.read_value("piezo_clamped_ion", default=None)
        piezo_relaxed = reader.read_value("piezo_relaxed_ion", default=None)
        d_piezo_relaxed = reader.read_value("d_tensor_relaxed_ion", default=None)

        # These are  [6, 3] tensors written by Fortran (need to transpose)
        g_piezo_relaxed = reader.read_value("g_tensor_relaxed_ion", default=None)
        if g_piezo_relaxed is not None:
            g_piezo_relaxed = g_piezo_relaxed.T.copy()
        h_piezo_relaxed = reader.read_value("h_tensor_relaxed_ion", default=None)
        if h_piezo_relaxed is not None:
            h_piezo_relaxed = h_piezo_relaxed.T.copy()

        return cls(structure=structure, elastic_clamped=elastic_clamped, elastic_relaxed=elastic_relaxed,
                   elastic_stress_corr=elastic_stress_corr, elastic_relaxed_fixed_D=elastic_relaxed_fixed_D,
                   piezo_clamped=piezo_clamped, piezo_relaxed=piezo_relaxed, d_piezo_relaxed=d_piezo_relaxed,
                   g_piezo_relaxed=g_piezo_relaxed, h_piezo_relaxed=h_piezo_relaxed)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representing with verbosity level `verbose`."""
        lines = []; app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        for tensor_type in ["elastic", "piezoelectric"]:
            name_tensor_list = self.name_tensor_list(tensor_type=tensor_type)
            if name_tensor_list:
                app("")
                app(marquee("Available %s tensors" % tensor_type, mark="="))
                for name, tensor in name_tensor_list:
                    ans = tensor.is_fit_to_structure(self.structure, tol=1e-2)
                    app("%s (fit_to_structure: %s)" % (name, ans))

        return "\n".join(lines)

    def name_tensor_list(self, tensor_names=None, tensor_type="all"):
        """
        List of (name, tensor) tuples. Only tensors stored in the object are returned.

        Args:
            tensor_names: List of tensor names to select. None means all.
            tensor_type: Select tensors by type. Must be in ["all", "elastic", "piezoelectric"].
        """
        l = []
        if tensor_names is None:
            for name in self.TYPE2NAMES[tensor_type]:
                tensor = getattr(self, name)
                if tensor is not None: l.append((name, tensor))
        else:
            for name in list_strings(tensor_names):
                tensor = getattr(self, name)
                if tensor is not None: l.append((name, tensor))

        return l

    def fit_to_structure(self, structure=None, symprec=0.1):
        """
        Return new ElasticData object with tensors that are invariant with respect to symmetry
        operations corresponding to `structure`.

        Args:
            structure (Structure): structure from which to generate symmetry operations
                If None, the internal structure is used.
            symprec (float): symmetry tolerance for the Spacegroup Analyzer
                used to generate the symmetry operations
        """
        structure = self.structure if structure is None else structure
        kwargs = {name: tensor.fit_to_structure(structure, symprec=symprec)
            for name, tensor in self.name_tensor_list()}

        return self.__class__(structure, **kwargs)

    def convert_to_ieee(self, structure=None, initial_fit=True, refine_rotation=True):
        """
        Return new set of tensors in IEEE format according to the 1987 IEEE standards.

        Args:
            structure (Structure): a structure associated with the
                tensor to be converted to the IEEE standard
                If None, the internal structure is used
            initial_fit (bool): flag to indicate whether initial
                tensor is fit to the symmetry of the structure.
                Defaults to true. Note that if false, inconsistent
                results may be obtained due to symmetrically
                equivalent, but distinct transformations
                being used in different versions of spglib.
            refine_rotation (bool): whether to refine the rotation
                produced by the ieee transform generator, default True
        """
        kwargs = {}
        structure = self.structure if structure is None else structure
        for name, tensor in self.name_tensor_list():
            # TODO: one should use the ieee stucture.
            kwargs[name] = tensor.convert_to_ieee(structure,
                initial_fit=initial_fit, refine_rotation=refine_rotation)

        return self.__class__(structure, **kwargs)

    def get_voigt_dataframe(self, tensor_names):
        """
        Return a |pandas-DataFrame| with voigt indices as colums.
        Useful to analyze the converge of individual elements of the tensor(s)

        Args:
            tensor_names: List of tensor names.
        """
        rows = []
        for name, tensor in self.name_tensor_list(tensor_names=tensor_names):
            voigt_map = tensor.get_voigt_dict(tensor.rank)
            row = {}
            for ind in voigt_map:
                row[voigt_map[ind]] = tensor[ind]
            row = OrderedDict(sorted(row.items(), key=lambda item: item[0]))
            row["tensor_name"] = name
            rows.append(row)


        return pd.DataFrame(rows, columns=list(rows[0].keys() if rows else None))

    #def get_average_elastic_dataframe(self, tensor_names="elastic_relaxed", fit_to_structure=False, symprec=0.1):
    #    """
    #    Return a |pandas-DataFrame|

    #    Args:
    #        tensor_names
    #    """
    #    schemes = ["voigt", "reuss", "vrh"]
    #    what_list = ["k", "g"] # "E", "v"]  # TODO
    #    rows = []
    #    for name, tensor in self.name_tensor_list(tensor_names=tensor_names):
    #        if fit_to_structure:
    #            tensor = tensor.fit_to_structure(self.structure, symprec=symprec)
    #        for scheme in schemes:
    #            row = OrderedDict()
    #            row["scheme"] = scheme
    #            row["tensor_name"] = name
    #            for what in what_list:
    #                # k_vrh
    #                aname = what + "_" + scheme
    #                row[what] = getattr(tensor, aname)
    #            rows.append(row)

    #    return pd.DataFrame(rows, columns=list(rows[0].keys()) if rows else None)

    def get_elast_properties_dataframe(self, tensor_names=["elastic_relaxed", "elastic_clamped"],
            include_base_props=True, ignore_errors=False, fit_to_structure=False, symprec=0.1):
        """
        Return a |pandas-DataFrame| with properties derived from the elastic tensor
        and the associated structure

        Args:
            tensor_names= ["elastic_relaxed", "elastic_clamped", "elastic_stress_corr", "elastic_relaxed_fixed_D"]
            include_base_props (bool): whether to include base properties, like k_vrh, etc.
            ignore_errors (bool): if set to true, will set problem properties
                that depend on a physical tensor to None, defaults to False
            fit_to_structure (bool): If True, properties are computed with the orginal tensors
                and symmetrized tensors. An additional column `fit_to_structure` is added to the dataframe.
            symprec (float): symmetry tolerance for the Spacegroup Analyzer
                used to generate the symmetry operations if `fit_to_structure`
        """
        do_fits = [False] if not fit_to_structure else [True, False]
        rows = []
        for name, tensor in self.name_tensor_list(tensor_names=tensor_names):
            for do_fit in do_fits:
                if do_fit:
                    tensor = tensor.fit_to_structure(self.structure, symprec=symprec)
                d = tensor.get_structure_property_dict(self.structure,
                        include_base_props=include_base_props, ignore_errors=ignore_errors)
                d.pop("structure")
                if len(do_fits) > 1:
                    # Add column telling whether fit has been performed
                    d["fit_to_structure"] = do_fit
                d["tensor_name"] = name
                rows.append(d)

        return pd.DataFrame(rows, columns=list(rows[0].keys() if rows else None))
