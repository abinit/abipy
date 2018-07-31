# coding: utf-8
"""
Objects to analyze elastic and piezoelectric tensors computed by anaddb.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import pandas as pd

from collections import OrderedDict
from monty.string import is_string, list_strings, marquee
from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.analysis.elasticity.tensors import Tensor
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.piezo import PiezoTensor
from abipy.core.mixins import Has_Structure
from abipy.flowtk.netcdf import ETSF_Reader


class MyElasticTensor(ElasticTensor):

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.to_voigt_dataframe()._repr_html_()

    def to_voigt_dataframe(self, tol=1e-5):
        tensor = self.zeroed(tol) if tol else self
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        #columns = ["1", "2", "3", "4", "5", "6"]
        rows = []
        for row in tensor.voigt:
            rows.append({k: v for k, v in zip(columns, row)})

        return pd.DataFrame(rows, index=columns, columns=columns)


class MyPiezoTensor(PiezoTensor):

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.to_voigt_dataframe()._repr_html_()

    def to_voigt_dataframe(self, tol=1e-5):
        tensor = self.zeroed(tol) if tol else self
        index = ["Px", "Py", "Pz"]
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        #index = ["P1", "P2", "P3"]
        #columns = ["1", "2", "3", "4", "5", "6"]
        rows = []
        for irow, row in enumerate(tensor.voigt):
            rows.append({k: v for k, v in zip(columns, row)})

        return pd.DataFrame(rows, index=index, columns=columns)


class ElasticData(Has_Structure):
    """
    Container with the different elastic and piezoelectric tensors
    computed by anaddb. Data is tored in pymatgen tensor objects.

    Provides methods to analyze/tabule data
    http://progs.coudert.name/elate/mp?query=mp-2172
    """

    ALL_ELASTIC_TENSOR_NAMES = (
        "elastic_clamped", "elastic_relaxed",
        "elastic_stress_corr", "elastic_relaxed_fixed_D",
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

    TENSOR_META = {
        "elastic_clamped": AttrDict(
            info="clamped-ion elastic tensor in Voigt notation (shape: (6, 6))",
            units="GPa"),
        "elastic_relaxed": AttrDict(
            info="relaxed-ion elastic tensor in Voigt notation (shape: (6, 6))",
            units="GPa"),
        "elastic_stress_corr": AttrDict(
            info="relaxed-ion elastic tensor considering the stress left inside cell in Voigt notation (shape: (6, 6))",
            units="GPa"),
        "elastic_relaxed_fixed_D": AttrDict(
            info="relaxed-ion elastic tensor at fixed displacement field in Voigt notation (shape: (6, 6))",
            units="GPa"),
        "piezo_clamped": AttrDict(
            info="clamped-ion piezoelectric tensor in Voigt notation (shape: (3, 6))",
            units="c/m^2"),
        "piezo_relaxed": AttrDict(
            info="relaxed-ion piezoelectric tensor in Voigt notation (shape: (3, 6))",
            units="c/m^2"),
        "d_piezo_relaxed": AttrDict(
            info="relaxed-ion piezoelectric d tensor in Voigt notation (shape: (3, 6))",
            units="pc/m^2"),
        "g_piezo_relaxed": AttrDict(
            info="relaxed-ion piezoelectric g tensor in Voigt notation (shape: (3, 6))",
            units="m^2/c"),
        "h_piezo_relaxed": AttrDict(
            info="relaxed-ion piezoelectric h tensor in Voigt notation (shape: (3, 6))",
            units="GN/c"),
    }

    def __init__(self, structure, params, elastic_clamped=None, elastic_relaxed=None, elastic_stress_corr=None,
                 elastic_relaxed_fixed_D=None, piezo_clamped=None, piezo_relaxed=None, d_piezo_relaxed=None,
                 g_piezo_relaxed=None, h_piezo_relaxed=None):
        """
        Args:
            structure: |Structure| object.
            params: Dictionary with input parameters.
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
        self.params = params
        self.elastic_clamped = self._define_variable(elastic_clamped, MyElasticTensor)
        self.elastic_relaxed = self._define_variable(elastic_relaxed, MyElasticTensor)
        self.elastic_stress_corr = self._define_variable(elastic_stress_corr, MyElasticTensor)
        self.elastic_relaxed_fixed_D = self._define_variable(elastic_relaxed_fixed_D, MyElasticTensor)
        self.piezo_clamped = self._define_variable(piezo_clamped, MyPiezoTensor)
        self.piezo_relaxed = self._define_variable(piezo_relaxed, MyPiezoTensor)
        self.d_piezo_relaxed = self._define_variable(d_piezo_relaxed, MyPiezoTensor)
        self.g_piezo_relaxed = self._define_variable(g_piezo_relaxed, MyPiezoTensor)
        self.h_piezo_relaxed = self._define_variable(h_piezo_relaxed, MyPiezoTensor)

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

        params = AttrDict(
            # NB: asr and chneut are always present in the new anaddb.nc file
            # Use -666 to support old formats.
            asr=int(reader.read_value("asr", default=-666)),
            chneut= int(reader.read_value("chneut", default=-666)),
            elaflag=int(reader.read_value("elaflag", default=0)),
            instrflag=int(reader.read_value("instrflag", default=0)),
            piezoflag=int(reader.read_value("piezoflag", default=0)),
            dieflag=int(reader.read_value("dieflag", default=0)),
        )

        ts = AttrDict({n: None for n in cls.ALL_TENSOR_NAMES})

        # [6, 6] symmetric tensors (written by Fortran, produced in ddb_elast)
        ts.elastic_clamped = reader.read_value("elastic_constants_clamped_ion", default=None)
        ts.elastic_relaxed = reader.read_value("elastic_constants_relaxed_ion", default=None)
        if params.elaflag == 5:
            ts.elastic_stress_corr = reader.read_value("elastic_constants_relaxed_ion_stress_corrected")

        # Written in ddb_piezo
        ts.elastic_relaxed_fixed_D = reader.read_value("elastic_tensor_relaxed_ion_fixed_D", default=None)

        # [3, 6] tensors (written by Fortran, produced in ddb_piezo).
        ts.piezo_clamped = reader.read_value("piezo_clamped_ion", default=None)
        ts.piezo_relaxed = reader.read_value("piezo_relaxed_ion", default=None)
        ts.d_piezo_relaxed = reader.read_value("d_tensor_relaxed_ion", default=None)

        # These are [6, 3] tensors written by Fortran (need to transpose).
        ts.g_piezo_relaxed = reader.read_value("g_tensor_relaxed_ion", default=None)
        if ts.g_piezo_relaxed is not None:
            ts.g_piezo_relaxed = ts.g_piezo_relaxed.T.copy()
        ts.h_piezo_relaxed = reader.read_value("h_tensor_relaxed_ion", default=None)
        if ts.h_piezo_relaxed is not None:
            ts.h_piezo_relaxed = ts.h_piezo_relaxed.T.copy()

        return cls(structure, params, **ts)

        #return cls(structure, params, elastic_clamped=elastic_clamped, elastic_relaxed=elastic_relaxed,
        #           elastic_stress_corr=elastic_stress_corr, elastic_relaxed_fixed_D=elastic_relaxed_fixed_D,
        #           piezo_clamped=piezo_clamped, piezo_relaxed=piezo_relaxed, d_piezo_relaxed=d_piezo_relaxed,
        #           g_piezo_relaxed=g_piezo_relaxed, h_piezo_relaxed=h_piezo_relaxed)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String represention with verbosity level `verbose`."""
        lines = []; app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(marquee("Parameters", mark="="))
        import json
        app(json.dumps(self.params, indent=2, sort_keys=True))

        for tensor_type in ("elastic", "piezoelectric"):
            name_tensor_list = self.name_tensor_list(tensor_type=tensor_type)
            if name_tensor_list:
                app("")
                app(marquee("%s tensors available" % tensor_type, mark="="))
                for name, tensor in name_tensor_list:
                    meta = self.TENSOR_META[name]
                    is_fit = tensor.is_fit_to_structure(self.structure, tol=1e-2)
                    # TODO
                    tol = dict(elastic=1e-3, piezoelectric=1e-5)[tensor_type]
                    app("[%s]" % name.upper())
                    app("%s" % meta.info)
                    app("Units: %s, set to zero below: %s, fit_to_structure: %s" % (meta.units, tol, is_fit))
                    app("")
                    if tensor_type == "elastic":
                        app(self.get_elastic_tensor_dataframe(tensor_name=name, tol=tol).to_string())
                    elif tensor_type == "piezoelectric":
                        app(self.get_piezoelectric_tensor_dataframe(tensor_name=name, tol=tol).to_string())
                    app("")

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
                if name not in self.TYPE2NAMES[tensor_type]:
                    raise ValueError("tensor name %s does not belong to type: `%s`" % (name, tensor_type))
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

        return self.__class__(structure, self.params, **kwargs)

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

        return self.__class__(structure, self.params, **kwargs)

    def get_elastic_tensor_dataframe(self, tensor_name="elastic_relaxed", tol=1e-3):
        """
        Args:
            tol: returns the matrix with all entries below a certain threshold
            (i.e. tol) set to zero
        """
        tensor = getattr(self, tensor_name)
        if tensor is None: return pd.DataFrame()
        if tol: tensor = tensor.zeroed(tol)
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        #columns = ["1", "2", "3", "4", "5", "6"]
        rows = []
        for row in tensor.voigt:
            rows.append({k: v for k, v in zip(columns, row)})

        return pd.DataFrame(rows, index=columns, columns=columns)

    def get_piezoelectric_tensor_dataframe(self, tensor_name="piezo_relaxed", tol=1e-5):
        """
        Args:
            tol: returns the matrix with all entries below a certain threshold
            (i.e. tol) set to zero
        """
        tensor = getattr(self, tensor_name)
        if tensor is None: return pd.DataFrame()
        if tol: tensor = tensor.zeroed(tol)
        index = ["Px", "Py", "Pz"]
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        #index = ["P1", "P2", "P3"]
        #columns = ["1", "2", "3", "4", "5", "6"]
        rows = []
        for row in tensor.voigt:
            rows.append({k: v for k, v in zip(columns, row)})

        return pd.DataFrame(rows, index=index, columns=columns)

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
                # Add column telling whether fit has been performed
                if len(do_fits) > 1: d["fit_to_structure"] = do_fit
                d["tensor_name"] = name
                rows.append(d)

        return pd.DataFrame(rows, columns=list(rows[0].keys() if rows else None))
