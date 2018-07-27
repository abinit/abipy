# coding: utf-8
"""
AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from monty.functools import lazy_property
from abipy.core.mixins import Has_Structure
from abipy.flowtk.netcdf import ETSF_Reader
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.tensors import Tensor
from pymatgen.analysis.piezo import PiezoTensor


class ElasticData(Has_Structure):
    """
    Handles all kind of elastic data return from abinit/anaddb relying on pymatgen tensor objects.
    """

    def __init__(self, structure, elastic_clamped=None, elastic_relaxed=None, elastic_stress_corr=None,
                 elastic_relaxed_fixed=None, piezo_clamped=None, piezo_relaxed=None, d_piezo_relaxed=None,
                 g_piezo_relaxed=None, h_piezo_relaxed=None):
        """
        Args:
            structure: the structure
            elastic_clamped: the values of a clamped-ion elastic tensor in Voigt notation (shape (6,3)) in GPa.
            elastic_relaxed: the values of a relaxed-ion elastic tensor in Voigt notation (shape (6,3)) in GPa.
            elastic_stress_corr: the values of a relaxed-ion elastic tensor considering the stress left inside cell
                in Voigt notation (shape (6,3)) in GPa.
            elastic_relaxed_fixed: the values of a relaxed-ion elastic tensor at fixed displacement field
                in Voigt notation (shape (6,3)) in GPa.
            piezo_clamped: the values of a clamped-ion piezoelectric tensor in Voigt notation (shape (6,3)) in c/m^2.
            piezo_relaxed: the values of a relaxed-ion piezoelectric tensor in Voigt notation (shape (6,3)) in c/m^2.
            d_piezo_relaxed: the values of a relaxed-ion piezoelectric d tensor in Voigt notation (shape (6,3)) in pc/m^2.
            g_piezo_relaxed: the values of a relaxed-ion piezoelectric g tensor in Voigt notation (shape (6,3)) in m^2/c.
            h_piezo_relaxed: the values of a relaxed-ion piezoelectric h tensor in Voigt notation (shape (6,3)) in GN/c.
        """
        self._structure = structure
        self.elastic_clamped = self._define_variable(elastic_clamped, ElasticTensor)
        self.elastic_relaxed = self._define_variable(elastic_relaxed, ElasticTensor)
        self.elastic_stress_corr = self._define_variable(elastic_stress_corr, ElasticTensor)
        self.elastic_relaxed_fixed = self._define_variable(elastic_relaxed_fixed, ElasticTensor)
        self.piezo_clamped = self._define_variable(piezo_clamped, PiezoTensor)
        self.piezo_relaxed = self._define_variable(piezo_relaxed, PiezoTensor)
        self.d_piezo_relaxed = self._define_variable(d_piezo_relaxed, PiezoTensor)
        self.g_piezo_relaxed = self._define_variable(g_piezo_relaxed, Tensor)
        self.h_piezo_relaxed = self._define_variable(h_piezo_relaxed, Tensor)

    def _define_variable(self, tensor_voigt, tensor_class):
        """
        Helper function to set values of a variable
        """
        if not tensor_voigt:
            return None
        else:
            return tensor_class.from_voigt(tensor_voigt)

    @property
    def structure(self):
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

        elastic_clamped = reader.read_value("elastic_constants_clamped_ion", default=None)
        elastic_relaxed = reader.read_value("elastic_constants_relaxed_ion", default=None)
        elastic_stress_corr = reader.read_value("elastic_constants_relaxed_ion_stress_corrected",default=None)
        elastic_relaxed_fixed = reader.read_value("elastic_tensor_relaxed_ion_fixed_D", default=None)

        piezo_clamped = reader.read_value("piezo_clamped_ion", default=None)
        if piezo_clamped is not None:
            piezo_clamped = piezo_clamped.T.copy()
        piezo_relaxed = reader.read_value("piezo_relaxed_ion", default=None)
        if piezo_relaxed is not None:
            piezo_relaxed = piezo_relaxed.T.copy()
        d_piezo_relaxed = reader.read_value("d_tensor_relaxed_ion", default=None)
        if d_piezo_relaxed is not None:
            d_piezo_relaxed = d_piezo_relaxed.T.copy()
        g_piezo_relaxed = reader.read_value("g_tensor_relaxed_ion", default=None)
        h_piezo_relaxed = reader.read_value("h_tensor_relaxed_ion", default=None)

        return cls(structure=structure, elastic_clamped=elastic_clamped, elastic_relaxed=elastic_relaxed,
                   elastic_stress_corr=elastic_stress_corr, elastic_relaxed_fixed=elastic_relaxed_fixed,
                   piezo_clamped=piezo_clamped, piezo_relaxed=piezo_relaxed, d_piezo_relaxed=d_piezo_relaxed,
                   g_piezo_relaxed=g_piezo_relaxed, h_piezo_relaxed=h_piezo_relaxed)
