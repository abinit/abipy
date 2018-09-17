# coding: utf-8
"""
This modules provides subclasses of pymatgen tensor objects.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd

try:
    from pymatgen.core.tensors import Tensor, SquareTensor
except ImportError:
    # Can be removed in v2019.1.1
    from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.stress import Stress as pmg_Stress
from pymatgen.analysis.piezo import PiezoTensor
from abipy.iotools import ETSF_Reader


class _Tensor33(object):

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.get_dataframe()._repr_html_()

    def get_dataframe(self, tol=1e-3):
        """Return |pandas-Dataframe| with tensor elements set to zero below `tol`."""
        tensor = self.zeroed(tol=tol)
        return pd.DataFrame({"x": tensor[:,0], "y": tensor[:,1], "z": tensor[:,2]}, index=["x", "y", "z"])

    def get_voigt_dataframe(self, tol=1e-3):
        """
        Return |pandas-DataFrame| with Voigt indices as colums.
        Elements below tol are set to zero.

        Useful to analyze the converge of individual elements.
        """
        tensor = self.zeroed(tol=tol)
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        d = {k: v for k, v in zip(columns, tensor.voigt)}
        return pd.DataFrame(d, index=[0], columns=columns)


class Stress(pmg_Stress, _Tensor33):
    """
    Stress tensor. rank2 symmetric tensor with shape [3, 3].

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Stress
    """


class DielectricTensor(SquareTensor, _Tensor33):
    """
    Subclass of |pmg-Tensor| describing a dielectric tensor.
    rank2 symmetric tensor with shape [3, 3].

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DielectricTensor
    """

    def reflectivity(self, n1=1, tol=1e-6):
        """
        If the tensor is diagonal (with off diagonal elements smaller than tol)
        returns the three components of the reflectivity

            :math:`|n1 - n2| / | n1 + n2 |`
        """
        d = np.diag(self)

        if np.max(np.abs(self - np.diag(d))) > tol:
            raise ValueError("The tensor is not diagonal.")

        n2 = np.sqrt(d)

        return np.abs((n1 - n2) / (n1 + n2)) ** 2


class ZstarTensor(SquareTensor, _Tensor33):
    """
    Born effective charge tensor (for a single atom).

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZstarTensor
    """


class NLOpticalSusceptibilityTensor(Tensor):
    """
    Subclass of |pmg-Tensor| containing the non-linear optical susceptibility tensor.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: NLOpticalSusceptibilityTensor
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Creates the tensor from an anaddb.nc netcdf file containing ``dchide``.
        This requires to run anaddb with ``tnlflag`` > 0
        """
        with ETSF_Reader(filepath) as reader:
            try:
                return cls(reader.read_value("dchide"))
            except Exception as exc:
                import traceback
                msg = traceback.format_exc()
                msg += ("Error while trying to read from file.\n"
                        "Verify that nlflag > 0 in anaddb\n")
                raise ValueError(msg)
