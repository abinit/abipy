# coding: utf-8
"""
This modules provides tensors objects extracted from dfpt calculations.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor
from abipy.iotools import ETSF_Reader


class NLOpticalSusceptibilityTensor(Tensor):
    """
    Subclass of :class:`pymatgen.analysis.elasticity.tensors.Tensor` containing the
    non-linear optical susceptibility tensor.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Creates the tensor from a anaddb.nc netcdf file containing ``dchide``.
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


class DielectricTensor(SquareTensor):
    """
    Subclass of :class:`pymatgen.analysis.elasticity.tensors.SquareTensor`
    describing a dielectric tensor.
    """

    def reflectivity(self, tol=1e-6):
        """
        If the tensor is diagonal (with off diagonal elements smaller than tol)
        returns the three components of the reflectivity.
        """

        d = np.diag(self)

        if np.max(np.abs(self-np.diag(d))) > tol:
            raise ValueError("The tensor is not diagonal.")

        d = np.sqrt(d)

        return np.abs((1-d)/(1+d))**2
