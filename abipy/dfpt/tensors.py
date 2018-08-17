# coding: utf-8
"""
This modules provides tensors objects produced by DFPT calculations.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd

from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor
from abipy.iotools import ETSF_Reader


class NLOpticalSusceptibilityTensor(Tensor):
    """
    Subclass of |pmg-Tensor| containing the non-linear optical susceptibility tensor.
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


class DielectricTensor(SquareTensor):
    """
    Subclass of |pmg-Tensor| describing a dielectric tensor.
    rank2 symmetric tensor with shape [3, 3].
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

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.get_dataframe()._repr_html_()

    def get_dataframe(self, tol=1e-3):
        """
        Return |pandas-Dataframe| with tensor elements set to zero below `tol`.
        """
        tensor = self.zeroed(tol=tol)
        return pd.DataFrame({"x": tensor[:,0], "y": tensor[:,1], "z": tensor[:,2]}, index=["x", "y", "z"])

    def get_voigt_dataframe(self, tol=1e-3):
        """
        Return |pandas-DataFrame| with Voigt indices as colums (C-indexing starting from 0).
        Useful to analyze the converge of individual elements of the tensor(s)
        Elements below tol are set to zero.
        """
        tensor = self.zeroed(tol=tol)
        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        d = {k: v for k, v in zip(columns, tensor.voigt)}
        return pd.DataFrame(d, index=[0], columns=columns)
