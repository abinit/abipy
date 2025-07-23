# coding: utf-8
"""
This modules provides subclasses of pymatgen tensor objects.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from pymatgen.core.tensors import Tensor, SquareTensor
from pymatgen.analysis.elasticity.elastic import ElasticTensor  # noqa: F401
from pymatgen.analysis.elasticity.stress import Stress as pmg_Stress
from pymatgen.analysis.piezo import PiezoTensor # noqa: F401
from abipy.iotools import ETSF_Reader


class _Tensor33:

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.get_dataframe()._repr_html_()

    def get_dataframe(self, tol=1e-3, cmode=None) -> pd.DataFrame:
        """
        Return |pandas-Dataframe| with tensor elements set to zero below `tol`.

        Args:
            cmode: "real" or "imag" to include only the real/imaginary part.
        """
        tensor = self.zeroed(tol=tol)
        if cmode == "real": tensor = tensor.real
        if cmode == "imag": tensor = tensor.imag

        return pd.DataFrame({"x": tensor[:,0], "y": tensor[:,1], "z": tensor[:,2]}, index=["x", "y", "z"])

    def get_voigt_dataframe(self, tol=1e-3) -> pd.DataFrame:
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

    def reflectivity(self, n1=1, tol=1e-6) -> pd.DataFrame:
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


class DielectricDataList(list):
    """
    A list of tuples with (DielectricTensor, structure, params)
    Useful for convergence studies.

    Example:

        diel_data = DielectricDataList()
        diel_data.append((eps0, structure0, params0))
        diel_data.append((eps1, structure1, params1))

        df = diel_data.get_dataframe()
    """

    def append(self, obj) -> None:
        """Extend append method with validation logic."""
        if not isinstance(obj, (list, tuple)):
            raise TypeError(f"Expecting list or tuple  but got {type(obj)=}")

        if len(obj) != 3:
            raise TypeError(f"Expecting two items but got {len(obj)=}")

        if not isinstance(obj[0], DielectricTensor):
            if isinstance(obj[0], np.ndarray):
                obj = list(obj)
                obj[0] = DielectricTensor(obj[0])
            else:
                raise TypeError(f"Expecting DielectricTensor instance but got {type(obj[0])=}")

        from abipy.core.structure import Structure
        if not isinstance(obj[1], Structure):
            raise TypeError(f"Expecting Structure instance but got {type(obj[1])=}")

        if not isinstance(obj[2], dict):
            raise TypeError(f"Expecting dict instance but got {type(obj[2])=}")

        return super().append(obj)

    @property
    def eps_list(self) -> list:
        return [obj[0] for obj in self]

    @property
    def structures(self) -> list:
        return [obj[1] for obj in self]

    @property
    def params_list(self) -> list[dict]:
        return [obj[2] for obj in self]

    def has_same_structure(self) -> bool:
        """True if all structures are equal."""
        if len(self) in (0, 1): return True
        structures = self.structures
        structure0 = structures[0]
        return all(structure0 == s for s in structures[1:])

    def __str__(self):
        return self.get_dataframe().to_string()

    def get_dataframe(self, with_params=True, with_geo=False, with_spglib=True, **kwargs) -> pd.DataFrame:
        """
        Dataframe with the components of eps_infinity.

        Args:
            with_params: True to add parameters.
            with_geo: True to add info on structure.
            with_params: True to add calculations parameters.
            kwargs: Optional kwargs passed to add_geo_and_params.
        """
        structures, eps_list, params_list = self.structures, self.eps_list, self.params_list

        comps2inds = {"xx": (0,0), "yy": (1,1), "zz": (2,2),
                      "xy": (0, 1), "xz": (0, 2), "yx": (1, 0), "yz": (1, 2), "zx": (2, 0), "zy": (2, 1)}

        rows = []
        for structure, eps, params in zip(structures, eps_list, params_list, strict=True):
            d = {}
            for k, ind in comps2inds.items():
                d[k] = eps[ind]

            if with_params:
                d.update(params)

            if with_geo:
                geo_dict = structure.get_dict4pandas(with_spglib=with_spglib, **kwargs)
                d.update(geo_dict)

            rows.append(d)

        return pd.DataFrame(rows)


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
    def from_file(cls, filepath: str) -> NLOpticalSusceptibilityTensor:
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
