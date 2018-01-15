# coding: utf-8
"""
AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from monty.functools import lazy_property
from monty.string import marquee
from monty.termcolor import cprint
from abipy.core.tensor import Tensor
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.iotools import ETSF_Reader
from abipy.dfpt.phonons import InteratomicForceConstants
from abipy.dfpt.ddb import Becs
from abipy.dfpt.tensors import NLOpticalSusceptibilityTensor


class AnaddbNcFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
    This object is usually instanciated with `abiopen("anaddb.nc")`.

    .. attribute:: structure

        |Structure| object.

    .. attribute:: emacro

        Macroscopic dielectric tensor. None if the file does not contain this information.

    .. attribute:: becs

        Born effective charges. None if the file does not contain this inf

    .. attribute:: ifc

        :class:`InteratomicForceConstants` object with the interatomic force constants calculated by anaddb.
        None, if the netcdf file does not contain the IFCs.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AnaddbNcFile
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(AnaddbNcFile, self).__init__(filepath)
        self.reader = ETSF_Reader(filepath)
        self._structure = self.reader.read_structure()

    def close(self):
        self.reader.close()

    @lazy_property
    def params(self):
        return {}

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """
        String representation

        Args:
            verbose: verbosity level.
        """
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        return "\n".join(lines)

    @property
    def structure(self):
        return self._structure

    @lazy_property
    def emacro(self):
        """
        Macroscopic dielectric tensor. None if the file does not contain this information.
        """
        try:
            return Tensor.from_cartesian_tensor(self.reader.read_value("emacro_cart"),
                                                self.structure.lattice, space="r")
        except Exception as exc:
            print(exc, "Returning None", sep="\n")
            return None

    @lazy_property
    def emacro_rlx(self):
        """
        Relaxed ion Macroscopic dielectric tensor. None if the file does not contain this information.
        """
        try:
            return Tensor.from_cartesian_tensor(self.reader.read_value("emacro_cart_rlx"),
                                                self.structure.lattice, space="r")
        except Exception as exc:
            print(exc, "Requires dieflag > 0", "Returning None", sep="\n")
            return None

    @lazy_property
    def becs(self):
        """
        Born effective charges. None if the file does not contain this information.
        """
        try:
            chneut = -666 # TODO: anaddb.nc should contain the input file.
            return Becs(self.reader.read_value("becs_cart"), self.structure, chneut=chneut, order="f")
        except Exception as exc:
            print(exc, "Returning None", sep="\n")
            return None

    @lazy_property
    def ifc(self):
        """
        The interatomic force constants calculated by anaddb.
        The following anaddb variables should be used in the run: ``ifcflag``, ``natifc``, ``atifc``, ``ifcout``.
        Return None, if the netcdf_ file does not contain the IFCs,
        """
        try:
            return InteratomicForceConstants.from_file(self.filepath)
        except Exception as exc:
            print(exc)
            cprint("Interatomic force constants have not been calculated. Returning None", "red")
            return None

    @lazy_property
    def dchide(self):
        """
        Non-linear optical susceptibility tensor.
        Returns a :class:`NLOpticalSusceptibilityTensor` or None if the file does not contain this information.
        """
        try:
            return NLOpticalSusceptibilityTensor(self.reader.read_value("dchide"))
        except Exception as exc:
            print(exc, "Requires nlflag > 0", "Returning None", sep="\n")
            return None

    @lazy_property
    def dchidt(self):
        """
        First-order change in the linear dielectric susceptibility.
        Returns a list of lists of 3x3 Tensor object with shape (number of atoms, 3).
        The [i][j] element of the list contains the Tensor representing the change due to the
        displacement of the ith atom in the jth direction.
        None if the file does not contain this information.
        """
        try:
            a = self.reader.read_value("dchidt").T
            dchidt = []
            for i in a:
                d = []
                for j in i:
                    d.append(Tensor.from_cartesian_tensor(j, self.structure.lattice, space="r"))
                dchidt.append(d)

            return dchidt
        except Exception as exc:
            print(exc, "Requires 0 < nlflag < 3", "Returning None", sep="\n")
            return None

    @lazy_property
    def oscillator_strength(self):
        """
        A complex |numpy-array| containing the oscillator strengths with shape (number of phonon modes, 3, 3),
        in a.u. (1 a.u.=253.2638413 m3/s2).
        None if the file does not contain this information.
        """
        try:
            return self.reader.read_value("oscillator_strength", cmode="c")
        except Exception as exc:
            print(exc, "Oscillator strengths require dieflag == 1, 3 or 4", "Returning None", sep="\n")
            raise
            return None

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ananc = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ananc)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)
