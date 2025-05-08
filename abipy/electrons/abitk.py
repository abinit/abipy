# coding: utf-8
"""Interface to the abitk Fortran executable."""
from __future__ import annotations

import os
import tempfile
import numpy as np

from functools import cached_property
from monty.string import marquee
from abipy.tools.numtools import is_diagonal
from abipy.core.mixins import AbinitNcFile, Has_Structure
from abipy.core.structure import Structure
from abipy.electrons.ebands import ElectronsReader
from abipy.flowtk.wrappers import Abitk


class KmeshFile(AbinitNcFile, Has_Structure):
    """
    Interface to the KMESH.nc file produced by abitk

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: KmeshFile
    """

    @classmethod
    def from_ngkpt_shifts(cls, ncpath_with_structure, ngkpt, shifts, kptopt=3, chksymbreak=0, verbose=0) -> Kmesh:
        """
        Args:
            ncpath_with_structure:
            ngkpt:
            shifts:
            kptopt:
            chksymbreak:
        """
        shiftk = np.reshape(shifts, (-1, 3))
        nshifk = len(shifts)

        def s(numbers):
            return ', '.join(str(n) for n in np.array(numbers).flatten())

        import shutil
        workdir = tempfile.mkdtemp()
        base = os.path.basename(ncpath_with_structure)
        target_path = os.path.join(workdir, base)
        shutil.copy(ncpath_with_structure, os.path.join(workdir, base))
        #print("workdir", workdir)

        exec_args = ["ibz", ncpath_with_structure,
                     f"--ngkpt {s(ngkpt)}",
                     f"--shiftk {s(shiftk)}",
                     f"--kptopt {kptopt}",
                     f"--chksymbreak {chksymbreak}"
                    ]
        abitk = Abitk(verbose=verbose)
        abitk.run(exec_args, workdir=workdir)

        return cls(os.path.join(workdir, "KMESH.nc"))

    @classmethod
    def from_file(cls, filepath: str) -> KmeshFile:
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = r = ElectronsReader(filepath)

        # Read dimensions
        self.nkibz = r.read_dimvalue("nkibz")
        self.nkbz = r.read_dimvalue("nkbz")

        # Read IBZ, IBZ
        #self.old_kptrlatt = r.read_value("kptrlatt")
        #self.old_shiftk = r.read_value("shiftk")
        self.kptrlatt = r.read_value("new_kptrlatt")
        self.shiftk = r.read_value("new_shiftk")

        self.kptopt = r.read_value("kptopt")
        self.kibz = r.read_value("kibz")
        self.kbz = r.read_value("kbz")
        self.wtk = r.read_value("wtk")

        # Read mapping and convert indices from Fortran to python
        self.bz2ibz = r.read_value("bz2ibz")
        self.bz2ibz[0] -= 1
        self.bz2ibz[1] -= 1

    def ngkpt_and_shifts(self):
        ngkpt = None if not is_diagonal(self.kptrlatt) else np.diag(self.kptrlatt)
        return ngkpt, self.shiftk

    @cached_property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.r.read_structure()

    @cached_property
    def params(self) -> dict:
        return {}

    def close(self) -> None:
        """Called at the end of the ``with`` context manager."""
        return self.r.close()

    def __str__(self) -> str:
        """String representation"""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app("")

        app(f"kptopt: {self.kptopt}")
        app(f"nkibz: {self.nkibz}")
        app(f"nkbz: {self.nkbz}")
        ngkpt, _ = self.ngkpt_and_shifts()
        if ngkpt is not None:
            app(f"ngkpt: {ngkpt}")
        else:
            app(f"kptrlatt: {np.array_str(self.kptrlatt, precision=3)}")

        app(f"shiftk: {self.shiftk}")

        if verbose > 1:
            app("")
            app(f"kibz: {self.kibz}")
            app(f"kbz: {self.kbz}")
            app(f"wtk: {np.array_str(self.wtk, precision=3)}")
            app(f"bz2ibz: {self.bz2ibz}")

        return "\n".join(lines)
