# coding: utf-8
"""Work subclasses related to effective mass calculations."""

import numpy as np

from abipy.core.kpoints import Kpoint
from .works import Work


class EffectiveMassLineWork(Work):
    """
    Work for the computation of effective masses via finite differences along a k-line.
    """

    @classmethod
    def from_scf_input(cls, scf_input, kpoint, step, npts=9, denpath=None, red_dirs=None, manager=None):
        """
        Build a Work from an |AbinitInput| representing a SCF-GS calculation.

        Args:
            scf_input: |AbinitInput| for SCF-GS used as template to generate the other inputs.
            kpoint: Reduced coordinates of the k-point where effective masses are wanted.
            step:
            npts: Number of points for central finite difference (3, 5, 7, 9).
            red_dirs:
            manager: |TaskManager| instance. Use default if None.
        """
        allowed_npts = {3, 5, 7, 9}
        if npts not in allowed_npts:
            raise ValueError("Number of points: `%s` should be in: `%s`" % (npts, str(allowed_npts)))

        reciprocal_lattice = scf_input.structure.lattice.reciprocal_lattice
        kpoint = Kpoint.as_kpoint(kpoint, reciprocal_lattice)
        red_dirs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] if red_dirs is None else np.reshape(red_dirs, (-1, 3))

        kpts = []
        for rdir in red_dirs:
            bvers = reciprocal_lattice.matrix.T @ rdir
            bvers /= np.sqrt(np.dot(bvers, bvers))
            k0 = kpoint.cart_coords - bvers * (npts // 2) * step
            for ii in range(npts):
                kc = k0 + ii * bvers * step
                kpts.append(kc)

        kpts = reciprocal_lattice.get_fractional_coords(np.reshape(kpts, (-1, 3)))
        print("List of kpoints:\n", kpts)

        nscf_input = scf_input.make_nscf_kptopt0(kpts)
        new = cls(manager=manager)

        if denpath is None:
            scf_task = new.register_scf_task(scf_input)
        else:
            scf_task = os.path.abspath(denpath)

        new.register_nscf_task(nscf_input, deps={scf_task: "DEN"})

        return new
