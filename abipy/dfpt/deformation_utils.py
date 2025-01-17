# deformation_utils.py
from __future__ import annotations

import numpy as np

from pymatgen.core import Lattice
#from abipy.core.structure import Structure
from abipy.core.symmetries import AbinitSpaceGroup
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def generate_deformations_volumic(structure, eps_V=0.02, scales=None):
    if scales is None:
        scales = [-1, 0, 1, 2, 3]
    rprim = structure.lattice.matrix
    structures_new = {}

    for i in scales:
        rprim2 = np.copy(rprim)
        rprim2[:, :] = rprim[:, :] * (1.00 + eps_V * i)**(1/3.)

        structure2 = structure.copy()
        structure2.lattice = Lattice(rprim2)
        #structure2.scale_lattice(structure2.volume*(1.00 + eps_V * i))
        namei = int(round(1000 * (1.00 + eps_V * i)))
        formatted_namei = f"{namei:04d}"
        structures_new[formatted_namei] = structure2


    return structures_new


#def generate_deformations(structure, eps=0.005) -> tuple:
def generate_deformations(structure, eps: float) -> tuple:
    """
    """
    spgrp = AbinitSpaceGroup.from_structure(structure)
    #print(spgrp)

    spgrp_number = spgrp.spgid
    rprim = structure.lattice.matrix

    rprim2 = np.copy(rprim)
    structures_new = {}
    strain_inds = []

    def _add(name, new_rprim, i, j, k, l, m, n) -> None:
        """Helper function to register a new structure in internal dict."""
        new_structure = structure.copy()
        new_structure.lattice = Lattice(new_rprim)
        structures_new[name] = new_structure
        strain_inds.append([i, j, k, l, m, n])

    i, j, k, l, m, n = 6 * [0]

    if 1 <= spgrp_number <= 2:
        disp=[[1,1,1,1,1,1],  [0,1,1,1,1,1],  [2,1,1,1,1,1],  [1,0,1,1,1,1],  [1,2,1,1,1,1],  [1,1,0,1,1,1],
              [1,1,2,1,1,1],  [1,1,1,0,1,1],  [1,1,1,2,1,1],  [1,1,1,1,0,1],  [1,1,1,1,2,1],  [1,1,1,1,1,0],
              [1,1,1,1,1,2],  [0,0,1,1,1,1],  [1,0,0,1,1,1],  [1,1,0,0,1,1],  [1,1,1,0,0,1],  [1,1,1,1,0,0],
              [0,1,0,1,1,1],  [0,1,1,0,1,1],  [0,1,1,1,0,1],  [0,1,1,1,1,0],  [1,0,1,0,1,1],  [1,0,1,1,0,1],
              [1,0,1,1,1,0],  [1,1,0,1,0,1],  [1,1,0,1,1,0],  [1,1,1,0,1,0] , [0 ,0,0,0,0,0]]
        if abs(rprim[1, 0]) > 1e-9 or abs(rprim[2, 0]) > 1e-9 or abs(rprim[2, 1]) > 1e-9:
            print("Warning: The lattice is oriented such that xz = xy = yz = 0.")
        rprim0 = np.copy(rprim)
        a=rprim[0, :]
        b=rprim[1, :]
        c=rprim[2, :]
        norm_a = np.linalg.norm(a)
        norm_b = np.linalg.norm(b)
        norm_c = np.linalg.norm(c)

        # Compute angles between vectors
        cos_ab = np.dot(a, b) / (norm_a * norm_b)
        cos_ac = np.dot(a, c) / (norm_a * norm_c)
        cos_bc = np.dot(b, c) / (norm_b * norm_c)

        rprim0[0,0] = 1.0
        rprim0[0,1] = 0.0
        rprim0[0,2] = 0.0
        rprim0[1,0] = cos_ab
        rprim0[1,1] = np.sqrt(1-cos_ab**2)
        rprim0[1,2] = 0.0
        rprim0[2,0] = cos_ac
        rprim0[2,1] = (cos_bc-rprim0[1,0]*rprim0[2,0])/rprim0[1,1]
        rprim0[2,2] = np.sqrt(1.0-rprim0[2,0]**2-rprim0[2,1]**2)
        rprim0[0,:] = rprim0[0,:]*norm_a
        rprim0[1,:] = rprim0[1,:]*norm_b
        rprim0[2,:] = rprim0[2,:]*norm_c
        print("Old rprim:")
        print(rprim)
        print("New rprim:")
        print(rprim0)

        for pair in disp:
            i, j, k, l, m, n = pair
            rprim2[ :,0] = rprim0[ :,0] * (1.00 + eps * i) + rprim0[ :,1] * (eps * l) +rprim0[ :,2] * (eps * m)
            rprim2[ :,1] = rprim0[ :,1] * (1.00 + eps * j) + rprim0[ :,2] * (eps * n)
            rprim2[ :,2] = rprim0[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            namel = int(round(1000 * (1.00 + eps * l)))
            namem = int(round(1000 * (1.00 + eps * m)))
            namen = int(round(1000 * (1.00 + eps * n)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}_{namel:04d}_{namem:04d}_{namen:04d}"

            _add(formatted_namei, rprim2, i, j, k, l, m, n)

    elif 3 <= spgrp_number <= 15:
        disp=[[1,1,1,1], [0,1,1,1], [2,1,1,1], [1,0,1,1], [1,2,1,1], [1,1,0,1], [1,1,2,1], [1,1,1,0],
              [1,1,1,2], [0,0,1,1], [1,0,0,1], [1,1,0,0], [0,1,0,1], [1,0,1,0], [0,1,1,0]]
        if abs(rprim[1, 0]) > 1e-9 or abs(rprim[0, 1]) > 1e-9 or abs(rprim[2, 1]) > 1e-9 or abs(rprim[1, 2]) > 1e-9:
            print("Error: Monoclinic structure with yx=xy=0 and yz=zy=0 lattice required.")
        elif abs(rprim[0, 2]) > 1e-9 :
            print("Warning: The lattice is oriented such that xz = 0.")
            rprim0 = np.copy(rprim)
            a=rprim[0, :]
            b=rprim[1, :]
            c=rprim[2, :]
            norm_a = np.linalg.norm(a)
            norm_b = np.linalg.norm(b)
            norm_c = np.linalg.norm(c)

            # Compute angles between vectors
            cos_ab = np.dot(a, b) / (norm_a * norm_b)
            cos_ac = np.dot(a, c) / (norm_a * norm_c)
            cos_bc = np.dot(b, c) / (norm_b * norm_c)

            rprim0[0,0] = norm_a
            rprim0[0,2] = 0.0
            rprim0[1,1] = norm_b
            rprim0[2,0] = norm_c*cos_ac
            rprim0[2,2] = norm_c*np.sqrt(1-cos_ac**2)
        print("Old rprim:")
        print(rprim)
        print("New rprim:")
        print(rprim0)

        for pair in disp:
            i, j, k, l = pair
            rprim2[ :,0] = rprim0[ :,0] * (1.00 + eps * i) +rprim0[ :,2] * (eps * l)
            rprim2[ :,1] = rprim0[ :,1] * (1.00 + eps * j)
            rprim2[ :,2] = rprim0[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            namel = int(round(1000 * (1.00 + eps * l)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}_{namel:04d}"

            _add(formatted_namei, rprim2, i, j, k, l, m, n)

    elif 16 <= spgrp_number <= 74:
        disp=[[0,0,1],[0,1,0],[1,0,0],[1,1,1],[0,1,1],[2,1,1],[1,0,1],[1,2,1],[1,1,0],[1,1,2]]
        for pair in disp:
            i, j, k = pair
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * j)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}"

            _add(formatted_namei, rprim2, i, j, k, l, m, n)

    elif 75 <= spgrp_number <= 194:
        disp=[[0,0],[1,1],[0,1],[2,1],[1,0],[1,2]]
        for pair in disp:
            i, k = pair
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * i)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namek = int(round(1000 * (1.00 + eps * k)))
            formatted_namei = f"{namei:04d}_{namek:04d}"

            _add(formatted_namei, rprim2, i, j, k, l, m, n)

    elif 195 <= spgrp_number <= 230:
        for i in range(3):
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * i)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * i)
            namei = int(round(1000 * (1.00 + eps * i)))
            formatted_namei = f"{namei:04d}"

            _add(formatted_namei, rprim2, i, j, k, l, m, n)

    else:
        raise ValueError(f"Invalid {spgrp_number=}")

    return structures_new, np.array(strain_inds, dtype=int), spgrp_number

