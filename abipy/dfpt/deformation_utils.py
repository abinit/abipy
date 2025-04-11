# deformation_utils.py
from __future__ import annotations

import numpy as np

from pymatgen.core import Lattice
from abipy.core.symmetries import AbinitSpaceGroup
from abipy import abilab


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
def generate_deformations(structure, eps: float, str_type ='BO', eps_ref = [0.005, 0.005 ,0.005], mode = "TEC" )-> tuple:
    """
    Generates deformed structures by applying strain to the input structure's lattice.

    Parameters:
    structure : pymatgen.Structure
        The input crystal structure.
    eps : float
        Strain magnitude to be applied to the lattice. Strains will be applied to the reference lattice. 
    str_type : str, optional
        Structure types :
        -  ref': The reference structure for the Taylor expansion method.
        - 'BO': BO structure, applies strain scaling using eps_ref to generate a
                deformed reference structure for the Taylor expansion method.
    eps_ref : list of float, optional
        A list of strain values corresponding to normal strains along xx, yy, and zz
        directions ([eps_xx, eps_yy, eps_zz]). Default: [0.005, 0.005, 0.005](not used in 'ref' structure type).
    mode : str, optional
        Determines the purpose of deformation:
        - 'TEC': Generates new structures needed to compute the thermal expansion coefficient.
        - 'ECs': Generates structures for both elastic constants and thermal expansion coefficient calculations.

    Returns:
    tuple
        A tuple containing the modified structures and associated strain indices.
    """
    spgrp = AbinitSpaceGroup.from_structure(structure)

    spgrp_number = spgrp.spgid
    rprim = np.copy(structure.lattice.matrix)
    #TOFIX: The primitive lattice is not a standard one. 
    #primitive_structure = structure.get_primitive_structure(tolerance=0.01, use_site_props=False, constrain_latt=False)
    #spgrp = AbinitSpaceGroup.from_structure(primitive_structure)
    #print(spgrp)

    angdeg    = structure.lattice.angles
    lattice_a = structure.lattice.abc[0]
    lattice_b = structure.lattice.abc[1]
    lattice_c = structure.lattice.abc[2]
    tol8 = 1.0e-8
    # Rotate lattice parameters to follow Abinit conventions.
    # Keep the angles unchanged if the lattice is orthogonal (90°)
    # or if it belongs to a cubic system with a primitive cell having 60° angles.
    if (not((abs(angdeg[0] - 60) + abs(angdeg[1] - 60) + abs(angdeg[2] - 60)) < tol8 ) and 
        not((abs(angdeg[0] - 90) + abs(angdeg[1] - 90) + abs(angdeg[2] - 90)) < tol8 )):
        if (abs(angdeg[0] - angdeg[1]) < tol8 and abs(angdeg[1] - angdeg[2]) < tol8):
            # Trigonal symmetry case
            cosang = np.cos(np.pi * angdeg[0] / 180.0)
            a2 = (2.0 / 3.0) * (1.0 - cosang)
            aa = np.sqrt(a2)
            cc = np.sqrt(1.0 - a2)

            rprim0 = np.array([
                [ aa,          0.0,          cc],
                [-0.5 * aa,  np.sqrt(3.0) * 0.5 * aa, cc],
                [-0.5 * aa, -np.sqrt(3.0) * 0.5 * aa, cc]
            ])
        else:
            # General case
            rprim0 = np.zeros((3, 3))
            rprim0[0, 0] = 1.0
            rprim0[1, 0] = np.cos(np.pi * angdeg[2] / 180.0)
            rprim0[1, 1] = np.sin(np.pi * angdeg[2] / 180.0)
            rprim0[2, 0] = np.cos(np.pi * angdeg[1] / 180.0)
            rprim0[2, 1] = (np.cos(np.pi * angdeg[0] / 180.0) - rprim0[1, 0] * rprim0[2, 0]) / rprim0[1, 1]
            rprim0[2, 2] = np.sqrt(1.0 - rprim0[2, 0]**2 - rprim0[2, 1]**2)
        rprim0[0,:] = rprim0[0,:]*lattice_a
        rprim0[1,:] = rprim0[1,:]*lattice_b
        rprim0[2,:] = rprim0[2,:]*lattice_c
        print("Old rprim:")
        print(rprim)
        print("New rprim:")
        print(rprim0)
    else:
        rprim0=rprim

    if str_type == 'BO':
        rprim_BO = np.copy(rprim)
        # Scale each lattice vector by the corresponding strain component
        # to generate the reference structure.
        rprim[ :,0] *= (1.00 + eps_ref[0])
        rprim[ :,1] *= (1.00 + eps_ref[1])
        rprim[ :,2] *= (1.00 + eps_ref[2])
    elif str_type != 'ref':
        raise ValueError("Invalid method. Choose 'ref' or 'BO'.")

    rprim2 = np.copy(rprim)
    structures_new = {}
    strain_inds = []
    rprim0 = np.copy(rprim)
    rprim=rprim0

    def _add(name, new_rprim, i, j, k, l, m, n) -> None:
        """Helper function to register a new structure in internal dict."""
        new_structure = structure.copy()
        new_structure.lattice = Lattice(new_rprim)
        structures_new[name] = new_structure
        strain_inds.append([i, j, k, l, m, n])

    # Initialize strain components in Voigt notation: xx, yy, zz, yz, xz, xy
    i, j, k, l, m, n = 6 * [0]

    if 1 <= spgrp_number <= 2: # Check for triclinic crystal systems
        # Define strain configurations in Voigt notation
        disp=[[0,0,0,0,0,0], [-1,0,0,0,0,0], [1,0,0,0,0,0], [0,-1,0,0,0,0], [0,1,0,0,0,0], [0,0,-1,0,0,0],
              [0,0,1,0,0,0], [0,0,0,-1,0,0], [0,0,0,1,0,0], [0,0,0,0,-1,0], [0,0,0,0,1,0], [0,0,0,0,0,-1],
              [0,0,0,0,0,1], [-1,-1,0,0,0,0], [0,-1,-1,0,0,0], [0,0,-1,-1,0,0], [0,0,0,-1,-1,0], [0,0,0,0,-1,-1],
              [-1,0,-1,0,0,0], [-1,0,0,-1,0,0], [-1,0,0,0,-1,0], [-1,0,0,0,0,-1], [0,-1,0,-1,0,0], [0,-1,0,0,-1,0],
              [0,-1,0,0,0,-1], [0,0,-1,0,-1,0], [0,0,-1,0,0,-1], [0,0,0,-1,0,-1]]

    elif 3 <= spgrp_number <= 15: # Check for triclinic crystal systems
        # Define strain configurations in Voigt notation
        disp=[[0,0,0,0,0,0], [-1,0,0,0,0,0], [1,0,0,0,0,0], [0,-1,0,0,0,0], [0,1,0,0,0,0], [0,0,-1,0,0,0],
              [0,0,1,0,0,0], [0,0,0,0,-1,0], [0,0,0,0,1,0], [-1,-1,0,0,0,0], [0,-1,-1,0,0,0], [0,0,-1,0,-1,0],
              [-1,0,-1,0,0,0], [0,-1,0,0,-1,0], [-1,0,0,0,-1,0]]
        if mode=="ECs":
            disp.extend([[0,0,0,-1,0,0], [0,0,0,0,0,-1], [0,0,0,-1,0,-1]])

    elif 16 <= spgrp_number <= 74: # Check for orthorhombic crystal systems
        disp=[[0,0,0,0,0,0], [-1,0,0,0,0,0], [1,0,0,0,0,0], [0,-1,0,0,0,0], [0,1,0,0,0,0], [0,0,-1,0,0,0],
              [0,0,1,0,0,0], [-1,-1,0,0,0,0], [0,-1,-1,0,0,0], [-1,0,-1,0,0,0]]
        if mode=="ECs":
            disp.extend([[0,0,0,1,0,0], [0,0,0,2,0,0], [0,0,0,0,1,0], [0,0,0,0,2,0], [0,0,0,0,0,1], [0,0,0,0,0,2]])

    elif 75 <= spgrp_number <= 194: # Check for uniaxial crystal systems
        if mode=="TEC":
            disp=[[0,0,0,0,0,0], [0,0,-1,0,0,0], [0,0,1,0,0,0], [-1,-1,0,0,0,0], [1,1,0,0,0,0], [-1,-1,-1,0,0,0]]
        else:
            disp=[[0,0,0,0,0,0], [-1,0,0,0,0,0], [1,0,0,0,0,0], [0,0,-1,0,0,0],
                  [0,0,1,0,0,0], [-1,-1,0,0,0,0], [-1,0,-1,0,0,0]]
            if    75  <= spgrp_number <= 142: # Check for Tetragonal crystal systems
                disp.extend([[0,0,0,1,0,0], [0,0,0,2,0,0], [0,0,0,0,0,1], [0,0,0,0,0,2]])
            elif  143 <= spgrp_number <= 167: # Check for trigonal systems
                disp.extend([[0,0,0,1,0,0], [0,0,0,2,0,0],[-1,0,0,1,0,0]])
            elif  168 <= spgrp_number <= 194: # Check for hexagonal crystal systems
                disp.extend([[0,0,0,1,0,0], [0,0,0,2,0,0]])

    elif 195 <= spgrp_number <= 230: # Check for cubic crystal systems
        if mode=="TEC":
            disp=[[0,0,0,0,0,0], [-1,-1,-1,0,0,0], [1,1,1,0,0,0]]
        else:
            disp=[[0,0,0,0,0,0], [-1,0,0,0,0,0], [1,0,0,0,0,0], [-1,-1,0,0,0,0], [0,0,0,1,0,0], [0,0,0,2,0,0]]

    else:
        raise ValueError(f"Invalid {spgrp_number=}")

    for pair in disp:
        i, j, k, l, m, n = pair
        rprim2[ :,0] = rprim0[ :,0] * (1.00 + eps * i) + rprim0[ :,1] * (eps * n) +rprim0[ :,2] * (eps * m)
        rprim2[ :,1] = rprim0[ :,1] * (1.00 + eps * j) + rprim0[ :,2] * (eps * l)
        rprim2[ :,2] = rprim0[ :,2] * (1.00 + eps * k)

        namei = ( eps * i)
        namej = ( eps * j)
        namek = ( eps * k)
        namel = ( eps * l)
        namem = ( eps * m)
        namen = ( eps * n)
        formatted_namei = f"{namei:.3f}_{namej:.3f}_{namek:.3f}_{namel:.3f}_{namem:.3f}_{namen:.3f}"

        _add(formatted_namei, rprim2, i, j, k, l, m, n)
    return structures_new, np.array(strain_inds, dtype=int), spgrp_number

