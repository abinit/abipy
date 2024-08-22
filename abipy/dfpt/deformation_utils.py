# deformation_utils.py

import numpy as np
from pymatgen.core import Structure, Lattice, Element
from abipy.core.symmetries import AbinitSpaceGroup
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import re


def generate_deformations_volumic(structure, eps_V=0.02, scales=[-1, 0, 1, 2, 3]):
    spgrp = AbinitSpaceGroup.from_structure(structure)
    spgrp_number = spgrp.spgid
    rprim = structure.lattice.matrix

    structures_new = {}

    for i in scales:
        rprim2 = np.copy(rprim)
        rprim2[:, :] = rprim[:, :] * (1.00 + eps_V * i)**(1/3.)
        namei = int(round(1000 * (1.00 + eps_V * i)))
        formatted_namei = f"{namei:04d}"

        structure2 = structure.copy()
        structure2.lattice = Lattice(rprim2)
        structures_new[formatted_namei] = structure2
        print(formatted_namei)
        print(rprim2)

    return structures_new

def generate_deformations(structure , eps=0.005):
    spgrp = AbinitSpaceGroup.from_structure(structure ) 
    spgrp_number=spgrp.spgid
    rprim= structure.lattice.matrix

    rprim2 = np.copy(rprim)
    rprim_new = {}
    structures_new = {}

    if 1 <= spgrp_number <= 2:
        disp=[[1,1,1,1,1,1],  [0,1,1,1,1,1],  [2,1,1,1,1,1],  [1,0,1,1,1,1],  [1,2,1,1,1,1],  [1,1,0,1,1,1],  
              [1,1,2,1,1,1],  [1,1,1,0,1,1],  [1,1,1,2,1,1],  [1,1,1,1,0,1],  [1,1,1,1,2,1],  [1,1,1,1,1,0],  
              [1,1,1,1,1,2],  [0,0,1,1,1,1],  [1,0,0,1,1,1],  [1,1,0,0,1,1],  [1,1,1,0,0,1],  [1,1,1,1,0,0],  
              [0,1,0,1,1,1],  [0,1,1,0,1,1],  [0,1,1,1,0,1],  [0,1,1,1,1,0],  [1,0,1,0,1,1],  [1,0,1,1,0,1],  
              [1,0,1,1,1,0],  [1,1,0,1,0,1],  [1,1,0,1,1,0],  [1,1,1,0,1,0] , [0 ,0,0,0,0,0]]  
        #if abs(rprim[1, 0]) > 1e-9 or abs(rprim[2, 0]) > 1e-9 or abs(rprim[2, 1]) > 1e-9:
        print("Warning: The lattice is oriented such that xz =xy =yz =0 .")
        a=rprim[0, :]
        b=rprim[1, :]
        c=rprim[2, :]
        print(a,b,c)
        norm_a = np.linalg.norm(a)
        norm_b = np.linalg.norm(b)
        norm_c = np.linalg.norm(c)

        # Compute angles between vectors
        cos_ab = np.dot(a, b) / (norm_a * norm_b)
        cos_ac = np.dot(a, c) / (norm_a * norm_c)
        cos_bc = np.dot(b, c) / (norm_b * norm_c)

        rprim[0,0] = 1.0 
        rprim[0,1] = 0.0 
        rprim[0,2] = 0.0 
        rprim[1,0] = cos_ab
        rprim[1,1] = np.sqrt(1-cos_ab**2)
        rprim[1,2] = 0.0 
        rprim[2,0] = cos_ac
        rprim[2,1] = (cos_bc-rprim[1,0]*rprim[2,0])/rprim[1,1]
        rprim[2,2] = np.sqrt(1.0-rprim[2,0]**2-rprim[2,1]**2)
        rprim[0,:] = rprim[0,:]*norm_a
        rprim[1,:] = rprim[1,:]*norm_b
        rprim[2,:] = rprim[2,:]*norm_c
        print("New rprim:")
        print(rprim)

        for pair in disp:
            i,j,k,l,m,n = pair
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i) + rprim[ :,1] * (eps * l) +rprim[ :,2] * (eps * m)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * j) + rprim[ :,2] * (eps * n)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            namel = int(round(1000 * (1.00 + eps * l)))
            namem = int(round(1000 * (1.00 + eps * m)))
            namen = int(round(1000 * (1.00 + eps * n)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}_{namel:04d}_{namem:04d}_{namen:04d}"

            structure2=structure
            structure2.lattice=Lattice(rprim2)
            structures_new[formatted_namei] = structure2 
            print (formatted_namei)
            print (rprim2)

        return structures_new
    elif 3 <= spgrp_number <= 15:
        disp=[[1,1,1,1], [0,1,1,1], [2,1,1,1], [1,0,1,1], [1,2,1,1], [1,1,0,1], [1,1,2,1], [1,1,1,0],
              [1,1,1,2], [0,0,1,1], [1,0,0,1], [1,1,0,0], [0,1,0,1], [1,0,1,0], [0,1,1,0]]
        if abs(rprim[1, 0]) > 1e-9 or abs(rprim[0, 1]) > 1e-9 or abs(rprim[2, 1]) > 1e-9 or abs(rprim[1, 2]) > 1e-9:
            print("Error: Monoclinic structure with yx=xy=0 and yz=zy=0 lattice required.")
        elif abs(rprim[0, 2]) > 1e-9 :
            print("Warning: The lattice is oriented such that xz = 0.")
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

            rprim[0,0] = norm_a
            rprim[0,2] = 0.0 
            rprim[1,1] = norm_b
            rprim[2,0] = norm_c*cos_ac
            rprim[2,2] = norm_c*np.sqrt(1-cos_ac**2)
        print("New rprim:")
        print(rprim)

        for pair in disp:
            i,j,k,l = pair
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i) +rprim[ :,2] * (eps * l)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * j)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            namel = int(round(1000 * (1.00 + eps * l)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}_{namel:04d}"

            structure2=structure
            structure2.lattice=Lattice(rprim2)
            structures_new[formatted_namei] = structure2 
            print (formatted_namei)
            print (rprim2)

        return structures_new
    elif 16 <= spgrp_number <= 74:
        disp=[[0,0,1],[0,1,0],[1,0,0],[1,1,1],[0,1,1],[2,1,1],[1,0,1],[1,2,1],[1,1,0],[1,1,2]]
        for pair in disp:
            i,j,k = pair
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * j)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * k)

            namei = int(round(1000 * (1.00 + eps * i)))
            namej = int(round(1000 * (1.00 + eps * j)))
            namek = int(round(1000 * (1.00 + eps * k)))
            formatted_namei = f"{namei:04d}_{namej:04d}_{namek:04d}"

            structure2=structure
            structure2.lattice=Lattice(rprim2)
            structures_new[formatted_namei] = structure2 
            print (formatted_namei)
            print (rprim2)

        return structures_new
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
            rprim_new[formatted_namei] = rprim2

            structure2=structure
            structure2.lattice=Lattice(rprim2)
            structures_new[formatted_namei] = structure2 
           # print (formatted_namei)
           # print (rprim2)

        return structures_new
    elif 195 <= spgrp_number <= 230:
        for i in range(3):
            rprim2[ :,0] = rprim[ :,0] * (1.00 + eps * i)
            rprim2[ :,1] = rprim[ :,1] * (1.00 + eps * i)
            rprim2[ :,2] = rprim[ :,2] * (1.00 + eps * i)
            namei = int(round(1000 * (1.00 + eps * i)))
            formatted_namei = f"{namei:04d}"

            structure2=structure
            structure2.lattice=Lattice(rprim2)
            structures_new[formatted_namei] = structure2 
           # print (formatted_namei)
           # print (rprim2)
        return structures_new

