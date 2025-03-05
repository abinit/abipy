"""
Computes thermal stress using EinfVib2EinfVib2 for specific configurations
across various crystallographic structures, from cubic to triclinic.
"""
from __future__ import annotations

import numpy as np
import os
#import abc
import math
import abipy.core.abinit_units as abu

from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.tools.serialization import mjson_load, HasPickleIO
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhdosFile
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.dfpt.vzsisa import anaget_phdoses_with_gauss


class QHA_ZSISA(HasPickleIO):
    """
    Abstract class for the quasi-harmonic approximation analysis.
    Provides some basic methods and plotting utils, plus a converter to write input files for phonopy-qha or to
    generate an instance of phonopy.qha.QHA. These can be used to obtain other quantities and plots.
    Does not include electronic entropic contributions for metals.
    """

    @classmethod
    def from_json_file(cls,
                       filepath: str,
                       nqsmall_or_qppa,
                       anaget_kwargs: dict | None = None,
                       smearing_ev: float | None = None,
                       verbose: int = 0) -> QHA_ZSISA:
        """
        Build an instance from a json file produced by ...

        Args:
            filepath: path to the json file
            anaget_kwargs: kwargs passed to anaget_phdoses_with_gauss.
            smearing_ev:
            verbose: Verbosity level
        """
        data = mjson_load(filepath)
        ddb_relax_paths = data["ddb_relax_paths"]
        gsr_relax_paths = data["gsr_relax_paths"]

        phdos_paths, phbands_paths = anaget_phdoses_with_gauss(nqsmall_or_qppa, smearing_ev,
                                                               ddb_relax_paths, anaget_kwargs, verbose)

        # Create a 6D array with shape (3, 3, 3, 3, 3, 3) initialized with None.
        gsr_paths_6d = np.full((3, 3, 3, 3, 3, 3), None, dtype=object)
        phdos_paths_6d = np.full((3, 3, 3, 3, 3, 3), None, dtype=object)

        strain_inds = data["strain_inds"]
        #print(strain_inds)
        for gsr_path, phdos_path, inds in zip(gsr_relax_paths, phdos_paths, strain_inds, strict=True):
            gsr_paths_6d[inds] = gsr_path
            phdos_paths_6d[inds] = phdos_path

        new = cls.from_files(gsr_paths_6d, phdos_paths_6d, gsr_guess, model='ZSISA')
        #new.pickle_dump(workdir, basename=None)
        return new

    @classmethod
    def from_files(cls, gsr_path_ref, phdos_paths_6D, gsr_path_guess, gsr_path_bo, model = 'ZSISA'):
        """
        Creates an instance of QHA from 6D lists of GSR files and PHDOS.nc files.
        The lists should have the same size and the volumes should match.

        Args:
            gsr_paths_ref: 6D list of paths to GSR files.
            phdos_paths_6D: 6D list of paths to PHDOS.nc files.
            gsr_path_guess: 6D list of paths to GSR files for initial guess.

        Returns: A new instance of QHA
        """
        for gp in gsr_path_guess:
            if os.path.exists(gp):
                with GsrFile.from_file(gp) as g:
                    structure_guess = g.structure
                    stress = g.cart_stress_tensor
            else:
                for gp in gsr_path_ref:
                    with DdbFile.from_file(gp) as g:
                        structure_guess = g.structure
                        stress = g.cart_stress_tensor
        stress_guess = stress/29421.02648438959

        for gp in gsr_path_bo:
            if os.path.exists(gp):
                with DdbFile.from_file(gp) as g:
                    structure_bo = g.structure

        phdos_paths_6D = np.array(phdos_paths_6D)

        current_shape = phdos_paths_6D.shape
        dims_to_add = 6 - len(current_shape)

        if dims_to_add > 0:
            new_shape = current_shape + (1,) * dims_to_add
            phdos_paths_6D = phdos_paths_6D.reshape(new_shape)

        dim = phdos_paths_6D.shape

        structures = []
        phdoses = []

        # Looping through each element in the 6D matrix
        for dim1_idx, dim1_list in enumerate(phdos_paths_6D):
            dim1_doses = []
            dim1_structures = []
            dim1_energies = []
            for dim2_idx, dim2_list in enumerate(dim1_list):
                dim2_doses = []
                dim2_structures = []
                dim2_energies = []
                for dim3_idx, dim3_list in enumerate(dim2_list):
                    dim3_doses = []
                    dim3_structures = []
                    dim3_energies = []
                    for dim4_idx, dim4_list in enumerate(dim3_list):
                        dim4_doses = []
                        dim4_structures = []
                        dim4_energies = []
                        for dim5_idx, dim5_list in enumerate(dim4_list):
                            dim5_doses = []
                            dim5_structures = []
                            dim5_energies = []
                            for path_idx, path in enumerate(dim5_list):
                                if os.path.exists(path):
                                    with PhdosFile(path) as p:
                                        dim5_doses.append(p.phdos)
                                        dim5_structures.append(p.structure)
                                else:
                                    # Handle the case when the file does not exist
                                    dim5_doses.append(None)
                                    dim5_structures.append(None)
                            dim4_doses.append(dim5_doses)
                            dim4_structures.append(dim5_structures)
                        dim3_doses.append(dim4_doses)
                        dim3_structures.append(dim4_structures)
                    dim2_doses.append(dim3_doses)
                    dim2_structures.append(dim3_structures)
                dim1_doses.append(dim2_doses)
                dim1_structures.append(dim2_structures)
            phdoses.append(dim1_doses)
            structures.append(dim1_structures)

        print("dim = ",dim)
        if list(dim) == [3, 3, 1, 1, 1, 1]:
            new_shape = (3, 3, 3, 1, 1, 1)
            dim=[3,3,3,1,1,1]
            structures2 = np.empty(new_shape, dtype=object)  # Use dtype=object for lists
            phdoses2 = np.empty(new_shape, dtype=object)

            # Copy data: Fill the new dim3 by copying dim2 data
            for i in range(3):  # Loop over dim1
                for j in range(3):  # Loop over dim2
                    structures2[i][i][j][0][0][0] = structures[i][j][0][0][0][0]
                    phdoses2[i][i][j][0][0][0] = phdoses[i][j][0][0][0][0]

            structures = structures2
            phdoses = phdoses2
        if list(dim) == [3, 1, 1, 1, 1, 1]:
            new_shape = (3, 3, 3, 1, 1, 1)
            dim=[3,3,3,1,1,1]
            structures2 = np.empty(new_shape, dtype=object)  # Use dtype=object for lists
            phdoses2 = np.empty(new_shape, dtype=object)

            # Copy data: Fill the new dim3 by copying dim2 data
            for i in range(3):  # Loop over dim1
                structures2[i][i][i][0][0][0] = structures[i][0][0][0][0][0]
                phdoses2[i][i][i][0][0][0] = phdoses[i][0][0][0][0][0]

            structures = structures2
            phdoses = phdoses2

        gsr_center = structure_bo
        spgrp = AbinitSpaceGroup.from_structure(structure_bo)
        spgrp_number = spgrp.spgid
        print (spgrp)
        if 1 <= spgrp_number <= 2:
            print ("Triclinic")
            case = "ZSISA_triclinic"
            sym = "triclinic"
        elif 3 <= spgrp_number <= 15: # Check for triclinic crystal systems
            print ("Monoclinic")
            sym = "monoclinic"
        elif 16 <= spgrp_number <= 74: # Check for orthorhombic crystal systems
            print ("Orthorhombic")
            sym = "orthorhombic"
        elif 75  <= spgrp_number <= 142: # Check for Tetragonal crystal systems
            print ("Tetragonal")
            sym = "tetragonal"
        elif 143 <= spgrp_number <= 167: # Check for trigonal systems
            print ("Trigonal")
            sym = "trigonal"
        elif  168 <= spgrp_number <= 194: # Check for hexagonal crystal systems
            print ("Hexagonal")
            sym = "hexagonal"
        elif 195 <= spgrp_number <= 230: # Check for cubic crystal systems
            print ("Cubic")
            sym = "cubic"
        else :
            raise RuntimeError("unknown method")

        return cls(structures, phdoses, dim, structure_guess, stress_guess, structure_bo, sym)

    def __init__(self, structures, phdoses, dim, structure_guess, stress_guess, structure_bo, sym):
        """
        Args:
            structures: list of structures at different volumes.
        """
        self.structures = structures
        self.sym = sym
        self.phdoses = phdoses
        self.dim = dim
        self.HaBohr3_eVA3 = abu.HaBohr3_GPa/abu.eVA3_GPa
        self.eVA3_HaBohr3 = abu.eVA3_GPa/abu.HaBohr3_GPa

        def extract_attribute(structures, attribute_func) -> np.ndarray:
            return np.array([[[[[[attribute_func(s) if s is not None else None for s in col] for col in row] for row in d1]
                            for d1 in d2] for d2 in d3] for d3 in structures])

        self.volumes = extract_attribute(structures, lambda s: s.volume)
        self.lattice_a = extract_attribute(structures, lambda s: s.lattice.abc[0])
        self.lattice_b = extract_attribute(structures, lambda s: s.lattice.abc[1])
        self.lattice_c = extract_attribute(structures, lambda s: s.lattice.abc[2])
        self.alpha = extract_attribute(structures, lambda s: s.lattice.angles[0])
        self.beta = extract_attribute(structures, lambda s: s.lattice.angles[1])
        self.gamma = extract_attribute(structures, lambda s: s.lattice.angles[2])
        self.ax = extract_attribute(structures, lambda s: s.lattice.matrix[0, 0])
        self.ay = extract_attribute(structures, lambda s: s.lattice.matrix[0, 1])
        self.az = extract_attribute(structures, lambda s: s.lattice.matrix[0, 2])
        self.bx = extract_attribute(structures, lambda s: s.lattice.matrix[1, 0])
        self.by = extract_attribute(structures, lambda s: s.lattice.matrix[1, 1])
        self.bz = extract_attribute(structures, lambda s: s.lattice.matrix[1, 2])
        self.cx = extract_attribute(structures, lambda s: s.lattice.matrix[2, 0])
        self.cy = extract_attribute(structures, lambda s: s.lattice.matrix[2, 1])
        self.cz = extract_attribute(structures, lambda s: s.lattice.matrix[2, 2])
        self.ave_x = np.zeros(dim)
        self.ave_y = np.zeros(dim)
        self.ave_z = np.zeros(dim)

        mask = self.ax != None  # Create a mask where self.ax is not None

        # Compute averages using the mask
        self.ave_x[mask] = (abs(self.ax[mask]) + abs(self.bx[mask]) + abs(self.cx[mask])) / 3.0
        self.ave_y[mask] = (abs(self.ay[mask]) + abs(self.by[mask]) + abs(self.cy[mask])) / 3.0
        self.ave_z[mask] = (abs(self.az[mask]) + abs(self.bz[mask]) + abs(self.cz[mask])) / 3.0

        self.structure_guess = structure_guess
        self.volume_guess = structure_guess.volume
        self.lattice_a_guess = structure_guess.lattice.abc[0]
        self.lattice_b_guess = structure_guess.lattice.abc[1]
        self.lattice_c_guess = structure_guess.lattice.abc[2]
        self.matrix_guess = structure_guess.lattice.matrix
        self.stress_guess = stress_guess
        self.frac_coords_guess = structure_guess.frac_coords
        self.angles_guess = structure_guess.lattice.angles
        self.ax_guess = structure_guess.lattice.matrix[0,0]
        self.ay_guess = structure_guess.lattice.matrix[0,1]
        self.az_guess = structure_guess.lattice.matrix[0,2]
        self.bx_guess = structure_guess.lattice.matrix[1,0]
        self.by_guess = structure_guess.lattice.matrix[1,1]
        self.bz_guess = structure_guess.lattice.matrix[1,2]
        self.cx_guess = structure_guess.lattice.matrix[2,0]
        self.cy_guess = structure_guess.lattice.matrix[2,1]
        self.cz_guess = structure_guess.lattice.matrix[2,2]
        self.ave_x_guess = (abs(self.ax_guess)+abs(self.bx_guess)+abs(self.cx_guess))/3.
        self.ave_y_guess = (abs(self.ay_guess)+abs(self.by_guess)+abs(self.cy_guess))/3.
        self.ave_z_guess = (abs(self.az_guess)+abs(self.bz_guess)+abs(self.cz_guess))/3.
        self.ax_bo = structure_bo.lattice.matrix[0,0]
        self.ay_bo = structure_bo.lattice.matrix[0,1]
        self.az_bo = structure_bo.lattice.matrix[0,2]
        self.bx_bo = structure_bo.lattice.matrix[1,0]
        self.by_bo = structure_bo.lattice.matrix[1,1]
        self.bz_bo = structure_bo.lattice.matrix[1,2]
        self.cx_bo = structure_bo.lattice.matrix[2,0]
        self.cy_bo = structure_bo.lattice.matrix[2,1]
        self.cz_bo = structure_bo.lattice.matrix[2,2]
        self.ave_x_bo = (abs(self.ax_bo)+abs(self.bx_bo)+abs(self.cx_bo))/3.
        self.ave_y_bo = (abs(self.ay_bo)+abs(self.by_bo)+abs(self.cy_bo))/3.
        self.ave_z_bo = (abs(self.az_bo)+abs(self.bz_bo)+abs(self.cz_bo))/3.

    def stress_v_ZSISA(self, temp, pressure) -> tuple:

        e,S = self.get_vib_free_energies(temp)

        v = self.volume_guess
        dv = self.volumes[0,0,0,0,0,0]-self.volumes[1,0,0,0,0,0]

        if (self.dim[0] == 3):
            v0 = self.volumes[1,0,0,0,0,0]
            dF_dV = (e[0,0,0,0,0,0]-e[2,0,0,0,0,0])/(2*dv)
            d2F_dV2 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+e[2,0,0,0,0,0])/(dv)**2
            dfdv = dF_dV + (v-v0)*d2F_dV2
        elif (self.dim[0] == 5):
            v0 = self.volumes[2,0,0,0,0,0]
            dF_dV = (-e[0,0,0,0,0,0]+ 8*e[1,0,0,0,0,0]-8*e[3,0,0,0,0,0]+e[4,0,0,0,0,0])/(12*dv)
            d2F_dV2 = (-e[0,0,0,0,0,0]+16*e[1,0,0,0,0,0]-30*e[2,0,0,0,0,0]+16*e[3,0,0,0,0,0]-e[4,0,0,0,0,0])/(12*dv**2)
            d3F_dV3 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+2*e[3,0,0,0,0,0]-e[4,0,0,0,0,0])/(2*dv**3)
            d4F_dV4 = (e[0,0,0,0,0,0]-4*e[1,0,0,0,0,0]+6*e[2,0,0,0,0,0]-4*e[3,0,0,0,0,0]+e[4,0,0,0,0,0])/(dv**4)
            dfdv = dF_dV + (v-v0)*d2F_dV2+ 0.5*(v-v0)**2*d3F_dV3+ 1/6.0*(v-v0)**3*d4F_dV4

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_a = -dfdv* self.eVA3_HaBohr3

        stress[0] = stress_a -pressure
        stress[1] = stress_a -pressure
        stress[2] = stress_a -pressure

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        return  dtol, stress

    def stress_ZSISA_1DOF(self, temp, pressure) -> tuple:

        e,S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,0,0,0,0,0]
        X1 = self.ave_x[1,1,1,0,0,0]
        XBO = self.ave_x_bo

        dexx = (X0-X1)/XBO
        exx0 = X1/XBO-1

        dF_dX = (e[0,0,0,0,0,0]-e[2,2,2,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,0,0,0,0,0]-2*e[1,1,1,0,0,0]+e[2,2,2,0,0,0])/(dexx)**2

        dF_dV = (e[0,0,0,0,0,0]-e[2,2,2,0,0,0])/(2*(X0-X1)) /(3*X1**2)
        d2F_dV2 = (e[0,0,0,0,0,0]-2*e[1,1,1,0,0,0]+e[2,2,2,0,0,0])/(X0-X1)**2  /(3*X1**2)**2

        dS_dX = (S[0,0,0,0,0,0]-S[2,2,2,0,0,0])/(2*dexx)
        d2S_dX2 = (S[0,0,0,0,0,0]-2*S[1,1,1,0,0,0]+S[2,2,2,0,0,0])/(dexx)**2

        x = self.ave_x_guess
        v = self.volume_guess
        exx_n = x/XBO-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2
        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_xx = -dfdx/v*(exx_n+1)/3.0 * self.eVA3_HaBohr3
        print (x/XBO, x/XBO, x/XBO)
        print (stress_xx , stress_xx,stress_xx )
        stress[0] = stress_xx -pressure
        stress[1] = stress_xx -pressure
        stress[2] = stress_xx -pressure

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        therm = [0,0,0,0,0,0]
        if all(dtol[i] < 1e-8 for i in range(6)):
            if os.path.exists("elastic_constant.txt"):
                matrix_elastic = self.elastic_constants("elastic_constant.txt")
                matrix_elastic = np.array(matrix_elastic)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                M = np.array([[ d2F_dX2/3*(exx0+1)**2  , 0          , 0         ,0,0,0],
                              [    0       , d2F_dX2/3*(exx0+1)**2  , 0         ,0,0,0],
                              [    0       , 0          , d2F_dX2/3*(exx0+1)**2 ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                M = M + matrix_elastic
                S1 = dsdx*(exx_n+1)/3
                S = np.array([S1,S1,S1,0,0,0])
                dstrain_dt = np.linalg.inv(M) @ S
                therm=[dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(exx_n+1),dstrain_dt[2]*(exx_n+1),0,0,0]
                print ("therm")
                print (therm)

        return dtol, stress, therm

    def stress_ZSISA_2DOF(self, temp, pressure) -> tuple:

        e,S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,0,1,0,0,0]
        Z0 = self.ave_z[1,1,0,0,0,0]

        X1 = self.ave_x[1,1,1,0,0,0]
        Z1 = self.ave_z[1,1,1,0,0,0]

        X2 = self.ave_x[2,2,1,0,0,0]
        Z2 = self.ave_z[1,1,2,0,0,0]

        XBO = self.ave_x_bo
        ZBO = self.ave_z_bo

        V = self.volume_guess

        dexx = (X0-X1)/XBO
        dezz = (Z0-Z1)/ZBO

        exx0 = X1/XBO-1
        ezz0 = Z1/ZBO-1

        dF_dX = (e[0,0,1,0,0,0]-e[2,2,1,0,0,0])/(2*dexx)
        dF_dZ = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dezz)

        d2F_dX2 = (e[0,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,2,1,0,0,0])/(dexx)**2
        d2F_dZ2 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dezz)**2
        d2F_dXdZ = (e[1,1,1,0,0,0] - e[0,0,1,0,0,0] - e[1,1,0,0,0,0] + e[0,0,0,0,0,0]) / (dexx *dezz)

        dS_dX = (S[0,0,1,0,0,0]-S[2,2,1,0,0,0])/(2*dexx)
        dS_dZ = (S[1,1,0,0,0,0]-S[1,1,2,0,0,0])/(2*dezz)
        d2S_dX2 = (S[0,0,1,0,0,0]-2*S[1,1,1,0,0,0]+S[2,2,1,0,0,0])/(dexx)**2
        d2S_dZ2 = (S[1,1,0,0,0,0]-2*S[1,1,1,0,0,0]+S[1,1,2,0,0,0])/(dezz)**2
        d2S_dXdZ = (S[1,1,1,0,0,0] - S[0,0,1,0,0,0] - S[1,1,0,0,0,0] + S[0,0,0,0,0,0]) / (dexx *dezz)

        x = self.ave_x_guess
        z = self.ave_z_guess
        v = self.volume_guess

        exx_n = x/XBO-1
        ezz_n = z/ZBO-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(ezz_n-ezz0)*d2F_dXdZ
        dfdz = dF_dZ + (ezz_n-ezz0)*d2F_dZ2+(exx_n-exx0)*d2F_dXdZ

        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2+(ezz_n-ezz0)*d2S_dXdZ
        dsdz = dS_dZ + (ezz_n-ezz0)*d2S_dZ2+(exx_n-exx0)*d2S_dXdZ

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_xx = -dfdx/V*(exx_n+1)*0.5 * self.eVA3_HaBohr3
        stress_zz = -dfdz/V*(ezz_n+1)     * self.eVA3_HaBohr3

        stress[0] = stress_xx -pressure
        stress[1] = stress_xx -pressure
        stress[2] = stress_zz -pressure
        print (x/XBO, x/XBO, z/ZBO)
        print ("stress",pressure)
        print (stress[0],stress[2])

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        therm = [0,0,0,0,0,0]
        if all(dtol[i] < 1e-8 for i in range(6)):
            if os.path.exists("elastic_constant.txt"):
                matrix_elastic = self.elastic_constants("elastic_constant.txt")
                matrix_elastic = np.array(matrix_elastic)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                M = np.array([[ d2F_dX2/2 *(exx0+1)**2       , 0                            , d2F_dXdZ/2*(exx0+1)*(ezz0+1) ,0,0,0],
                              [    0                         , d2F_dX2/2 *(exx0+1)**2       , d2F_dXdZ/2*(exx0+1)*(ezz0+1) ,0,0,0],
                              [ d2F_dXdZ/2*(exx0+1)*(ezz0+1) , d2F_dXdZ/2*(exx0+1)*(ezz0+1) , d2F_dZ2*(ezz0+1)**2          ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                M = M + matrix_elastic+P/self.eVA3_HaBohr3
                S1 = dsdx*(exx_n+1)*0.5
                S3 = dsdz*(ezz_n+1)
                S = np.array([S1,S1,S3,0,0,0])
                dstrain_dt = np.linalg.inv(M) @ S
                therm = [dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(exx_n+1),dstrain_dt[2]*(ezz_n+1),0,0,0]
                print (therm)


        return dtol, stress, therm

    def stress_ZSISA_3DOF(self, temp, pressure, mode) -> tuple:

        e,S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,1,1,0,0,0]
        Y0 = self.ave_y[1,0,1,0,0,0]
        Z0 = self.ave_z[1,1,0,0,0,0]

        X1 = self.ave_x[1,1,1,0,0,0]
        Y1 = self.ave_y[1,1,1,0,0,0]
        Z1 = self.ave_z[1,1,1,0,0,0]

        XBO = self.ave_x_bo
        YBO = self.ave_y_bo
        ZBO = self.ave_z_bo

        dexx = (X0-X1)/XBO
        exx0 = X1/XBO-1
        if (self.sym == "cubic"):
            deyy = dezz = dexx
            eyy0 = ezz0 = exx0
        elif (self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal"):
            deyy = dexx
            eyy0 = exx0
            dezz = (Z0-Z1)/ZBO
            ezz0 = Z1/ZBO-1
        else :
            dexx = (X0-X1)/XBO
            dezz = (Z0-Z1)/ZBO
            deyy = (Y0-y1)/YBO

            exx0 = X1/XBO-1
            eyy0 = Y1/YBO-1
            ezz0 = Z1/ZBO-1



        dF_dX = (e[0,1,1,0,0,0]-e[2,1,1,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,1,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,1,1,0,0,0])/(dexx)**2
        dS_dX = (S[0,1,1,0,0,0]-S[2,1,1,0,0,0])/(2*dexx)
        d2S_dX2 = (S[0,1,1,0,0,0]-2*S[1,1,1,0,0,0]+S[2,1,1,0,0,0])/(dexx)**2
        if (self.sym == "cubic" or self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal"):
            d2F_dXdY = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[0,1,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
            d2S_dXdY = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[0,1,1,0,0,0] + S[0,0,1,0,0,0]) / (dexx *deyy)
        if (self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal" or self.sym == "orthorhombic"):
            dF_dZ = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dezz)
            d2F_dZ2 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dezz)**2
            d2F_dXdZ = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,1,0,0,0,0] + e[0,1,0,0,0,0]) / (dexx *dezz)
            dS_dZ = (S[1,1,0,0,0,0]-S[1,1,2,0,0,0])/(2*dezz)
            d2S_dZ2 = (S[1,1,0,0,0,0]-2*S[1,1,1,0,0,0]+S[1,1,2,0,0,0])/(dezz)**2
            d2S_dXdZ = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[1,1,0,0,0,0] + S[0,1,0,0,0,0]) / (dexx *dezz)


        if (self.sym == "cubic"):
            dF_dY = dF_dX
            dF_dZ = dF_dX
            d2F_dY2 = d2F_dX2
            d2F_dZ2 = d2F_dX2
            d2F_dXdZ = d2F_dXdY
            d2F_dYdZ = d2F_dXdY
            dS_dY = dS_dX
            dS_dZ = dS_dX
            d2S_dY2 = d2S_dX2
            d2S_dZ2 = d2S_dX2
            d2S_dXdZ = d2S_dXdY
            d2S_dYdZ = d2S_dXdY
        elif (self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal"):
            dF_dY = dF_dX
            d2F_dY2 = d2F_dX2
            d2F_dYdZ = d2F_dXdZ
            dS_dY = dS_dX
            d2S_dY2 = d2S_dX2
            d2S_dYdZ = d2S_dXdZ
        elif (self.sym == "orthorhombic"):
            dF_dY = (e[1,0,1,0,0,0]-e[1,2,1,0,0,0])/(2*deyy)
            d2F_dY2 = (e[1,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[1,2,1,0,0,0])/(deyy)**2
            d2F_dXdY = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,0,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
            d2F_dYdZ = (e[1,1,1,0,0,0] - e[1,0,1,0,0,0] - e[1,1,0,0,0,0] + e[1,0,0,0,0,0]) / (deyy *dezz)

            dS_dY = (S[1,0,1,0,0,0]-S[1,2,1,0,0,0])/(2*deyy)
            d2S_dY2 = (S[1,0,1,0,0,0]-2*S[1,1,1,0,0,0]+S[1,2,1,0,0,0])/(deyy)**2
            d2S_dXdY = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[1,0,1,0,0,0] + S[0,0,1,0,0,0]) / (dexx *deyy)
            d2S_dYdZ = (S[1,1,1,0,0,0] - S[1,0,1,0,0,0] - S[1,1,0,0,0,0] + S[1,0,0,0,0,0]) / (deyy *dezz)
        d2F_dXY2  = 0.0
        d2F_dXdYZ = 0.0
        d2F_dYZ2  = 0.0
        d2F_dXZ2  = 0.0
        if mode == 'ECs':
            deyz= (self.ave_y[1,1,1,0,0,0] - self.ave_y[1,1,1,1,0,0])/ZBO
            #deyz= (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,1,0,0])/ZBO
            Fyz=[e[1,1,1,2,0,0],e[1,1,1,1,0,0],e[1,1,1,1,0,0],e[1,1,1,2,0,0]]
            Dyz=[2*deyz,deyz,-deyz,-2*deyz]
            param = np.polyfit(Dyz,Fyz,2)
            d2F_dYZ2 = 2*param[0]
            if self.sym == "trigonal":
                d2F_dXdYZ= (e[1,1,1,1,0,0] - e[0,1,1,1,0,0] - e[1,1,1,0,0,0] + e[0,1,1,0,0,0]) / (dexx *deyy)
            if self.sym == "tetragonal" or self.sym == "orthorhombic":
                dexy= (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,0,0,1])/YBO
                Fxy=[e[1,1,1,0,0,2],e[1,1,1,0,0,1],e[1,1,1,0,0,1],e[1,1,1,0,0,2]]
                Dxy=[2*dexy,dexy,-dexy,-2*dexy]
                param = np.polyfit(Dxy,Fxy,2)
                d2F_dXY2 = 2*param[0]
            if  self.sym == "orthorhombic":
                dexz= (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,0,1,0])/ZBO
                Fxz=[e[1,1,1,0,2,0],e[1,1,1,0,1,0],e[1,1,1,0,1,0],e[1,1,1,0,2,0]]
                Dxz=[2*dexz,dexz,-dexz,-2*dexz]
                param = np.polyfit(Dxz,Fxz,2)
                d2F_dXZ2 = 2*param[0]


        x = self.ave_x_guess
        y = self.ave_y_guess
        z = self.ave_z_guess
        v = self.volume_guess

        scale_x_guess = x/XBO
        scale_y_guess = y/YBO
        scale_z_guess = z/ZBO

        exx_n = scale_x_guess-1
        eyy_n = scale_y_guess-1
        ezz_n = scale_z_guess-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(eyy_n-eyy0)*d2F_dXdY+(ezz_n-ezz0)*d2F_dXdZ
        dfdy = dF_dY + (eyy_n-eyy0)*d2F_dY2+(exx_n-exx0)*d2F_dXdY+(ezz_n-ezz0)*d2F_dYdZ
        dfdz = dF_dZ + (ezz_n-ezz0)*d2F_dZ2+(exx_n-exx0)*d2F_dXdZ+(eyy_n-eyy0)*d2F_dYdZ

        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2+(eyy_n-eyy0)*d2S_dXdY+(ezz_n-ezz0)*d2S_dXdZ
        dsdy = dS_dY + (eyy_n-eyy0)*d2S_dY2+(exx_n-exx0)*d2S_dXdY+(ezz_n-ezz0)*d2S_dYdZ
        dsdz = dS_dZ + (ezz_n-ezz0)*d2S_dZ2+(exx_n-exx0)*d2S_dXdZ+(eyy_n-eyy0)*d2S_dYdZ

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_xx = -dfdx/v*(exx_n+1)* self.eVA3_HaBohr3
        stress_yy = -dfdy/v*(eyy_n+1)* self.eVA3_HaBohr3
        stress_zz = -dfdz/v*(ezz_n+1)* self.eVA3_HaBohr3

        if os.path.exists("elastic_constant.txt"):
            matrix_elastic = self.elastic_constants("elastic_constant.txt")
            matrix_elastic = np.array(matrix_elastic)
            matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
            scale_xx=(exx0+1)*(exx0+1)
            scale_yy=(eyy0+1)*(eyy0+1)
            scale_zz=(ezz0+1)*(ezz0+1)
            scale_xz=(exx0+1)*(ezz0+1)
            scale_yz=(eyy0+1)*(ezz0+1)
            scale_xy=(exx0+1)*(eyy0+1)
            M = np.array([[ d2F_dX2*scale_xx, d2F_dXdY*scale_xy ,d2F_dXdZ*scale_xz ,d2F_dXdYZ*scale_xz  ,0.0               ,0.0              ],
                          [d2F_dXdY*scale_xy, d2F_dY2 *scale_yy ,d2F_dYdZ*scale_yz ,0.0                 ,0.0               ,0.0              ],
                          [d2F_dXdZ*scale_xz, d2F_dYdZ*scale_yz ,d2F_dZ2 *scale_zz ,0.0                 ,0.0               ,0.0              ],
                          [d2F_dXdYZ*scale_xz,0.0               ,0.0               ,d2F_dYZ2*scale_zz   ,0.0               ,0.0              ],
                          [0.0              , 0.0               ,0.0               ,0.0                 ,d2F_dXZ2*scale_zz ,0.0              ],
                          [0.0              , 0.0               ,0.0               ,0.0                 ,0.0               ,d2F_dXY2*scale_yy]])
            P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                          [pressure*v ,0.0,pressure*v  ,0,0,0],
                          [pressure*v ,pressure*v,0.0  ,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0]])
            M = M + matrix_elastic+P/ self.eVA3_HaBohr3
            S1 = dsdx*(exx_n+1)
            S2 = dsdy*(eyy_n+1)
            S3 = dsdz*(ezz_n+1)
            S = np.array([S1,S2,S3,0,0,0])
            dstrain_dt = np.linalg.inv(M) @ S
            therm=[dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),0,0,0]
            print ("therm")
            print (therm)
            M=M/v*abu.eVA3_GPa

            with open("elastic.txt", "w") as f:
                f.write("Elastic [GPa] \n")
                if mode == 'ECs':
                    if self.sym == "cubic":
                        f.write(f"    {'C11':12s} {'C12':12s} {'C44':12s}\n")
                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[3,3]:12.6f} \n")
                    elif self.sym == "hexagonal":
                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C44':12s}\n")
                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} \n")
                    elif self.sym == "trigonal":
                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C14':12s} {'C44':12s}\n")
                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[0,3]:12.6f} {M[3,3]:12.6f} \n")
                    elif self.sym == "tetragonal":
                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C44':12s} {'C66':12s}\n")
                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} {M[5,5]:12.6f} \n")
                    if  self.sym == "orthorhombic":
                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C22':12s} {'C23':12s} {'C33':12s} {'C44':12s} {'C55':12s} {'C66':12s}\n")
                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[1,1]:12.6f} {M[1,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} {M[4,4]:12.6f} {M[5,5]:12.6f} \n")
                else:
                    f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C22':12s} {'C23':12s} {'C33':12s} \n")
                    f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[1,1]:12.6f} {M[1,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f}  \n")
        else:
            therm = [0,0,0,0,0,0]

        stress[0] = stress_xx -pressure
        stress[1] = stress_yy -pressure
        stress[2] = stress_zz -pressure
        print (scale_x_guess,scale_y_guess,scale_z_guess)
        print ("stress")
        print (stress[0],stress[1],stress[2])

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        return dtol, stress, therm

    def stress_ZSISA_monoclinic(self, temp, pressure, mode) -> tuple:

        e,S = self.get_vib_free_energies(temp)
        AxBO =self.ax_bo
        ByBO =self.by_bo
        CzBO =self.cz_bo
        BxBO =self.bx_bo
        CxBO =self.cx_bo
        CyBO =self.cy_bo

        Ax1 = self.ax[1,1,1,1,1,1]
        Bx1 = self.bx[1,1,1,1,1,1]
        By1 = self.by[1,1,1,1,1,1]
        Cx1 = self.cx[1,1,1,1,1,1]
        Cy1 = self.cy[1,1,1,1,1,1]
        Cz1 = self.cz[1,1,1,1,1,1]

        Ax0 =self.ax[0,1,1,1,1,1]
        By0 =self.by[1,0,1,1,1,1]
        Cz0 =self.cz[1,1,0,1,1,1]
        Bx0 =self.bx[1,1,1,1,1,0]
        Cx0 =self.cx[0,1,1,1,0,1]
        Cy0 =self.cy[1,1,1,0,1,1]

        dexx= (Ax0-Ax1)/AxBO
        deyy= (By0-By1)/ByBO
        dezz= (Cz0-Cz1)/CzBO
        dexz= (AxBO*(Cx0-Cx1)-(Ax0-Ax1)*CxBO)/(AxBO*CzBO)

        exx0= Ax1/AxBO-1
        eyy0= By1/ByBO-1
        ezz0= Cz1/CzBO-1
        exz0= (AxBO*Cx1-Ax1*CxBO)/(AxBO*CzBO)


        dF_dA1   = (e[0,1,1,1,1,1]-e[2,1,1,1,1,1])/(2*dexx)
        dF_dB2   = (e[1,0,1,1,1,1]-e[1,2,1,1,1,1])/(2*deyy)
        dF_dC3   = (e[1,1,0,1,1,1]-e[1,1,2,1,1,1])/(2*dezz)
        dF_dC1   = (e[1,1,1,1,0,1]-e[1,1,1,1,2,1])/(2*dexz)


        d2F_dA12 = (e[0,1,1,1,1,1]-2*e[1,1,1,1,1,1]+e[2,1,1,1,1,1])/(dexx)**2
        d2F_dB22 = (e[1,0,1,1,1,1]-2*e[1,1,1,1,1,1]+e[1,2,1,1,1,1])/(deyy)**2
        d2F_dC32 = (e[1,1,0,1,1,1]-2*e[1,1,1,1,1,1]+e[1,1,2,1,1,1])/(dezz)**2
        d2F_dC12 = (e[1,1,1,1,0,1]-2*e[1,1,1,1,1,1]+e[1,1,1,1,2,1])/(dexz)**2

        d2F_dA1dB2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,0,1,1,1,1] + e[0,0,1,1,1,1]) / (dexx *deyy)
        d2F_dA1dC3 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,0,1,1,1] + e[0,1,0,1,1,1]) / (dexx *dezz)
        d2F_dA1dC1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,0,1] + e[0,1,1,1,0,1]) / (dexx *dexz)

        d2F_dB2dC3 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,0,1,1,1] + e[1,0,0,1,1,1]) / (deyy *dezz)
        d2F_dB2dC1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,0,1] + e[1,0,1,1,0,1]) / (deyy *dexz)

        d2F_dC3dC1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,0,1] + e[1,1,0,1,0,1]) / (dezz *dexz)


        dS_dA1   = (S[0,1,1,1,1,1]-S[2,1,1,1,1,1])/(2*dexx)
        dS_dB2   = (S[1,0,1,1,1,1]-S[1,2,1,1,1,1])/(2*deyy)
        dS_dC3   = (S[1,1,0,1,1,1]-S[1,1,2,1,1,1])/(2*dezz)
        dS_dC1   = (S[1,1,1,1,0,1]-S[1,1,1,1,2,1])/(2*dexz)

        d2S_dA12 = (S[0,1,1,1,1,1]-2*S[1,1,1,1,1,1]+S[2,1,1,1,1,1])/(dexx)**2
        d2S_dB22 = (S[1,0,1,1,1,1]-2*S[1,1,1,1,1,1]+S[1,2,1,1,1,1])/(deyy)**2
        d2S_dC32 = (S[1,1,0,1,1,1]-2*S[1,1,1,1,1,1]+S[1,1,2,1,1,1])/(dezz)**2
        d2S_dC12 = (S[1,1,1,1,0,1]-2*S[1,1,1,1,1,1]+S[1,1,1,1,2,1])/(dexz)**2

        d2S_dA1dB2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,0,1,1,1,1] + S[0,0,1,1,1,1]) / (dexx *deyy)
        d2S_dA1dC3 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,0,1,1,1] + S[0,1,0,1,1,1]) / (dexx *dezz)
        d2S_dA1dC1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,0,1] + S[0,1,1,1,0,1]) / (dexx *dexz)

        d2S_dB2dC3 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,0,1,1,1] + S[1,0,0,1,1,1]) / (deyy *dezz)
        d2S_dB2dC1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,0,1] + S[1,0,1,1,0,1]) / (deyy *dexz)

        d2S_dC3dC1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,0,1] + S[1,1,0,1,0,1]) / (dezz *dexz)

        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess
        v = self.volume_guess

        ax= a
        by= b
        cz= c*math.sin(math.pi*self.angles_guess[1]/180)
        cx= c*math.cos(math.pi*self.angles_guess[1]/180)

        exx_n= ax/AxBO-1
        eyy_n= by/ByBO-1
        ezz_n= cz/CzBO-1
        exz_n= (AxBO*cx-ax*CxBO)/(AxBO*CzBO)

        dfda1= dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dA1dC3+(exz_n-exz0)*d2F_dA1dC1
        dfdb2= dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dB2dC3+(exz_n-exz0)*d2F_dB2dC1
        dfdc3= dF_dC3 + (ezz_n-ezz0)*d2F_dC32+(exx_n-exx0)*d2F_dA1dC3+(eyy_n-eyy0)*d2F_dB2dC3+(exz_n-exz0)*d2F_dC3dC1
        dfdc1= dF_dC1 + (exz_n-exz0)*d2F_dC12+(exx_n-exx0)*d2F_dA1dC1+(eyy_n-eyy0)*d2F_dB2dC1+(ezz_n-ezz0)*d2F_dC3dC1

        dsda1= dS_dA1 + (exx_n-exx0)*d2S_dA12+(eyy_n-eyy0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dA1dC3+(exz_n-exz0)*d2S_dA1dC1
        dsdb2= dS_dB2 + (eyy_n-eyy0)*d2S_dB22+(exx_n-exx0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dB2dC3+(exz_n-exz0)*d2S_dB2dC1
        dsdc3= dS_dC3 + (ezz_n-ezz0)*d2S_dC32+(exx_n-exx0)*d2S_dA1dC3+(eyy_n-eyy0)*d2S_dB2dC3+(exz_n-exz0)*d2S_dC3dC1
        dsdc1= dS_dC1 + (exz_n-exz0)*d2S_dC12+(exx_n-exx0)*d2S_dA1dC1+(eyy_n-eyy0)*d2S_dB2dC1+(ezz_n-ezz0)*d2S_dC3dC1

        dtol   =np.zeros(6)
        stress =np.zeros(6)

        stress_a1= -dfda1/v*(exx_n+1)* self.eVA3_HaBohr3
        stress_b2= -dfdb2/v*(eyy_n+1)* self.eVA3_HaBohr3
        stress_c3= -dfdc3/v*(ezz_n+1)* self.eVA3_HaBohr3
        stress_c1= -1.0/v*(dfdc1*(ezz_n+1)+dfda1*exz_n) * self.eVA3_HaBohr3
        print (ax/AxBO, by/ByBO, cz/CzBO)
        print (stress_a1,stress_b2,stress_c3,stress_c1)

        if os.path.exists("elastic_constant.txt"):
            if mode == 'ECs':
                dexy= (Bx0-Bx1)/(ByBO)
                deyz= (Cy0-Cy1)/(CzBO)
                d2F_dB12 = (e[1,1,1,1,1,0]-2*e[1,1,1,1,1,1]+e[1,1,1,1,1,0])/(deyz)**2
                d2F_dC22 = (e[1,1,1,0,1,1]-2*e[1,1,1,1,1,1]+e[1,1,1,0,1,1])/(dexy)**2
                d2F_dB1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,1,0] + e[1,1,1,0,1,0]) / (dexy *deyz)
            else:
                d2F_dB12 =  0.0
                d2F_dC22 =  0.0
                d2F_dB1dC2= 0.0
            matrix_elastic=self.elastic_constants("elastic_constant.txt")
            matrix_elastic = np.array(matrix_elastic)
            matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
            scale_xx=(exx0+1)*(exx0+1)
            scale_yy=(eyy0+1)*(eyy0+1)
            scale_zz=(ezz0+1)*(ezz0+1)
            scale_xz=(exx0+1)*(ezz0+1)
            scale_yz=(eyy0+1)*(ezz0+1)
            scale_xy=(exx0+1)*(eyy0+1)
            M = np.array([[ d2F_dA12 *scale_xx, d2F_dA1dB2*scale_xy ,d2F_dA1dC3*scale_xz ,0.0                 ,d2F_dA1dC1*scale_xz ,0.0                ],
                          [d2F_dA1dB2*scale_xy, d2F_dB22  *scale_yy ,d2F_dB2dC3*scale_yz ,0.0                 ,d2F_dB2dC1*scale_yz ,0.0                ],
                          [d2F_dA1dC3*scale_xz, d2F_dB2dC3*scale_yz ,d2F_dC32  *scale_zz ,0.0                 ,d2F_dC3dC1*scale_zz ,0.0                ],
                          [0.0                , 0.0                 ,0.0                 ,d2F_dC22  *scale_zz ,0.0                 ,d2F_dB1dC2*scale_yz],
                          [d2F_dA1dC1*scale_xz, d2F_dB2dC1*scale_yz ,d2F_dC3dC1*scale_zz ,0.0                 ,d2F_dC12  *scale_zz ,0.0                ],
                          [0.0                , 0.0                 ,0.0                 ,d2F_dB1dC2*scale_yz ,0.0                 ,d2F_dB12  *scale_yy]])
            P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                          [pressure*v ,0.0,pressure*v  ,0,0,0],
                          [pressure*v ,pressure*v,0.0  ,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0]])
            M = M + matrix_elastic+P/ self.eVA3_HaBohr3
            S1= dsda1*(exx_n+1)
            S2= dsdb2*(eyy_n+1)
            S3= dsdc3*(ezz_n+1)
            S4= 0.0
            S5= dsdc1*(ezz_n+1)+dsda1*exz_n
            S6= 0.0
            S = np.array([S1,S2,S3,S4,S5,S6])
            dstrain_dt = np.linalg.inv(M) @ S
            #therm=[dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),0,dstrain_dt[4]*(ezz_n+1)+dstrain_dt[0]*exz_n,0]
            therm=[dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1)-(dstrain_dt[4]*(ezz_n+1)+dstrain_dt[0]*exz_n)/cz,0,-dstrain_dt[4],0]

            M=M/v*abu.eVA3_GPa
            with open("elastic.txt", "w") as f:
                f.write("Elastic [GPa] \n")
                if mode != 'ECs':
                    f.write(f" Warning: C44, C46, and C66 do not include the free energy contribution (only BO energy).\n")
                f.write(f" \t   xx\t\tyy\t\tzz\t\tyz\t\txz\t\txy\n")
                f.write(f" xx {M[0,0]:14.8f}  {M[0,1]:14.8f}  {M[0,2]:14.8f}  {M[0,3]:14.8f}  {M[0,4]:14.8f}  {M[0,5]:14.8f}\n")
                f.write(f" yy {M[1,0]:14.8f}  {M[1,1]:14.8f}  {M[1,2]:14.8f}  {M[1,3]:14.8f}  {M[1,4]:14.8f}  {M[1,5]:14.8f}\n")
                f.write(f" zz {M[2,0]:14.8f}  {M[2,1]:14.8f}  {M[2,2]:14.8f}  {M[2,3]:14.8f}  {M[2,4]:14.8f}  {M[2,5]:14.8f}\n")
                f.write(f" yz {M[3,0]:14.8f}  {M[3,1]:14.8f}  {M[3,2]:14.8f}  {M[3,3]:14.8f}  {M[3,4]:14.8f}  {M[3,5]:14.8f}\n")
                f.write(f" xz {M[4,0]:14.8f}  {M[4,1]:14.8f}  {M[4,2]:14.8f}  {M[4,3]:14.8f}  {M[4,4]:14.8f}  {M[4,5]:14.8f}\n")
                f.write(f" xy {M[5,0]:14.8f}  {M[5,1]:14.8f}  {M[5,2]:14.8f}  {M[5,3]:14.8f}  {M[5,4]:14.8f}  {M[5,5]:14.8f}\n")
            print ("therm")
            print (therm)
        else:
            therm = [0,0,0,0,0,0]

        stress[0] = stress_a1 -pressure
        stress[1] = stress_b2 -pressure
        stress[2] = stress_c3 -pressure
        stress[4] = stress_c1

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])
        dtol[4] = abs(stress[4]-self.stress_guess[2,0])

        return dtol, stress , therm

    def stress_ZSISA_triclinic(self, temp, pressure, mode) -> tuple:

        e,S = self.get_vib_free_energies(temp)

        AxBO =self.ax_bo
        ByBO =self.by_bo
        CzBO =self.cz_bo
        BxBO =self.bx_bo
        CxBO =self.cx_bo
        CyBO =self.cy_bo

        Ax0 = self.ax[0,1,1,1,1,1]
        Bx0 = self.bx[0,1,1,1,1,0]
        By0 = self.by[1,0,1,1,1,1]
        Cx0 = self.cx[0,1,1,1,1,0]+0.5*(self.cx[1,1,1,1,0,1]-self.cx[1,1,1,1,2,1])
        Cy0 = self.cy[1,0,1,0,1,1]
        Cz0 = self.cz[1,1,0,1,1,1]

        Ax1 = self.ax[1,1,1,1,1,1]
        Bx1 = self.bx[1,1,1,1,1,1]
        By1 = self.by[1,1,1,1,1,1]
        Cx1 = self.cx[1,1,1,1,1,1]
        Cy1 = self.cy[1,1,1,1,1,1]
        Cz1 = self.cz[1,1,1,1,1,1]

        V = self.volume_guess

        dexx = (Ax0-Ax1)/AxBO
        deyy = (By0-By1)/ByBO
        dezz = (Cz0-Cz1)/CzBO
        dexy = (AxBO*(Bx0-Bx1)-(Ax0-Ax1)*BxBO)/(AxBO*ByBO)
        deyz = (ByBO*(Cy0-Cy1)-(By0-By1)*CyBO)/(ByBO*CzBO)
        dexz = (AxBO*(ByBO*(Cx0-Cx1)-(Bx0-Bx1)*CyBO)-(Ax0-Ax1)*(ByBO*CxBO-BxBO*CyBO))/(AxBO*ByBO*CzBO)

        exx0 = Ax1/AxBO-1
        eyy0 = By1/ByBO-1
        ezz0 = Cz1/CzBO-1
        exy0 = (Ax0*Bx1-Ax1*Bx0)/(AxBO*ByBO)
        eyz0 = (By0*Cy1-By1*Cy0)/(ByBO*CzBO)
        exz0 = (AxBO*(ByBO*Cx1-Bx1*CyBO)-Ax1*(ByBO*CxBO-BxBO*CyBO))/(AxBO*ByBO*CzBO)

        dF_dA1 = (e[0,1,1,1,1,1]-e[2,1,1,1,1,1])/(2*dexx)
        dF_dB2 = (e[1,0,1,1,1,1]-e[1,2,1,1,1,1])/(2*deyy)
        dF_dC3 = (e[1,1,0,1,1,1]-e[1,1,2,1,1,1])/(2*dezz)
        dF_dC2 = (e[1,1,1,0,1,1]-e[1,1,1,2,1,1])/(2*deyz)
        dF_dC1 = (e[1,1,1,1,0,1]-e[1,1,1,1,2,1])/(2*dexz)
        dF_dB1 = (e[1,1,1,1,1,0]-e[1,1,1,1,1,2])/(2*dexy)

        d2F_dA12 = (e[0,1,1,1,1,1]-2*e[1,1,1,1,1,1]+e[2,1,1,1,1,1])/(dexx)**2
        d2F_dB22 = (e[1,0,1,1,1,1]-2*e[1,1,1,1,1,1]+e[1,2,1,1,1,1])/(deyy)**2
        d2F_dC32 = (e[1,1,0,1,1,1]-2*e[1,1,1,1,1,1]+e[1,1,2,1,1,1])/(dezz)**2
        d2F_dC22 = (e[1,1,1,0,1,1]-2*e[1,1,1,1,1,1]+e[1,1,1,2,1,1])/(deyz)**2
        d2F_dC12 = (e[1,1,1,1,0,1]-2*e[1,1,1,1,1,1]+e[1,1,1,1,2,1])/(dexz)**2
        d2F_dB12 = (e[1,1,1,1,1,0]-2*e[1,1,1,1,1,1]+e[1,1,1,1,1,2])/(dexy)**2

        d2F_dA1dB2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,0,1,1,1,1] + e[0,0,1,1,1,1]) / (dexx *deyy)
        d2F_dA1dC3 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,0,1,1,1] + e[0,1,0,1,1,1]) / (dexx *dezz)
        d2F_dA1dC2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,0,1,1] + e[0,1,1,0,1,1]) / (dexx *deyz)
        d2F_dA1dC1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,0,1] + e[0,1,1,1,0,1]) / (dexx *dexz)
        d2F_dA1dB1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,1,0] + e[0,1,1,1,1,0]) / (dexx *dexy)

        d2F_dB2dC3 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,0,1,1,1] + e[1,0,0,1,1,1]) / (deyy *dezz)
        d2F_dB2dC2 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,0,1,1] + e[1,0,1,0,1,1]) / (deyy *deyz)
        d2F_dB2dC1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,0,1] + e[1,0,1,1,0,1]) / (deyy *dexz)
        d2F_dB2dB1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,1,0] + e[1,0,1,1,1,0]) / (deyy *dexy)

        d2F_dC3dC2 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,0,1,1] + e[1,1,0,0,1,1]) / (dezz *deyz)
        d2F_dC3dC1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,0,1] + e[1,1,0,1,0,1]) / (dezz *dexz)
        d2F_dC3dB1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,1,0] + e[1,1,0,1,1,0]) / (dezz *dexy)

        d2F_dC1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,0,1] + e[1,1,1,0,0,1]) / (deyz *dexz)
        d2F_dB1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,1,0] + e[1,1,1,0,1,0]) / (deyz *dexy)

        d2F_dB1dC1 = (e[1,1,1,1,1,1] - e[1,1,1,1,0,1] - e[1,1,1,1,1,0] + e[1,1,1,1,0,0]) / (dexz *dexy)

        dS_dA1 = (S[0,1,1,1,1,1]-S[2,1,1,1,1,1])/(2*dexx)
        dS_dB2 = (S[1,0,1,1,1,1]-S[1,2,1,1,1,1])/(2*deyy)
        dS_dC3 = (S[1,1,0,1,1,1]-S[1,1,2,1,1,1])/(2*dezz)
        dS_dC2 = (S[1,1,1,0,1,1]-S[1,1,1,2,1,1])/(2*deyz)
        dS_dC1 = (S[1,1,1,1,0,1]-S[1,1,1,1,2,1])/(2*dexz)
        dS_dB1 = (S[1,1,1,1,1,0]-S[1,1,1,1,1,2])/(2*dexy)

        d2S_dA12 = (S[0,1,1,1,1,1]-2*S[1,1,1,1,1,1]+S[2,1,1,1,1,1])/(dexx)**2
        d2S_dB22 = (S[1,0,1,1,1,1]-2*S[1,1,1,1,1,1]+S[1,2,1,1,1,1])/(deyy)**2
        d2S_dC32 = (S[1,1,0,1,1,1]-2*S[1,1,1,1,1,1]+S[1,1,2,1,1,1])/(dezz)**2
        d2S_dC22 = (S[1,1,1,0,1,1]-2*S[1,1,1,1,1,1]+S[1,1,1,2,1,1])/(deyz)**2
        d2S_dC12 = (S[1,1,1,1,0,1]-2*S[1,1,1,1,1,1]+S[1,1,1,1,2,1])/(dexz)**2
        d2S_dB12 = (S[1,1,1,1,1,0]-2*S[1,1,1,1,1,1]+S[1,1,1,1,1,2])/(dexy)**2

        d2S_dA1dB2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,0,1,1,1,1] + S[0,0,1,1,1,1]) / (dexx *deyy)
        d2S_dA1dC3 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,0,1,1,1] + S[0,1,0,1,1,1]) / (dexx *dezz)
        d2S_dA1dC2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,0,1,1] + S[0,1,1,0,1,1]) / (dexx *deyz)
        d2S_dA1dC1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,0,1] + S[0,1,1,1,0,1]) / (dexx *dexz)
        d2S_dA1dB1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,1,0] + S[0,1,1,1,1,0]) / (dexx *dexy)

        d2S_dB2dC3 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,0,1,1,1] + S[1,0,0,1,1,1]) / (deyy *dezz)
        d2S_dB2dC2 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,0,1,1] + S[1,0,1,0,1,1]) / (deyy *deyz)
        d2S_dB2dC1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,0,1] + S[1,0,1,1,0,1]) / (deyy *dexz)
        d2S_dB2dB1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,1,0] + S[1,0,1,1,1,0]) / (deyy *dexy)

        d2S_dC3dC2 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,0,1,1] + S[1,1,0,0,1,1]) / (dezz *deyz)
        d2S_dC3dC1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,0,1] + S[1,1,0,1,0,1]) / (dezz *dexz)
        d2S_dC3dB1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,1,0] + S[1,1,0,1,1,0]) / (dezz *dexy)

        d2S_dC1dC2 = (S[1,1,1,1,1,1] - S[1,1,1,0,1,1] - S[1,1,1,1,0,1] + S[1,1,1,0,0,1]) / (deyz *dexz)
        d2S_dB1dC2 = (S[1,1,1,1,1,1] - S[1,1,1,0,1,1] - S[1,1,1,1,1,0] + S[1,1,1,0,1,0]) / (deyz *dexy)

        d2S_dB1dC1 = (S[1,1,1,1,1,1] - S[1,1,1,1,0,1] - S[1,1,1,1,1,0] + S[1,1,1,1,0,0]) / (dexz *dexy)

        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess

        cos_ab = math.cos(math.pi*self.angles_guess[2]/180)
        cos_ac = math.cos(math.pi*self.angles_guess[1]/180)
        cos_bc = math.cos(math.pi*self.angles_guess[0]/180)

        ax = 1.0
        ay = 0.0
        az = 0.0
        bx = cos_ab
        by = np.sqrt(1-cos_ab**2)
        bz = 0.0
        cx = cos_ac
        cy = (cos_bc-bx*cx)/by
        cz = np.sqrt(1.0-cx**2-cy**2)
        ax = ax*a
        bx = bx*b
        by = by*b
        cx = cx*c
        cy = cy*c
        cz = cz*c
        v = self.volume_guess

        exx_n = ax/Ax0-1
        eyy_n = by/By0-1
        ezz_n = cz/Cz0-1
        exy_n = (Ax0*bx-ax*Bx0)/(Ax0*By0)
        eyz_n = (By0*cy-by*Cy0)/(By0*Cz0)
        exz_n = (Ax0*(By0*cx-bx*Cy0)-ax*(By0*Cx0-Bx0*Cy0))/(Ax0*By0*Cz0)

        dfda1 = dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dA1dC3+(exy_n-exy0)*d2F_dA1dB1+(exz_n-exz0)*d2F_dA1dC1+(eyz_n-eyz0)*d2F_dA1dC2
        dfdb2 = dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dB2dC3+(exy_n-exy0)*d2F_dB2dB1+(exz_n-exz0)*d2F_dB2dC1+(eyz_n-eyz0)*d2F_dB2dC2
        dfdc3 = dF_dC3 + (ezz_n-ezz0)*d2F_dC32+(exx_n-exx0)*d2F_dA1dC3+(eyy_n-eyy0)*d2F_dB2dC3+(exy_n-exy0)*d2F_dC3dB1+(exz_n-exz0)*d2F_dC3dC1+(eyz_n-eyz0)*d2F_dC3dC2
        dfdb1 = dF_dB1 + (exy_n-exy0)*d2F_dB12+(exx_n-exx0)*d2F_dA1dB1+(eyy_n-eyy0)*d2F_dB2dB1+(ezz_n-ezz0)*d2F_dC3dB1+(exz_n-exz0)*d2F_dB1dC1+(eyz_n-eyz0)*d2F_dB1dC2
        dfdc1 = dF_dC1 + (exz_n-exz0)*d2F_dC12+(exx_n-exx0)*d2F_dA1dC1+(eyy_n-eyy0)*d2F_dB2dC1+(ezz_n-ezz0)*d2F_dC3dC1+(exy_n-exy0)*d2F_dB1dC1+(eyz_n-eyz0)*d2F_dC1dC2
        dfdc2 = dF_dC2 + (eyz_n-eyz0)*d2F_dC22+(exx_n-exx0)*d2F_dA1dC2+(eyy_n-eyy0)*d2F_dB2dC2+(ezz_n-ezz0)*d2F_dC3dC2+(exy_n-exy0)*d2F_dB1dC2+(exz_n-exz0)*d2F_dC1dC2

        dsda1 = dS_dA1 + (exx_n-exx0)*d2S_dA12+(eyy_n-eyy0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dA1dC3+(exy_n-exy0)*d2S_dA1dB1+(exz_n-exz0)*d2S_dA1dC1+(eyz_n-eyz0)*d2S_dA1dC2
        dsdb2 = dS_dB2 + (eyy_n-eyy0)*d2S_dB22+(exx_n-exx0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dB2dC3+(exy_n-exy0)*d2S_dB2dB1+(exz_n-exz0)*d2S_dB2dC1+(eyz_n-eyz0)*d2S_dB2dC2
        dsdc3 = dS_dC3 + (ezz_n-ezz0)*d2S_dC32+(exx_n-exx0)*d2S_dA1dC3+(eyy_n-eyy0)*d2S_dB2dC3+(exy_n-exy0)*d2S_dC3dB1+(exz_n-exz0)*d2S_dC3dC1+(eyz_n-eyz0)*d2S_dC3dC2
        dsdb1 = dS_dB1 + (exy_n-exy0)*d2S_dB12+(exx_n-exx0)*d2S_dA1dB1+(eyy_n-eyy0)*d2S_dB2dB1+(ezz_n-ezz0)*d2S_dC3dB1+(exz_n-exz0)*d2S_dB1dC1+(eyz_n-eyz0)*d2S_dB1dC2
        dsdc1 = dS_dC1 + (exz_n-exz0)*d2S_dC12+(exx_n-exx0)*d2S_dA1dC1+(eyy_n-eyy0)*d2S_dB2dC1+(ezz_n-ezz0)*d2S_dC3dC1+(exy_n-exy0)*d2S_dB1dC1+(eyz_n-eyz0)*d2S_dC1dC2
        dsdc2 = dS_dC2 + (eyz_n-eyz0)*d2S_dC22+(exx_n-exx0)*d2S_dA1dC2+(eyy_n-eyy0)*d2S_dB2dC2+(ezz_n-ezz0)*d2S_dC3dC2+(exy_n-exy0)*d2S_dB1dC2+(exz_n-exz0)*d2S_dC1dC2

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_a1 = -dfda1/V*(exx_n+1) * self.eVA3_HaBohr3
        stress_b2 = -dfdb2/V*(eyy_n+1) * self.eVA3_HaBohr3
        stress_c3 = -dfdc3/V*(ezz_n+1) * self.eVA3_HaBohr3
        stress_b1 = -1.0/V*(dfdb1*(eyy_n+1)+dfda1*exy_n) * self.eVA3_HaBohr3
        stress_c2 = -1.0/V*(dfdc2*(ezz_n+1)+dfdb2*eyz_n) * self.eVA3_HaBohr3
        stress_c1 = -1.0/V*(dfdc1*(ezz_n+1)+dfdb1*eyz_n+dfda1*exz_n) * self.eVA3_HaBohr3

        if os.path.exists("elastic_constant.txt"):
            matrix_elastic = self.elastic_constants("elastic_constant.txt")
            matrix_elastic = np.array(matrix_elastic)
            matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
            M = np.array([[ d2F_dA12, d2F_dA1dB2, d2F_dA1dC3 , d2F_dA1dC2, d2F_dA1dC1, d2F_dA1dB1],
                          [d2F_dA1dB2, d2F_dB22, d2F_dB2dC3 , d2F_dB2dC2, d2F_dB2dC1, d2F_dB2dB1],
                          [d2F_dA1dC3, d2F_dB2dC3, d2F_dC32  , d2F_dC3dC2, d2F_dC3dC1, d2F_dC3dB1],
                          [d2F_dA1dC2, d2F_dB2dC2, d2F_dC3dC2, d2F_dC22 , d2F_dC1dC2, d2F_dB1dC2],
                          [d2F_dA1dC1, d2F_dB2dC1, d2F_dC3dC1, d2F_dC1dC2, d2F_dC12  , d2F_dB1dC1],
                          [d2F_dA1dB1, d2F_dB2dB1, d2F_dC3dB1, d2F_dB1dC2, d2F_dB1dC1, d2F_dB12  ]])
            P = np.array([[0.0,pressure*v,pressure*v,0,0,0],
                          [pressure*v,0.0,pressure*v,0,0,0],
                          [pressure*v,pressure*v,0.0,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0],
                          [0,0,0,0,0,0]])
            M = M + matrix_elastic+P/ self.eVA3_HaBohr3
            S1 = dsda1*(exx_n+1)
            S2 = dsdb2*(eyy_n+1)
            S3 = dsdc3*(ezz_n+1)
            S4 = (dsdc2*(ezz_n+1)+dsdb2*eyz_n)
            S5 = (dsdc1*(ezz_n+1)+dsdb1*eyz_n+dsda1*exz_n)
            S6 = (dsdb1*(eyy_n+1)+dsda1*exy_n)

            S = np.array([S1,S2,S3,S4,S5,S6])
            dstrain_dt = np.linalg.inv(M) @ S
            M=M/v*abu.eVA3_GPa
            with open("elastic.txt", "w") as f:
                f.write("Elastic [GPa] \n")
                f.write(f" \t   xx\t\tyy\t\tzz\t\tyz\t\txz\t\txy\n")
                f.write(f" xx {M[0,0]:14.8f}  {M[0,1]:14.8f}  {M[0,2]:14.8f}  {M[0,3]:14.8f}  {M[0,4]:14.8f}  {M[0,5]:14.8f}\n")
                f.write(f" yy {M[1,0]:14.8f}  {M[1,1]:14.8f}  {M[1,2]:14.8f}  {M[1,3]:14.8f}  {M[1,4]:14.8f}  {M[1,5]:14.8f}\n")
                f.write(f" zz {M[2,0]:14.8f}  {M[2,1]:14.8f}  {M[2,2]:14.8f}  {M[2,3]:14.8f}  {M[2,4]:14.8f}  {M[2,5]:14.8f}\n")
                f.write(f" yz {M[3,0]:14.8f}  {M[3,1]:14.8f}  {M[3,2]:14.8f}  {M[3,3]:14.8f}  {M[3,4]:14.8f}  {M[3,5]:14.8f}\n")
                f.write(f" xz {M[4,0]:14.8f}  {M[4,1]:14.8f}  {M[4,2]:14.8f}  {M[4,3]:14.8f}  {M[4,4]:14.8f}  {M[4,5]:14.8f}\n")
                f.write(f" xy {M[5,0]:14.8f}  {M[5,1]:14.8f}  {M[5,2]:14.8f}  {M[5,3]:14.8f}  {M[5,4]:14.8f}  {M[5,5]:14.8f}\n")
            therm = [dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),dstrain_dt[3]*(eyz_n),dstrain_dt[4]*(exz_n),dstrain_dt[5]*(exy_n)]
        else:
            therm = [0,0,0,0,0,0]


        stress[0] = stress_a1 -pressure
        stress[1] = stress_b2 -pressure
        stress[2] = stress_c3 -pressure
        stress[3] = stress_c2
        stress[4] = stress_c1
        stress[5] = stress_b1

        print ("stress")
        print (ax/Ax0, by/By0, cz/Cz0)
        print (stress)

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])
        dtol[3] = abs(stress[3]-self.stress_guess[2,1])
        dtol[4] = abs(stress[4]-self.stress_guess[2,0])
        dtol[5] = abs(stress[5]-self.stress_guess[1,0])

        return dtol, stress, therm
    # **************************************************************************************
    def stress_ZSISA_slab_1DOF(self, temp,pressure ):
        e,S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,0,0,0,0,0]
        X1 = self.ave_x[1,0,0,0,0,0]

        V = self.volume_guess

        dexx = (X0-X1)/X0
        exx0 = X1/X0-1

        dF_dX = (e[0,0,0,0,0,0]-e[2,0,0,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+e[2,0,0,0,0,0])/(dexx)**2

        x = self.ave_x_guess
        exx_n = x/X0-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_xx = -dfdx/V*(exx_n+1)*0.5 * self.eVA3_HaBohr3
        print (x/X0, x/X0, x/X0)
        print (stress_xx)

        stress[0] = stress_xx -pressure
        stress[1] = stress_xx -pressure

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])

        return dtol, stress

    def stress_ZSISA_slab_2DOF(self, temp, pressure) -> tuple:
        e,S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,1,0,0,0,0]
        Y0 = self.ave_y[1,0,0,0,0,0]

        X1 = self.ave_x[1,1,0,0,0,0]
        Y1 = self.ave_y[1,1,0,0,0,0]

        V = self.volume_guess

        dexx = (X0-X1)/X0
        deyy = (Y0-Y1)/Y0

        exx0 = X1/X0-1
        eyy0 = Y1/Y0-1

        dF_dX = (e[0,1,0,0,0,0]-e[2,1,0,0,0,0])/(2*dexx)
        dF_dY = (e[1,0,0,0,0,0]-e[1,2,0,0,0,0])/(2*deyy)

        d2F_dX2 = (e[0,1,0,0,0,0]-2*e[1,1,0,0,0,0]+e[2,1,0,0,0,0])/(dexx)**2
        d2F_dY2 = (e[1,0,0,0,0,0]-2*e[1,1,0,0,0,0]+e[1,2,0,0,0,0])/(deyy)**2
        d2F_dXdY = (e[1,1,0,0,0,0] - e[0,1,0,0,0,0] - e[1,0,0,0,0,0] + e[0,0,0,0,0,0]) / (dexx *deyy)

        x = self.ave_x_guess
        y = self.ave_y_guess

        exx_n = x/X0-1
        eyy_n = y/Y0-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(eyy_n-eyy0)*d2F_dXdY
        dfdy = dF_dY + (eyy_n-eyy0)*d2F_dY2+(exx_n-exx0)*d2F_dXdY

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_xx = -dfdx/V*(exx_n+1)* self.eVA3_HaBohr3
        stress_yy = -dfdy/V*(eyy_n+1)* self.eVA3_HaBohr3
        print (x/X0, x/X0, y/Y0)
        print (stress_xx,stress_yy)

        stress[0] = stress_xx -pressure
        stress[1] = stress_yy -pressure

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])

        return dtol, stress

    def stress_ZSISA_slab_3DOF(self, temp, pressure) -> tuple:

        e,S = self.get_vib_free_energies(temp)


        Ax0 = self.ax[0,1,1,0,0,0]
        Bx0 = self.bx[0,1,0,0,0,0]
        By0 = self.by[1,0,1,0,0,0]

        Ax1 = self.ax[1,1,1,0,0,0]
        Bx1 = self.bx[1,1,1,0,0,0]
        By1 = self.by[1,1,1,0,0,0]

        V = self.volume_guess

        dexx = (Ax0-Ax1)/Ax0
        deyy = (By0-By1)/By0
        dexy = (Ax0*(Bx0-Bx1)-(Ax0-Ax1)*Bx0)/(Ax0*By0)

        exx0 = Ax1/Ax0-1
        eyy0 = By1/By0-1
        exy0 = (Ax0*Bx1-Ax1*Bx0)/(Ax0*By0)

        dF_dA1 = (e[0,1,1,0,0,0]-e[2,1,1,0,0,0])/(2*dexx)
        dF_dB2 = (e[1,0,1,0,0,0]-e[1,2,1,0,0,0])/(2*deyy)
        dF_dB1 = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dexy)

        d2F_dA12 = (e[0,1,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,1,1,0,0,0])/(dexx)**2
        d2F_dB22 = (e[1,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[1,2,1,0,0,0])/(deyy)**2
        d2F_dB12 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dexy)**2

        d2F_dA1dB2 = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,0,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
        d2F_dA1dB1 = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,1,0,0,0,0] + e[0,1,0,0,0,0]) / (dexx *dexy)

        d2F_dB2dB1 = (e[1,1,1,0,0,0] - e[1,0,1,0,0,0] - e[1,1,0,0,0,0] + e[1,0,0,0,0,0]) / (deyy *dexy)

        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess

        cos_ab = math.cos(math.pi*self.angles_guess[2]/180)

        ax = a
        ay = 0.0
        az = 0.0
        bx = b*cos_ab
        by = b*np.sqrt(1-cos_ab**2)
        bz = 0.0
        cx = 0.0
        cy = 0.0
        cz = c

        exx_n = ax/Ax0-1
        eyy_n = by/By0-1
        exy_n = (Ax0*bx-ax*Bx0)/(Ax0*By0)

        dfda1 = dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(exy_n-exy0)*d2F_dA1dB1
        dfdb2 = dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(exy_n-exy0)*d2F_dB2dB1
        dfdb1 = dF_dB1 + (exy_n-exy0)*d2F_dB12+(exx_n-exx0)*d2F_dA1dB1+(eyy_n-eyy0)*d2F_dB2dB1

        dtol = np.zeros(6)
        stress = np.zeros(6)

        stress_a1 = -dfda1/V*(exx_n+1) * self.eVA3_HaBohr3
        stress_b2 = -dfdb2/V*(eyy_n+1) * self.eVA3_HaBohr3
        stress_b1 = -1.0/V*(dfdb1*(eyy_n+1)+dfda1*exy_n) * self.eVA3_HaBohr3
        print (ax/Ax0, by/By0)
        print (stress_a1,stress_b2,stress_b1)
        stress[0] = stress_a1 -pressure
        stress[1] = stress_b2 -pressure
        stress[5] = stress_b1

        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[5] = abs(stress[5]-self.stress_guess[1,0])

        return dtol, stress

    def cal_stress(self, temp, pressure = 0, mode = "TEC"):
        #Bohr2GPa = 29421.033
        pressure_gpa = pressure
        #pressure = pressure/Bohr2GPa
        pressure = pressure/abu.HaBohr3_GPa
        print ("Pressure = ", pressure_gpa,"GPa")
        print ("Temperature = ", temp,"K")

        if (self.sym == "v_ZSISA"):
            dtol,stress = self.stress_v_ZSISA(temp,pressure)
        elif (self.sym == "cubic" and mode == "TEC"):
            dtol,stress,therm = self.stress_ZSISA_1DOF(temp,pressure)

        elif (self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal") and mode == "TEC":
            dtol,stress,therm = self.stress_ZSISA_2DOF(temp,pressure)

        elif (self.sym == "cubic" or self.sym == "trigonal" or self.sym == "hexagonal" or self.sym == "tetragonal") and mode == "ECs":
            dtol,stress,therm = self.stress_ZSISA_3DOF(temp,pressure,mode)

        elif (self.sym == "orthorhombic"):
            dtol,stress,therm = self.stress_ZSISA_3DOF(temp,pressure,mode)

        elif (self.sym == "monoclinic"):
            dtol,stress,therm = self.stress_ZSISA_monoclinic(temp,pressure,mode)

        elif (self.sym == "triclinic"):
            dtol,stress,therm = self.stress_ZSISA_triclinic(temp,pressure,mode)

        elif (self.sym == "ZSISA_slab_1DOF"):
            dtol,stress = self.stress_ZSISA_slab_1DOF(temp,pressure)

        elif (self.sym == "ZSISA_slab_2DOF"):
            dtol,stress = self.stress_ZSISA_slab_2DOF(temp,pressure)

        elif (self.sym == "ZSISA_slab_3DOF"):
            dtol,stress = self.stress_ZSISA_slab_3DOF(temp,pressure)

        else:
            raise ValueError(f"Unknown {sym=}")

        if all(dtol[i] < 1e-8 for i in range(6)):
            with open("cell.txt", "a") as f:
                f.write(f"{temp} {pressure_gpa:.2f} {self.lattice_a_guess:.12f} {self.lattice_b_guess:.12f} {self.lattice_c_guess:.12f} {self.angles_guess[0]:.5f} {self.angles_guess[1]:.5f} {self.angles_guess[2]:.5f} {self.volume_guess:.12f} {self.ave_x_guess:.12f} {self.ave_y_guess:.12f} {self.ave_z_guess:.12f} \n")
            print("Converged !!!")
            with open("thermal.txt", "a") as f:
                f.write(f"{temp} {pressure_gpa:.2f} {therm[0]:.12e} {therm[1]:.12e} {therm[2]:.12e}  {therm[3]:.12e}  {therm[4]:.12e}  {therm[5]:.12e} \n")
            condition = True
        else :
            condition = False

        ang2bohr = 0.529177249
        with open("lat_info.txt", "w") as f:
            f.write("xred\n")
            for i in range(len(self.frac_coords_guess)):
                f.write(f"  {self.frac_coords_guess[i,0]:.12e}  {self.frac_coords_guess[i,1]:.12e}  {self.frac_coords_guess[i,2]:.12e}\n")
            f.write(f"#angdeg\n")
            f.write(f"#  {self.angles_guess[0]:.8e}  {self.angles_guess[1]:.8e}  {self.angles_guess[2]:.8e}\n")
            f.write(f"#acell \n")
            f.write(f"#  {self.lattice_a_guess/ang2bohr:.12e}  {self.lattice_b_guess/ang2bohr:.12e}  {self.lattice_c_guess/ang2bohr:.12e}\n")
            f.write(f"acell \n")
            f.write(f"  1.00 1.00 1.00\n")
            f.write(f"rprim \n")
            f.write(f"  {self.matrix_guess[0,0]/ang2bohr:.12e}  {self.matrix_guess[0,1]/ang2bohr:.12e}  {self.matrix_guess[0,2]/ang2bohr:.12e} \n")
            f.write(f"  {self.matrix_guess[1,0]/ang2bohr:.12e}  {self.matrix_guess[1,1]/ang2bohr:.12e}  {self.matrix_guess[1,2]/ang2bohr:.12e} \n")
            f.write(f"  {self.matrix_guess[2,0]/ang2bohr:.12e}  {self.matrix_guess[2,1]/ang2bohr:.12e}  {self.matrix_guess[2,2]/ang2bohr:.12e} \n")
            f.write(f"strtarget \n")
            f.write(f"  {stress[0]:.12e}  {stress[1]:.12e}  {stress[2]:.12e} \n")
            f.write(f"  {stress[3]:.12e}  {stress[4]:.12e}  {stress[5]:.12e} \n")
        return condition

    def get_vib_free_energies(self, temp) -> tuple:

        f = np.zeros((self.dim[0],self.dim[1],self.dim[2],self.dim[3],self.dim[4],self.dim[5]))
        entropy = np.zeros((self.dim[0],self.dim[1],self.dim[2],self.dim[3],self.dim[4],self.dim[5]))

        for i, dim0_list in enumerate(self.phdoses):
            for j, dim1_list in enumerate(dim0_list):
                for k, dim2_list in enumerate(dim1_list):
                    for l, dim3_list in enumerate(dim2_list):
                        for m, dim4_list in enumerate(dim3_list):
                            for n, dos in enumerate(dim4_list):
                               if dos is not None:
                                   f[i,j,k,l,m,n] = dos.get_free_energy(temp, temp, 1).values
                                   entropy[i,j,k,l,m,n] = dos.get_entropy(temp, temp, 1).values
                               else:
                                   f[i,j,k,l,m,n] = None
                                   entropy[i,j,k,l,m,n] = None
        return f, entropy

    def elastic_constants(self, file_name):
        with open(file_name, "rt") as file:
            lines = file.readlines()

        in_relaxed_section = False
        elastic_data = []

        for line in lines:
            if "[ELASTIC_RELAXED]" in line:
                in_relaxed_section = True
                continue
            if "[" in line and in_relaxed_section:  # Stop if another section starts
                break
            if in_relaxed_section:
                try:
                    # Extract numbers from each line
                    numbers = list(map(float, line.split()[1:7]))
                    if numbers:  # Add non-empty lines
                        elastic_data.append(numbers)
                except ValueError:
                    continue

        # Convert the list into a numpy array for easier indexing
        matrix_elastic = np.array(elastic_data)

        # Extract the desired values
        return matrix_elastic
