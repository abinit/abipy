"""
Tools for writing cube files.
See http://paulbourke.net/dataformats/cube/ and http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
"""

from __future__ import annotations

import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.units import bohr_to_angstrom

__all__ = [
    "cube_write_data",
    "cube_write_structure_mesh",
]


def cube_write_structure_mesh(file, structure, mesh) -> None:
    """
    Write structure and mesh information to a cube file.

    Args:
        file: File-like object.
        structure: Structure object.
        mesh: Mesh3D object.
    """
    fwrite = file.write
    fwrite("Density generated from abipy\n")
    fwrite("in the cube file format\n")
    dvx = mesh.dvx / bohr_to_angstrom
    dvy = mesh.dvy / bohr_to_angstrom
    dvz = mesh.dvz / bohr_to_angstrom
    fwrite(f"{len(structure):4d} {0.0:.6f} {0.0:.6f} {0.0:.6f}\n")
    fwrite(f"{mesh.nx:4d} {dvx[0]:.6f} {dvx[1]:.6f} {dvx[2]:.6f}\n")
    fwrite(f"{mesh.ny:4d} {dvy[0]:.6f} {dvy[1]:.6f} {dvy[2]:.6f}\n")
    fwrite(f"{mesh.nz:4d} {dvz[0]:.6f} {dvz[1]:.6f} {dvz[2]:.6f}\n")
    for site in structure:
        cc = site.coords / bohr_to_angstrom
        fwrite(f"{site.specie.Z:d} {0.0:.10f} {cc[0]:.10f} {cc[1]:.10f} {cc[2]:.10f}\n")


def cube_write_data(file, data, mesh) -> None:
    """
    Write data to a cube file.

    Args:
        file: File-like object.
        data: Numpy array with the data.
        mesh: Mesh3D object.
    """
    fwrite = file.write
    data_bohrs = data * (bohr_to_angstrom**3)
    for ix in range(mesh.nx):
        for iy in range(mesh.ny):
            for iz in range(mesh.nz):
                fwrite(f"{data_bohrs[ix, iy, iz]:.5e}\n")


def cube_read_structure_mesh_data(filepath: str) -> tuple:
    with open(filepath) as fh:
        # The two first lines are comments
        for ii in range(2):
            fh.readline()
        # Number of atoms
        natoms = int(fh.readline().split()[0])
        # The next three lines give the mesh and the vectors
        sp = fh.readline().split()
        nx = int(sp[0])
        dvx = np.array([float(sp[ii]) for ii in range(1, 4)]) * bohr_to_angstrom
        sp = fh.readline().split()
        ny = int(sp[0])
        dvy = np.array([float(sp[ii]) for ii in range(1, 4)]) * bohr_to_angstrom
        sp = fh.readline().split()
        nz = int(sp[0])
        dvz = np.array([float(sp[ii]) for ii in range(1, 4)]) * bohr_to_angstrom
        uc_matrix = np.array([nx * dvx, ny * dvy, nz * dvz])
        sites = []
        lattice = Lattice(uc_matrix)
        for ii in range(natoms):
            sp = fh.readline().split()
            cc = np.array([float(sp[ii]) for ii in range(2, 5)]) * bohr_to_angstrom
            sites.append(
                PeriodicSite(int(sp[0]), coords=cc, lattice=lattice, to_unit_cell=False, coords_are_cartesian=True)
            )
        data = np.zeros((nx, ny, nz))
        ii = 0
        for line in fh:
            for val in line.split():
                data[ii // (ny * nz), (ii // nz) % ny, ii % nz] = float(val)
                ii += 1
        data = data / (bohr_to_angstrom**3)
        if ii != nx * ny * nz:
            raise ValueError("Wrong number of data points ...")
        from abipy.core.structure import Structure

        structure = Structure.from_sites(sites=sites)
        from abipy.core.mesh3d import Mesh3D

        mesh = Mesh3D(shape=[nx, ny, nz], vectors=uc_matrix)

        return structure, mesh, data
