# coding: utf-8
"""
Tools for writing cube files.
See http://paulbourke.net/dataformats/cube/ and http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
"""

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.units import bohr_to_angstrom


__all__ = [
    "cube_write_structure_mesh",
    "cube_write_data",
]


def cube_write_structure_mesh(file, structure, mesh):
    fwrite = file.write
    fwrite("Density generated from abipy\n")
    fwrite("in the cube file format\n")
    dvx = mesh.dvx / bohr_to_angstrom
    dvy = mesh.dvy / bohr_to_angstrom
    dvz = mesh.dvz / bohr_to_angstrom
    fwrite('{:4d} {:.6f} {:.6f} {:.6f}\n'.format(len(structure), 0.0, 0.0, 0.0))
    fwrite('{:4d} {:.6f} {:.6f} {:.6f}\n'.format(mesh.nx, dvx[0], dvx[1], dvx[2]))
    fwrite('{:4d} {:.6f} {:.6f} {:.6f}\n'.format(mesh.ny, dvy[0], dvy[1], dvy[2]))
    fwrite('{:4d} {:.6f} {:.6f} {:.6f}\n'.format(mesh.nz, dvz[0], dvz[1], dvz[2]))
    for site in structure:
        cc = site.coords / bohr_to_angstrom
        fwrite('{:d} {:.10f} {:.10f} {:.10f} {:.10f}\n'.format(site.specie.Z, 0.0, cc[0], cc[1], cc[2]))


def cube_write_data(file, data, mesh):
    fwrite = file.write
    data_bohrs = data * (bohr_to_angstrom ** 3)
    for ix in range(mesh.nx):
        for iy in range(mesh.ny):
            for iz in range(mesh.nz):
                fwrite('{:.5e}\n'.format(data_bohrs[ix, iy, iz]))


def cube_read_structure_mesh_data(file):
    with open(file, 'r') as fh:
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
        uc_matrix = np.array([nx*dvx, ny*dvy, nz*dvz])
        sites = []
        lattice = Lattice(uc_matrix)
        for ii in range(natoms):
            sp = fh.readline().split()
            cc = np.array([float(sp[ii]) for ii in range(2, 5)]) * bohr_to_angstrom
            sites.append(PeriodicSite(int(sp[0]), coords=cc, lattice=lattice, to_unit_cell=False,
                                      coords_are_cartesian=True))
        data = np.zeros((nx, ny, nz))
        ii = 0
        for line in fh:
            for val in line.split():
                data[ii // (ny * nz), (ii // nz) % ny, ii % nz] = float(val)
                ii += 1
        data = data / (bohr_to_angstrom ** 3)
        if ii != nx*ny*nz:
            raise ValueError('Wrong number of data points ...')
        from abipy.core.structure import Structure
        structure = Structure.from_sites(sites=sites)
        from abipy.core.mesh3d import Mesh3D
        mesh = Mesh3D(shape=[nx, ny, nz], vectors=uc_matrix)

        return structure, mesh, data
