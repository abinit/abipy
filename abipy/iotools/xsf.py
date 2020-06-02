# coding: utf-8
"""Tools for writing Xcrysden files."""

import numpy as np

from pymatgen.core.units import Energy, EnergyArray #, ArrayWithUnit
from abipy.tools.numtools import transpose_last3dims, add_periodic_replicas


__all__ = [
    "xsf_write_structure_and_data_to_path",
    "xsf_write_structure",
    "xsf_write_data",
    "bxsf_write",
]


def xsf_write_structure(file, structures):
    """
    Write the crystalline structure in the Xcrysden format (XSF)

    Args:
        file: file-like object.
        structures: :class:`Structure` or list of :class:`Structure` objects.
    """
    animation = True
    if not isinstance(structures, (list, tuple)):
        structures = [structures]
        animation = False

    fwrite = file.write

    if animation:
        # axsf file.
        fwrite('ANIMSTEPS %s\n' % len(structures))

    fwrite('CRYSTAL\n')

    for n, struct in enumerate(structures):
        cell = struct.lattice.matrix

        fwrite('# Primitive lattice vectors in Angstrom\n')
        fwrite('PRIMVEC %d\n' % (n + 1))
        for i in range(3):
            fwrite(' %.14f %.14f %.14f\n' % tuple(cell[i]))

        cart_coords = struct.cart_coords
        atomic_numbers = struct.atomic_numbers

        # TODO
        cart_forces = None
        #if "cartesian_forces" in structure.site_properties:
        #cart_forces = ArrayWithUnit().to("Ha ang^-1")

        fwrite("# Cartesian coordinates in Angstrom.\n")
        fwrite('PRIMCOORD %d\n' % (n + 1))
        fwrite(' %d 1\n' % len(cart_coords))

        for a in range(len(cart_coords)):
            fwrite(' %2d' % atomic_numbers[a])
            fwrite(' %20.14f %20.14f %20.14f' % tuple(cart_coords[a]))
            if cart_forces is None:
                fwrite('\n')
            else:
                fwrite(' %20.14f %20.14f %20.14f\n' % tuple(cart_forces[a]))


def xsf_write_structure_and_data_to_path(filepath, structure, datar, **kwargs):
    """Simplified interface to xsf routines to write structure and data to filepath."""
    with open(filepath, mode="wt") as fh:
        xsf_write_structure(fh, structure)
        xsf_write_data(fh, structure, datar, **kwargs)


def xsf_write_data(file, structure, data, add_replicas=True, cplx_mode=None):
    """
    Write data in the Xcrysden format (XSF)

    Args:
        file: file-like object.
        structure: :class:`Structure` object.
        data: array-like object in C-order, i.e data[nx, ny, nz]
        add_replicas: If True, data is padded with redundant data points.
            in order to have a periodic 3D array of shape: (nx+1, ny+1, nz+1).
        cplx_mode: string defining the data to print when data is a complex array.
            Possible choices are (case-insensitive):

                - "re"  for real part.
                - "im" for imaginary part.
                - "abs" for the absolute value
    """
    fwrite = file.write

    # Check this one
    if add_replicas:
        data = add_periodic_replicas(data)

    if np.iscomplexobj(data):
        if cplx_mode is None:
            raise TypeError("cplx_mode must be specified when data is a complex array.")
        cplx_mode = cplx_mode.lower()
        if cplx_mode == "re":
            data = data.real
        elif cplx_mode == "im":
            data = data.imag
        elif cplx_mode == "abs":
            data = np.abs(data)
        else:
            raise ValueError("Wrong value for cplx_mode: %s" % cplx_mode)

    shape, ndim = data.shape, data.ndim

    if ndim == 3:
        ngrids = 1
        data = np.asarray([data])
    elif ndim == 4:
        ngrids = shape[0]
    else:
        raise ValueError("ndim %d is not supported" % ndim)

    # Xcrysden uses Fortran-order.
    # Transpose (...,x,y,z) --> (...,z,y,x) to speed up the write below.
    fdata = transpose_last3dims(data)
    fgrid = fdata.shape[-3:]

    cell = structure.lattice_vectors(space="r")
    origin = np.zeros(3)

    fwrite('BEGIN_BLOCK_DATAGRID_3D\n')
    fwrite(' data\n')

    for dg in range(ngrids):
        fwrite(" BEGIN_DATAGRID_3Dgrid#" + str(dg+1) + "\n")
        fwrite('%d %d %d\n' % shape[-3:])

        fwrite('%f %f %f\n' % tuple(origin))
        for i in range(3):
            fwrite('%f %f %f\n' % tuple(cell[i]))

        for z in range(fgrid[0]):
            for y in range(fgrid[1]):
                slice_x = fdata[dg,z,y]
                fwrite(' '.join(['%f' % d for d in slice_x]))
                fwrite('\n')
            fwrite('\n')

        fwrite(' END_DATAGRID_3D\n')
    fwrite('END_BLOCK_DATAGRID_3D\n')


def bxsf_write(file, structure, nsppol, nband, ndivs, ucdata_sbk, fermie, unit="eV"):
    """
    Write band structure data in the Xcrysden format (XSF)

    Args:
        file: file-like object.
        structure: :class:`Structure` object.
        nsppol: Number of spins.
        nband: Number of bands.
        ndivs: Number of divisions of the full k-mesh.
        ucdata_sbk: Array [nsppol, nband, ndivs[0], ndivs[1], mpdvis[2]] with energies
            in the unic cell mesh in unit `unit`.
        fermie: Fermi energy.
        unit=Unit of input `ucdata_sbk` and `fermie`. Energies will be converted to Hartree before writing.

    .. note::

        #. The k-points must span the reciprocal unit cell, not the Brillouin zone.

        #. The mesh must be closed and centered on Gamma.

        #. Energies are written in row-major (i.e. C) order.

        #  Energies are in Hartree.

    See also http://www.xcrysden.org/doc/XSF.html
    """
    # Xscryden uses Ha for energies.
    ucdata_sbk = EnergyArray(ucdata_sbk, unit).to("Ha")
    fermie = Energy(fermie, unit).to("Ha")

    ucdata_sbk = np.reshape(ucdata_sbk, (nsppol, nband, np.product(ndivs)))

    close_it = False
    if not hasattr(file, "write"):
        file = open(file, mode="w")
        close_it = True

    fw = file.write

    # Write the header.
    fw('BEGIN_INFO\n')
    fw('# Band-XCRYSDEN-Structure-File for Visualization of Fermi Surface generated by the ABINIT package\n')
    fw('# NOTE: the first band is relative to spin-up electrons,\n')
    fw('#       the second band to spin-down electrons (if any) and so on ...\n#\n')
    fw('# Launch as: xcrysden --bxsf\n#\n')
    fw(' Fermi Energy: %f\n' % fermie)
    fw('END_INFO\n\n')
    fw('BEGIN_BLOCK_BANDGRID_3D\n')
    fw(' band_energies\n')
    fw(' BEGIN_BANDGRID_3D\n')

    fw(str(nsppol * nband) + "\n")  # Number of bands written.
    fw("%d %d %d\n" % tuple(ndivs)) # Number of division in the full BZ mesh.
    fw("0 0 0\n")                   # Unshifted meshes are not supported.

    # Reciprocal lattice vectors in Ang^{-1}
    gcell = structure.lattice_vectors("g")
    for i in range(3):
        fw('%f %f %f\n' % tuple(gcell[i]))

    # Write energies on the full mesh for all spins and bands.
    idx = 0
    for band in range(nband):
        for spin in range(nsppol):
            idx += 1
            enes = ucdata_sbk[spin, band, :]
            fw(" BAND: %d\n" % idx)
            fw("\n".join("%.18e" % v for v in enes))
            fw("\n")

    fw(' END_BANDGRID_3D\n')
    fw('END_BLOCK_BANDGRID_3D\n')
    file.flush()

    if close_it:
        file.close()
