#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from abipy.core.kpoints import issamek

def slow_mesh2ibz(structure, bz, ibz):
    #print("bz",bz)
    #print("ibz",ibz)
    bz2ibz = -np.ones(len(bz), dtype=np.int)

    for ik_bz, kbz in enumerate(bz):
        #print("kbz",kbz)
        found = False

        for ik_ibz, kibz in enumerate(ibz):
            #print("kibz",kibz)
            if found: break
            for symmop in structure.spacegroup:
                krot = symmop.rotate_k(kibz.frac_coords)
                if issamek(krot, kbz):
                    bz2ibz[ik_bz] = ik_ibz
                    found = True
                    break

    return bz2ibz


def map_mesh2ibz(structure, mpdivs, shifts, ibz):
    """
    This function computes the mapping between the
    points in the full BZ and the points in the IBZ.

    Args:
        structure:
            pymatgen `Structure` instance.
        mpdivs
            The three MP divisions
        shifts:
            Array-like object with the MP shift.
        ibz:
            array-like object with the reduced coordinates of the
            points in the IBZ.

    Returns
        ndarray array that maps the full mesh onto ibz i.e.

            bzmap[idx_bz] = id_ibz

        where idx_bz is the index the BZ and idx_ibz the index in the IBZ.
    """
    # Build bz grid.
    mpdivs = np.asarray(mpdivs)

    shifts = np.reshape(shifts, (-1,3))
    num_shifts = len(shifts)

    from collections import namedtuple
    class KSymmetryTables(namedtuple("KSymmetryTables", "bz2ibz ktabi ktabo")):
        pass

    spg = structure.spacegroup

    timrev = 2 if spg.has_timerev else 1
    spg_fdata = spg.to_fortran_arrays()

    # Build mapping.
    ibz_fdata = ibz.to_fortran_arrays()

    tables = num_shifts * [None]

    for ish, shift in enumerate(shifts):
        kbz, count = np.empty((mpdivs.prod(), 3)), 0

        for i in range(mpdivs[0]):
            x = (i + shift[0]) / mpdivs[0]
            for j in range(mpdivs[1]):
                y = (j + shift[1]) / mpdivs[1]
                for k in range(mpdivs[2]):
                    z = (k + shift[2]) / mpdivs[2]
                    kbz[count] = (x, y, z)
                    count += 1

        bz2ibz, ktabi, ktabo, ierr = mod.map_bz2ibz(kibz=ibz_fdata.frac_coords,
                                                    kbz=np.asfortranarray(kbz.T),
                                                    timrev=timrev,
                                                    symrec=spg_fdata.symrec,
                                                    symafm=spg_fdata.symafm
                                                    )
        if ierr != 0:
            msg = "map_bz2ibz returned ierr %d, Kmesh is not symmetric" % ierr
            #raise KmeshNotSymmetricError(msg)
            raise ValueError(msg)

        #for (ik_bz, full_kpt) in enumerate(kbz):
        #    ik_ibz = bz2ibz[ik_bz]
        #    isym = ktabo[ik_bz]
        #    itime = ktabi[ik_bz]
        #    op = spg[isym + itime
        #    krot = op.rotate_k(ibz[ik_ibz].frac_coords, wrap_tows=True)
        #    assert full_kpt == krot

        tables[ish] = KSymmetryTables(bz2ibz=bz2ibz, ktabi=ktabi, ktabo=ktabo)

    return tuple(tables)


def main():
    #import m_kpoints
    #from m_kpoints import m_kpoints  as mod
    #print(m_kpoints.__doc__)
    #v1 = np.zeros(3)
    #v2 = v1 + 1.0
    #print(v1, v2)

    #ans, g0 = mod.isamek(v1, v2)
    #print(ans,g0)

    from abipy.core.kpoints import IrredZone, kmesh_from_mpdivs

    filename = "/Users/gmatteo/Coding/abipy/abipy/data/runs/data_si_ebands/outdata/si_scf_GSR.nc"
    structure = Structure.from_file(filename)

    ibz = IrredZone.from_file(filename)
    #print(ibz)

    mpdivs = [8, 8, 8]
    shifts = [0.0, 0.0, 0.0]

    #tables = map_mesh2ibz(structure, mpdivs, shifts, ibz)

    bz = kmesh_from_mpdivs(mpdivs, shifts)

    ktab = slow_map_mesh2ibz(structure, bz, ibz.frac_coords)

if __name__ == "__main__":
    main()
