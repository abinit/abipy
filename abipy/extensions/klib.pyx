#!python
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import print_function, division

import numpy as np
cimport numpy as np

cimport cython

cdef extern from "math.h":
    double floor(double x)

cdef inline int nint(double num):
    """Truncate to the nearest integer."""
    return int(floor(num + 0.5))


cdef inline timrotk(int time_sign, const int symrec[3][3], const double kpt[3], double krot[3]):
    """Rotate the kpoint kpt in reduced coordinates"""
    cdef int i, j

    for i in range(3):
        krot[i] = 0.0
        for j in range(3):
            krot[i] += symrec[i][j] * kpt[j]
        krot[i] *= time_sign


def map_bz2ibz(structure, 
               np.ndarray[dtype=np.double_t, ndim=2] bz, 
               np.ndarray[dtype=np.double_t, ndim=2] ibz,
               atol=1.e-8
               ):

    cdef:
        unsigned int ik_bz, ik_ibz, found, i, j
        int nsym = len(structure.spacegroup), len_bz=len(bz), len_ibz=len(ibz)
        int time_sign, isym
        double kdiff[3], kbz[3], kibz[3], krot[3]
        int int_kdiff[3] #, rot_g[3][3]
        double _atol = atol
        #np.ndarray np_rot_g = np.ones((3,3), dtype=np.int)
        #np.ndarray bz2ibz = -np.ones(len(bz), dtype=np.int)
        # FIXME
        int bz2ibz[100000000] 

        #FIXME
        int nsym_tr = len(structure.spacegroup)
        int rot_g[128][3][3]

    for isym, symmop in enumerate(structure.spacegroup):
        time_sign = symmop.time_sign
        rg = symmop.rot_g
        rot_g[isym][0][0] = time_sign * rg[0,0]
        rot_g[isym][0][1] = time_sign * rg[0,1]
        rot_g[isym][0][2] = time_sign * rg[0,2]
        rot_g[isym][1][0] = time_sign * rg[1,0]
        rot_g[isym][1][1] = time_sign * rg[1,1]
        rot_g[isym][1][2] = time_sign * rg[1,2]
        rot_g[isym][2][0] = time_sign * rg[2,0]
        rot_g[isym][2][1] = time_sign * rg[2,1]
        rot_g[isym][2][2] = time_sign * rg[2,2]

    for ik_bz in range(len_bz):
        kbz[0] = bz[ik_bz,0]
        kbz[1] = bz[ik_bz,1]
        kbz[2] = bz[ik_bz,2]

        found = 0
        for ik_ibz in range(len_ibz):
            if found != 0: break
            kibz[0] = ibz[ik_ibz,0]
            kibz[1] = ibz[ik_ibz,1]
            kibz[2] = ibz[ik_ibz,2]

            for isym in range(nsym_tr):

                #krot = np.dot(rot_g, kibz) * time_sign
                for i in range(3):
                    krot[i] = 0.0
                    for j in range(3):
                        krot[i] += rot_g[isym][i][j] * kibz[j]
                    krot[i] = krot[i] * time_sign

                kdiff[0] = krot[0] - kbz[0]
                kdiff[1] = krot[1] - kbz[1]
                kdiff[2] = krot[2] - kbz[2]

                #int_kdiff[0] = int(round(kdiff[0]))
                #int_kdiff[1] = int(round(kdiff[1]))
                #int_kdiff[2] = int(round(kdiff[2]))

                #int_kdiff[0] = int(floor(kdiff[0]+0.5))
                #int_kdiff[1] = int(floor(kdiff[1]+0.5))
                #int_kdiff[2] = int(floor(kdiff[2]+0.5))

                int_kdiff[0] = nint(kdiff[0])
                int_kdiff[1] = nint(kdiff[1])
                int_kdiff[2] = nint(kdiff[2])

                if (abs(int_kdiff[0] - kdiff[0]) < _atol and 
                    abs(int_kdiff[1] - kdiff[1]) < _atol and 
                    abs(int_kdiff[2] - kdiff[2]) < _atol):

                    bz2ibz[ik_bz] = ik_ibz
                    #ktabi[ik_bz] = isym
                    #ktabo[ik_bz] = itime
                    #umklp[ik_bz] = g0
                    found = 1
                    break

    #return bz2ibz


cdef _capi_mesh2ibz(const int len_bz,
                    const double bz_arr[][3], 
                    const int len_ibz,
                    const double ibz_arr[][3],
                    const int nsym,
                    const int timrev,
                    const int symrec[][3][3],
                    const double atol,
                    int bz2ibz[],
                    #int bz2ibz[len_bz],
                    #int bz2ibz[len_bz],
                    ):

    cdef:
        unsigned int ik_bz, ik_ibz, found, i, j, isym, itime
        int time_sign, 
        int int_kdiff[3] 
        double kdiff[3], kbz[3], kibz[3], krot[3]

    for ik_bz in range(len_bz):
        kbz[0] = bz_arr[ik_bz][0]
        kbz[1] = bz_arr[ik_bz][1]
        kbz[2] = bz_arr[ik_bz][2]

        found = 0
        for ik_ibz in range(len_ibz):
            if found != 0: break
            kibz[0] = ibz_arr[ik_ibz][0]
            kibz[1] = ibz_arr[ik_ibz][1]
            kibz[2] = ibz_arr[ik_ibz][2]

            for itime in range(timrev): 
                time_sign = +1
                if itime == 2: time_sign = -1

                for isym in range(nsym):
                    for i in range(3):
                        krot[i] = 0.0
                        for j in range(3):
                            krot[i] += symrec[isym][i][j] * kibz[j]
                        krot[i] = krot[i] * time_sign

                    kdiff[0] = krot[0] - kbz[0]
                    kdiff[1] = krot[1] - kbz[1]
                    kdiff[2] = krot[2] - kbz[2]

                    int_kdiff[0] = int(floor(kdiff[0]+0.5))
                    int_kdiff[1] = int(floor(kdiff[1]+0.5))
                    int_kdiff[2] = int(floor(kdiff[2]+0.5))

                    if (abs(int_kdiff[0] - kdiff[0]) < atol and 
                        abs(int_kdiff[1] - kdiff[1]) < atol and 
                        abs(int_kdiff[2] - kdiff[2]) < atol):

                        bz2ibz[ik_bz] = ik_ibz
                        #ktabo[ik_bz] = isym
                        #ktabi[ik_bz] = itime
                        #umklp[ik_bz,:] = g0
                        found = 1
                        break
