"""Database of unit cells in the ABINIT format."""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.abilab import Structure, ArrayWithUnit

__all__ = [
    "ucell_names",
    "ucell",
    "structure_from_ucell",
]


# Public API
def ucell_names():
    """List with the name of the entries."""
    return list(_UCELLS.keys())


def ucell(name):
    """Return the entry in the database with the given name."""
    return _UCELLS[name].copy()


def structure_from_ucell(name):
    """Returns a `Structure` from the name of entry in the database."""
    try:
        return Structure.from_abivars(ucell(name))
    except KeyError:
        raise KeyError("Cannot find key %s in:\n %s" % (name, list(_UCELLS.keys())))


_UCELLS = {
    "Si": dict(ntypat=1,
               natom=2,
               typat=[1, 1],
               znucl=14,
               acell=3*[10.217],
               rprim=[[0.0,  0.5,  0.5],
                      [0.5,  0.0,  0.5],
                      [0.5,  0.5,  0.0]],
               xred=[ [0.0 , 0.0 , 0.0],
                      [0.25, 0.25, 0.25]],
                    ),
    "Si-shifted": dict(ntypat=1,
               natom=2,
               typat=[1, 1],
               znucl=14,
               acell=3*[10.217],
               rprim=[[0.0,  0.5,  0.5],
                      [0.5,  0.0,  0.5],
                      [0.5,  0.5,  0.0]],
               xred=[ [1.0 , 1.0 , 1.0],
                      [0.25, 0.25, 0.25]],
                    ),
    "Al": dict(ntypat=1,
               natom=1,
               typat=[1],
               znucl=13,
               acell=3*[7.60],
               rprim=[[0.0,  0.5,  0.5],
                      [0.5,  0.0,  0.5],
                      [0.5,  0.5,  0.0]],
               xred=[ [0.0 , 0.0 , 0.0]],
                    ),

    "Al-negative-volume": dict(ntypat=1,
               natom=1,
               typat=[1],
               znucl=13,
               acell=3*[7.60],
               rprim=[[0.0,  0.5,  0.5],
                      [-0.5,  -0.0,  -0.5],
                      [0.5,  0.5,  0.0]],
               xred=[ [0.0 , 0.0 , 0.0]],
                    ),

    "ZnO": dict(ntypat=2,
                natom=2,
                typat=[1, 2],
                acell= 3*[8.6277],
                rprim= [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]],
                znucl=[30, 8],
                xred=[[.0, .0, .0], [.25,.25,.25]],
    ),

    "SiC": dict(ntypat=2,
                natom=2,
                typat=[1, 2],
                acell=3*[8.19],
                rprim=[[.0, .5, .5],
                       [.5, .0, .5],
                       [.5, .5, .0]],
                znucl=[6, 14],
                xred=[ [.0, .0, .0],
                       [.25,.25,.25] ]
                ),

    "AlAs": dict(natom=2,
                 typat=[1, 2],
                 acell=3*[10.61],
                 rprim=[[0.0,  0.5,  0.5],
                        [0.5,  0.0,  0.5],
                        [0.5,  0.5,  0.0]],
                 znucl=[13, 33],
                 xred=[[0.0,  0.0,  0.0],
                       [0.25, 0.25, 0.25]]
                ),

    "GaAs": dict(natom=2,
                 typat=[1, 2],
                 acell=3*[10.60],
                 rprim=[[0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5],
                        [0.5, 0.5, 0.0]],
                 znucl=[31, 33],
                       xred=[3 *[0.00],
                             3 *[0.25]]
                ),

    "NiO": dict(natom=4,
                typat=[1, 1, 2, 2],
                acell=3 * [7.92],
                rprim=[0.0, 1/2, 1/2,
                       1/2, 0.0, 1/2,
                       1.0, 1.0, 0.0],
                znucl=[28, 8],
                xred=[0.0, 0.0, 0.0,
                      0.0, 0.0, 0.5,
                      0.5, 0.5, 0.25,
                      0.5, 0.5, 0.75
                      ]
                ),
    "MgB2": dict(natom=3,
                 typat=[1, 2, 2],
                 acell=ArrayWithUnit([3.086, 3.086, 3.523], "ang").to("bohr"),
                 rprim= [ 0.866025403784439, 0.5, 0.0,
                         -0.866025403784439, 0.5, 0.0,
                          0.0              , 0.0, 1.0],
                 znucl=[12, 5],
                 xred=[0.0, 0.0, 0.0, # Mg
                       1/3, 2/3, 0.5, # B
                       2/3, 1/3, 0.5] # B
                ),
    "Fe-fm":  dict(natom=1,
                   typat=1,
                   acell=3*[5.42],
                   rprim=[-0.5,  0.5,  0.5,
                           0.5, -0.5,  0.5,
                           0.5,  0.5, -0.5],
                   znucl=26,
                   xred=[0.0, 0.0, 0.0],
                  ),
    "SiO2-alpha": dict(natom=9,
                  typat=3*[1] + 6*[2],
                  acell=ArrayWithUnit([4.91304,  4.91304, 5.40463], "ang").to("bohr"),
                  rprim=[ 5.0000000000e-01, -8.6602540378e-01,  0.0000000000e+00,
                          5.0000000000e-01,  8.6602540378e-01,  0.0000000000e+00,
                          0.0000000000e+00,  0.0000000000e+00,  1.0000000000e+00],
                  znucl=[14, 8],
                  # Experimental parameters (Wyckoff pag 312)
                  # u(Si)= 0.465
                  # x= 0.415 ; y= 0.272 ; z= 0.120
                  xred=[ 0.465,   0.000,   0.000             , #Si
                         0.000,   0.465,   2/3               , #Si
                        -0.465,  -0.465,   1/3               , #Si
                         0.415,   0.272,   0.120             , #O
                        -0.143,  -0.415,   0.4533333333333333, #O
                        -0.272,   0.143,   0.7866666666666666, #O
                         0.143,  -0.272,  -0.120             , #O
                         0.272,   0.415,   0.5466666666666666, #O
                        -0.415,  -0.143,   0.2133333333333333] #O
             ),
}

#dict(natom=,
#     typat=,
#     acell=,
#     rprim=,
#     znucl=,
#     xred=,
#     ),
