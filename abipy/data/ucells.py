"""Database of unit cells in the ABINIT format."""
from __future__ import division, print_function

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
    """Returnn the entry in the database with the given name."""
    return _UCELLS[name].copy()


def structure_from_ucell(name):
    """Returns a `Structure` from the name of entry in the database."""
    return Structure.from_abivars(ucell(name))


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
                 xred=[0.0, 0.0, 0.0, #Mg
                       1/3, 2/3, 0.5, #B
                       2/3, 1/3, 0.5] #B
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
}

#dict(natom=,                       
#     typat=,                       
#     acell=,
#     rprim=,
#     znucl=,
#     xred=,
#     ),
