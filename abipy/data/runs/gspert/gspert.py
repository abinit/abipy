#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.data as abidata  
import abipy.abilab as abilab

"""
# H2 molecule in a big box
#
# In this input file, the location of the information on this or that line
# is not important : a keyword is located by the parser, and the related
# information should follow. 
# The "#" symbol indicates the beginning of a comment : the remaining
# of the line will be skipped.

ndtset 2
ecut1 20
ecut2 40

nstep 40          # Maximal number of SCF cycles

getwfk2 -1
nnsclo2 1
nstep2 1
nline2 1
wfoptalg2 1
# nbdbuf 2
useria2 1
userib2 0
# userib = 0 : nothing else
# userib = 1 : normalization
# userib = 2 : ortho
# userib = 3 : RR
# iscf 17 # 7, default : no double counting. 17 : double counting
# nband 10


#Definition of the unit cell
acell 10 10 10    # The keyword "acell" refers to the
                  # lengths of the primitive vectors (in Bohr)
#rprim 1 0 0  0 1 0  0 0 1 # This line, defining orthogonal primitive vectors,
                           # is commented, because it is precisely the default value of rprim

#Definition of the atom types
ntypat 1          # There is only one type of atom
znucl 1           # The keyword "znucl" refers to the atomic number of the 
                  # possible type(s) of atom. The pseudopotential(s) 
                  # mentioned in the "files" file must correspond
                  # to the type(s) of atom. Here, the only type is Hydrogen.
                         

#Definition of the atoms
natom 2           # There are two atoms
typat 1 1         # They both are of type 1, that is, Hydrogen
xcart             # This keyword indicates that the location of the atoms
                  # will follow, one triplet of number for each atom
  -0.7 0.0 0.0    # Triplet giving the cartesian coordinates of atom 1, in Bohr
   0.7 0.0 0.0    # Triplet giving the cartesian coordinates of atom 2, in Bohr

#Definition of the planewave basis set

#Definition of the k-point grid
kptopt 0          # Enter the k points manually 
nkpt 1            # Only one k point is needed for isolated system,
                  # taken by default to be 0.0 0.0 0.0

#Definition of the SCF procedure
# toldfe 1.0d-14     # Will stop when, twice in a row, the difference 
tolwfr 1.0d-30     # Will stop when, twice in a row, the difference 
                  # between two consecutive evaluations of total energy 
                  # differ by less than toldfe (in Hartree) 
                  # This value is way too large for most realistic studies of materials
"""


def make_inputs(structure, pseudos, ecuts, nksmall, tsmear=None):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    assert len(ecuts) == 2
    structure = abilab.Structure.as_structure(structure)
    multi = abilab.MultiDataset(structure=structure, pseudos=pseudos, ndtset=len(ecuts))

    # Get the number of valence electrons, add 5 and use this value to compute nband
    nval = structure.num_valence_electrons(pseudos)
    nelec = nval + 5
    # TODO: Mind the smearing

    ksamp = structure.calc_ksampling(nksmall)
    multi.set_kmesh(ngkpt=ksamp.ngkpt, shiftk=ksamp.shiftk)
    #multi.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])

    # Set ecut.
    for i, ecut in enumerate(ecuts):
      multi[i].set_vars(ecut=ecut)

    # Global variables
    multi.set_vars(
      nstep=40,
      tolwfr=1e-20,
      nband=int(nval / 2),
      #iscf=17,
    )
    
    if tsmear is not None:
        multi.set_vars(tsmear=tsmear)

    # Dataset 1
    #multi[0].set_vars(
    #)

    # Dataset 2
    multi[1].set_vars(
      nnsclo=1,
      nstep=1,
      nline=1,
      wfoptalg=1,
      useria=1,
      userib=0,
    )

    # Generate and return the input files.
    return multi.split_datasets()


def build_flow(structure, pseudos, ecuts, nksmall):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    #workdir = options.workdir
    #if not options.workdir:
    #    workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Get the SCF inputs.
    inputs = make_inputs(structure, pseudos, ecuts, nksmall)

    for inp in inputs:
     print(inp)

    # Build the flow.
    flow = abilab.Flow(workdir="hello_pert")
    work = abilab.Work()
    task0 = work.register_scf_task(inputs[0])
    task0.set_name("gs_step")
    task1 = work.register_scf_task(inputs[1], deps={task0: "WFK"})
    task1.set_name("gspert_step")
    flow.register_work(work)

    return flow


if __name__ == "__main__":
    pseudos = abidata.pseudos("14si.pspnc")
    structure=abidata.cif_file("si.cif")
    ecuts = [5, 10]
    nksmall = 4

    #inputs = make_inputs(structure, pseudos, ecuts, nksmall)
    flow = build_flow(structure, pseudos, ecuts, nksmall)
    flow.build_and_pickle_dump()
    #flow.make_scheduler().start()
