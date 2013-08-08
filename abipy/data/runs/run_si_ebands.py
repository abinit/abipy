#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import abipy.data as data  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.task import RunMode
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.io.abinitio.launcher import SimpleResourceManager
from pymatgen.io.abinitio.calculations import bandstructure

#class TestResults(object):
#    def __init__(self):

def find_ncdata(top):
    """
    Find all netcdf files starting from the top-level directory top

    Returns:
        dictionary with mapping: basename --> absolute path.
    """
    ncdata = {}
    for dirpath, dirnames, filenames in os.walk(top):
        for basename in filenames:
            apath = os.path.join(dirpath, basename)
            if basename.endswith(".nc"):
                assert basename not in ncdata
                ncdata[basename] = apath 

    return ncdata


def abidiff(ref_path, new_path):
    with open(ref_path, "r") as ref, open(new_path,"r") as new:
        ref_lines = ref.readlines()
        new_lines = new.readlines()
        return new_path
        #return ref_lines != new_lines


def main(dry_run=False):
    structure = AbiStructure.asabistructure(data.cif_file("si.cif"))

    pseudos = PseudoTable(data.pseudos("14si.pspnc"))
    runmode = RunMode.sequential()

    kppa = scf_kppa = 40
    nscf_nband = 6
    ndivsm = 5
    dos_ngkpt = [4,4,4]
    dos_shiftk = [0.1, 0.2, 0.3]
    max_ncpus = 1

    extra_abivars = dict(
        ecut=12, 
        accesswff=3, 
        istwfk="*1",
    )

    refdir = "data_" + os.path.basename(__file__).replace(".py","")
    workdir = "tmp_" + refdir

    refdir = os.path.join(os.path.dirname(__file__), refdir)
    workdir = os.path.join(os.path.dirname(__file__), workdir)

    work = bandstructure(workdir, runmode, structure, pseudos, scf_kppa, nscf_nband, ndivsm, 
                         spin_mode="unpolarized", smearing=None, **extra_abivars)

    retcodes = SimpleResourceManager(work, max_ncpus, sleep_time=5).run()
    retcode = max(retcodes)

    if retcode !=0:
        return retcode

    # Remove all files except those matching these regular expression.
    work[0].rename("out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    work[1].rename("out_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")

    work.rmtree(exclude_wildcard="*.abin|*_WFK*")
    work.rm_indatadir()
    work.rm_tmpdatadir()

    if not os.path.exists(refdir):
        work.move(refdir)

    else:
        diffs = {}
        for dirpath, dirnames, filenames in os.walk(refdir):
            for fname in filenames:
                ref_path = os.path.join(dirpath, fname)
                new_path = os.path.join(workdir, os.path.relpath(ref_path, start=refdir))
                diffs[ref_path] = abidiff(ref_path, new_path)

        print(diffs)

    #dos_kppa = 10
    #bands = bandstructure("hello_dos", runmode, structure, pseudos, scf_kppa, nscf_nband,
    #                      ndivsm, accuracy="normal", spin_mode="polarized",
    #                      smearing="fermi_dirac:0.1 eV", charge=0.0, scf_solver=None,
    #                      dos_kppa=dos_kppa)

    return retcode 


if __name__ == "__main__":
    sys.exit(main())
