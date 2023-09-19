"""
"""
from __future__ import annotations

import os

from fnmatch import fnmatch
from monty.string import list_strings
from abipy.dynamics.hist import HistFile


def get_structures_labels_from_file(filepath):
    """
    """
    basename = os.path.basename(filepath)


    if basename.endswith("_HIST.nc"):
        with HistFile(filepath) as hist:
            structures = hist.reader.read_all_structures()
            # etotals per atom in eV.
            energies_per_atom = [float(e) / hist.r.natom for e in hist.etotals]
            cart_forces = hist.r.read_cart_forces(unit="eV ang^-1")
            # GPa units.
            stress_cart_tensors, pressures = hist.r.read_cart_stress_tensors()
            labels = {
                "energies": energies_per_atom,
                "forces": cart_forces,
                "stresses": stress_cart_tensors,
            }
            return structures, labels

    elif fnmatch(basename, "vasprun*.xml*"):
        #from pymatgen.io.vasp.outputs import Vasprun
        #with warnings.catch_warnings():
        #    warnings.simplefilter("ignore")
        #    vasprun = Vasprun(filepath)

        #num_steps = len(vasprun.ionic_steps)
        #for istep, step in enumerate(vasprun.ionic_steps):
        #    #print(step.keys())
        #    structure, forces, stress = step["structure"], step["forces"], step["stress"]
        #    ene = get_energy_step(step)
        raise NotImplementedError
        return structures, labels

    raise ValueError(f"Don't know how to extract data from: {filepath=}")


def get_structures_labels_from_files(filepaths):
    """
    """
    filepaths = list_strings(filepaths)
    for i, path in enumerate(filepaths):
        this_structures, this_labels = get_structures_labels_from_file(path)
        if i == 0:
            structures, labels = this_structures, this_labels
        else:
            structures += this_structures
            for k in labels:
                labels[k] += this_labels[key]

    return structures, labels

