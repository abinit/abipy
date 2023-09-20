"""
"""
from __future__ import annotations

import os
import warnings
import numpy as np

from fnmatch import fnmatch
from monty.string import list_strings
from abipy.core.structure import Structure
from abipy.dynamics.hist import HistFile


def get_energy_step(step: dict) -> float:
    """
    Copied from final_energy property in vasp.outputs.

    Addresses a bug in vasprun.xml. See https://www.vasp.at/forum/viewtopic.php?f=3&t=16942
    """
    final_istep = step
    total_energy = final_istep["e_0_energy"]
    final_estep = final_istep["electronic_steps"][-1]
    electronic_energy_diff = final_estep["e_0_energy"] - final_estep["e_fr_energy"]
    total_energy_bugfix = np.round(electronic_energy_diff + final_istep["e_fr_energy"], 8)
    if np.abs(total_energy - total_energy_bugfix) > 1e-7:
        return total_energy_bugfix
    return total_energy


def get_structures_labels_from_file(filepath: str) -> tuple[list[Structure], dict]:
    """
    Read energies, forces, stresses and magmoms from an external file. 
    Return list of structures and dict with results.
    """
    basename = os.path.basename(str(filepath))

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
                #"magmoms": None,
            }

    elif fnmatch(basename, "vasprun*.xml*"):
        from pymatgen.io.vasp.outputs import Vasprun
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vasprun = Vasprun(filepath)

        structures = vasprun.structures
        natom = len(structures[0])
        energies_per_atom, forces, stresses, magmoms = [], [], [], []
        for step in vasprun.ionic_steps:
            energies_per_atom.append(get_energy_step(step) / natom)
            forces.append(step["forces"])
            stresses.append(step["stress"])

        labels = {
            "energies": energies_per_atom,
            "forces": forces,
            "stresses": stresses,
            #"magmoms": None,
        }

    else:
        raise ValueError(f"Don't know how to extract data from: {filepath=}")

    return structures, labels


def get_structures_labels_from_files(filepaths) -> tuple[list[Structure], dict]:
    """
    Read energies, forces, stresses and magmoms from an external file. 
    Return list of structures and dict with results.
    """
    for i, path in enumerate(list_strings(filepaths)):
        this_structures, this_labels = get_structures_labels_from_file(path)
        if i == 0:
            structures, labels = this_structures, this_labels
        else:
            structures += this_structures
            for k in labels:
                if labels[k] is None: continue
                labels[k] += this_labels[k]

    #for s in structures: print(s)
    #for k, v in labels.items(): print(k, v)

    return structures, labels

