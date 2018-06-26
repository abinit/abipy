from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import glob
import numpy as np
from collections import defaultdict
from abipy.core.func1d import Function1D
from abipy.electrons.gsr import GsrFile
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.abinit.pseudos import Pseudo
from monty.collections import tree
from monty.io import zopen
from monty.functools import lazy_property


class Coxp(object):
    """
    Wrapper class for the crystal orbital projections produced from Lobster.
    Wraps both a COOP and a COHP.
    Can contain both the total and orbitalwise projections.
    """

    def __init__(self, energies, total=None, partial=None, averaged=None, fermie=None):
        """
        Args:
            energies: list of energies. Shifted such that the Fermi level lies at 0 eV.
            total:  a dictionary with the values of the total overlap, not projected over orbitals.
                The dictionary should have the following nested keys: a tuple with the index of the sites
                considered as a pair (0-based, e.g. (0,1)), the spin (i.e. 0 or 1), a "single"
                or "integrated" string indicating the value corresponding to the value of
                the energy or to the integrated value up to that energy.
            partial: a dictionary with the values of the partial crystal orbital projections.
                The dictionary should have the following nested keys: a tuple with the index of the sites
                considered as a pair (0-based, e.g. (0,1)), a tuple with the string representing the
                projected orbitals for that pair (e.g. ("4s", "4p_x")), the spin (i.e. 0 or 1),
                a "single" or "integrated" string indicating the value corresponding to the value of
                the energy or to the integrated value up to that energy. Each dictionary should contain a
                numpy array with a list of COP values with the same size as energies.
            averaged: a dictionary with the values of the partial crystal orbital projections
                averaged over all atom pairs specified. The main key should indicate the spin (0 or 1)
                and the nested dictionary should have "single" and "integrated" as keys.
            fermie: value of the fermi energy.
        """
        self.energies = energies
        self.partial = partial or {}
        self.averaged = averaged or {}
        self.fermie = fermie
        self.total = total

    @property
    def site_pairs_total(self):
        """
        List of site pairs available for the total COP
        """
        return list(self.total.keys())

    @property
    def site_pairs_partial(self):
        """
        List of site pairs available for the partial COP
        """
        return list(self.partial.keys())

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance of Coxp from the files produce by Lobster.
        Accepts gzipped files.

        Args:
            filepath: path to the COHPCAR.lobster or COOPCAR.lobster.

        Returns:
            A Coxp.
        """

        float_patt = r'-?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
        header_patt = re.compile(r'\s+(\d+)\s+(\d+)\s+(\d+)\s+('+float_patt+
                                 r')\s+('+float_patt+r')\s+('+float_patt+r')')
        pair_patt = re.compile(r'No\.\d+:([a-zA-Z]+)(\d+)(?:\[([a-z0-9_\-^]+)\])?->([a-zA-Z]+)(\d+)(?:\[([a-z0-9_\-^]+)\])?')

        with zopen(filepath, "rt") as f:
            #find the header
            for line in f:
                match = header_patt.match(line.rstrip())
                if match:
                    n_column_groups = int(match.group(1))
                    n_spin = int(match.group(2))
                    n_en_steps = int(match.group(3))
                    min_en = float(match.group(4))
                    max_en = float(match.group(5))
                    fermie = float(match.group(6))
                    break
            else:
                raise ValueError("Can't find the header in file {}".format(filepath))

            n_pairs = n_column_groups-1

            count_pairs = 0
            pairs_data = []
            #parse the pairs considered
            for line in f:
                match = pair_patt.match(line.rstrip())
                if match:
                    # adds a tuple: [type1, index1, orbital1, type2, index2, orbital2]
                    # with orbital1, orbital2 = None if the pair is not orbitalwise
                    type1, index1, orbital1, type2, index2, orbital2 = match.groups()
                    #0-based indexing
                    index1 = int(index1) - 1
                    index2 = int(index2) - 1
                    pairs_data.append([type1, index1, orbital1, type2, index2, orbital2])
                    count_pairs += 1
                    if count_pairs == n_pairs:
                        break

            spins = [0, 1][:n_spin]

            data = np.fromstring(f.read(), dtype=np.float, sep=' ').reshape([n_en_steps, 1+n_spin*n_column_groups*2])

            # initialize and fill results
            energies = data[:, 0]
            averaged = defaultdict(dict)
            total = tree()
            partial = tree()

            for i, s in enumerate(spins):
                base_index = 1+i*n_column_groups*2
                averaged[s]['single'] = data[:, base_index]
                averaged[s]['integrated'] = data[:, base_index+1]
                for j, p in enumerate(pairs_data):
                    index1 = p[1]
                    index2 = p[4]
                    # partial or total
                    if p[2] is not None:
                        single = data[:, base_index+2*(j+1)]
                        integrated = data[:, base_index+2*(j+1)+1]
                        partial[(index1, index2)][(p[2], p[5])][s]['single'] = single
                        partial[(index2, index1)][(p[5], p[2])][s]['single'] = single
                        partial[(index1, index2)][(p[2], p[5])][s]['integrated'] = integrated
                        partial[(index2, index1)][(p[5], p[2])][s]['integrated'] = integrated
                    else:
                        single = data[:, base_index+2*(j+1)]
                        integrated = data[:, base_index+2*(j+1)+1]
                        total[(index1, index2)][s]['single'] = single
                        total[(index2, index1)][s]['single'] = single
                        total[(index1, index2)][s]['integrated'] = integrated
                        total[(index2, index1)][s]['integrated'] = integrated

        return cls(energies=energies, total=total, partial=partial, averaged=averaged, fermie=fermie)

    @lazy_property
    def functions_pair_lorbitals(self):
        """
        Extracts a dictionary with keys pair, orbital, spin and containing a Function1D object resolved
        for l orbitals.
        """
        if not self.partial:
            raise RuntimeError("Partial orbitals not calculated.")

        results = tree()
        for pair, pair_data in self.partial.items():
            #check if the symmetric has already been calculated
            if (pair[1], pair[0]) in results:
                for orbs, orbs_data in results[(pair[1], pair[0])].items():
                    results[pair][(orbs[1], orbs[0])] = orbs_data
                continue

            #for each look at all orbital possibilities
            for orbs, orbs_data in pair_data.items():
                k=(orbs[0].split("_")[0], orbs[1].split("_")[0])
                if k in results[pair]:
                    for spin in orbs_data.keys():
                        results[pair][k][spin]=results[pair][k][spin]+Function1D(self.energies,orbs_data[spin]['single'])
                else:
                    for spin in orbs_data.keys():
                        results[pair][k][spin]=Function1D(self.energies, orbs_data[spin]['single'])

        return results

    @lazy_property
    def functions_pair_morbitals(self):
        """
        Extracts a dictionary with keys pair, orbital, spin and containing a Function1D object resolved
        for l and m orbitals.
        """

        if not self.partial:
            raise RuntimeError("Partial orbitals not calculated.")
        results = tree()
        for pair, pair_data in self.partial.items():
            for orbs, orbs_data in pair_data.items():
                for spin in orbs_data.keys():
                    results[pair][orbs][spin]=Function1D(self.energies, orbs_data[spin]['single'])
        return results

    @lazy_property
    def functions_pair(self):
        """
        Extracts a dictionary with keys pair, spin and containing a Function1D object for the total COP.
        """

        results = tree()
        for pair, pair_data in self.total.items():
            for spin in pair_data.keys():
                results[pair][spin] = Function1D(self.energies, pair_data[spin]['single'])
        return results


class ICoxp(object):
    """
    Wrapper class for the integrated crystal orbital projections up to the Fermi energy
    produced from Lobster.
    May contain the output stored in ICOHPLIST.lobster and ICOOPLIST.lobster
    """

    def __init__(self, values):
        """
        Args:
            values: a dictionary with the following keys: a tuple with the index of the sites
                considered as a pair (0-based, e.g. (0,1)), the spin (i.e. 0 or 1)
        """
        self.values = values

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance of ICoxp from the files produce by Lobster.
        Accepts gzipped files.

        Args:
            filepath: path to the ICOHPLIST.lobster or ICOOPLIST.lobster.

        Returns:
            A ICoxp.
        """
        float_patt = r'-?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
        header_patt = re.compile(r'.*?(over+\s#\s+bonds)?\s+for\s+spin\s+(\d).*')
        data_patt = re.compile(r'\s+\d+\s+([a-zA-Z]+)(\d+)\s+([a-zA-Z]+)(\d+)\s+('+
                               float_patt+r')\s+('+float_patt+r')(\d+)?')
        values = tree()
        spin = None
        avg_num_bonds = False
        with zopen(filepath, "rt") as f:
            for line in f:
                match = header_patt.match(line.rstrip())
                if match:
                    spin = [0, 1][int(match.group(2))-1]
                    avg_num_bonds = match.group(1) is not None
                match = data_patt.match(line.rstrip())
                if match:
                    type1, index1, type2, index2, dist, avg, n_bonds = match.groups()
                    #0-based indexing
                    index1 = int(index1) - 1
                    index2 = int(index2) - 1
                    avg_data = {'average': float(avg), 'distance': dist, 'n_bonds': int(n_bonds) if n_bonds else None}
                    values[(index1, index2)][spin] = avg_data
                    values[(index2, index1)][spin] = avg_data

        return cls(values)


class LobsterDos(object):
    """
    Total and partial dos extracted from lobster DOSCAR.
    The fermi energy is always at the zero value.
    """

    def __init__(self, energies, total_dos, pdos):
        """

        Args:
            energies: list of energies. Shifted such that the Fermi level lies at 0 eV.
            total_dos: a dictionary with spin as a key (i.e. 0 or 1) and containing the values of the total DOS.
                Should have the same size as energies.
            pdos: a dictionary with the values of the projected DOS.
                The dictionary should have the following nested keys: the index of the site (0-based),
                the string representing the projected orbital (e.g. "4p_x"), the spin (i.e. 0 or 1).
                Each dictionary should contain a numpy array with a list of DOS values with the
                same size as energies.
        """

        self.energies = energies
        self.total_dos = total_dos
        self.pdos = pdos

    @classmethod
    def from_file(self, filepath):
        """
        Generates an instance from the DOSCAR.lobster file.
        Accepts gzipped files.

        Args:
            filepath: path to the DOSCAR.lobster.

        Returns:
            A LobsterDos.
        """

        with zopen(filepath, "rt") as f:
            dos_data = f.readlines()

        n_sites = int(dos_data[0].split()[0])
        n_energies = int(dos_data[5].split()[2])

        n_spin = 1 if len(dos_data[6].split()) == 3 else 2
        spins = [0, 1][:n_spin]

        # extract np array for total dos
        tdos_data = np.fromiter((d for l in dos_data[6:6+n_energies] for d in l.split()),
                                dtype=np.float).reshape((n_energies, 1+2*n_spin))

        energies = tdos_data[:,0]
        total_dos = {}
        for i_spin, spin in enumerate(spins):
            total_dos[spin] = tdos_data[:,1+2*i_spin]

        pdos = tree()
        # read partial doses
        for i_site in range(n_sites):
            i_first_line = 5+(n_energies+1)*(i_site+1)

            # read orbitals
            orbitals = dos_data[i_first_line].rsplit(';', 1)[-1].split()

            # extract np array for partial dos
            pdos_data = np.fromiter((d for l in dos_data[i_first_line+1:i_first_line+1+n_energies] for d in l.split()),
                dtype=np.float).reshape((n_energies, 1+n_spin*len(orbitals)))

            for i_orb, orb in enumerate(orbitals):
                for i_spin, spin in enumerate(spins):
                    pdos[i_site][orb][spin] = pdos_data[:, i_spin+n_spin*i_orb+1]

        return LobsterDos(energies=energies, total_dos=total_dos, pdos=pdos)


class LobsterInput(object):
    """
    This object stores the basic variables for a Lobster input and generates the lobsterin file.
    """

    accepted_basis_sets = {"bunge", "koga", "pbevaspfit2015"}

    available_advanced_options = {"basisRotation", "writeBasisFunctions", "onlyReadVasprun.xml", "noMemoryMappedFiles",
                                  "skipPAWOrthonormalityTest", "doNotIgnoreExcessiveBands", "doNotUseAbsoluteSpilling",
                                  "skipReOrthonormalization", "doNotOrthogonalizeBasis", "forceV1HMatrix",
                                  "noSymmetryCorrection", "symmetryDetectionPrecision", "useOriginalTetrahedronMethod",
                                  "useDecimalPlaces", "forceEnergyRange"}

    def __init__(self, basis_set=None, basis_functions=None, include_orbitals=None, atom_pairs=None, dist_range=None,
                 orbitalwise=True, start_en=None, end_en=None, en_steps=None, gaussian_smearing=None,
                 bwdf=None, advanced_options=None):
        """
        Args
            basis_set: String containing one of the possible basis sets available: bunge, koga, pbevaspfit2015
            basis_functions: list of strings giving the symbol of each atom and the basis functions: "Ga 4s 4p"
            include_orbitals: string containing which types of valence orbitals to use. E.g. "s p d"
            atom_pairs: list of tuples containing the couple of elements for which the COHP analysis will be
             performed. Index is 1-based.
            dist_range: list of tuples, each containing the minimum and maximum distance (in Angstrom) used to
             automatically generate atom pairs. Each tuple can also contain two atomic symbol to restric the match to
             the specified elements. examples: (0.5, 1.5) or (0.5, 1.5, 'Zn', 'O')
            start_en: starting energy with respect to the Fermi level (in eV)
            end_en: ending energy with respect to the Fermi level (in eV)
            en_steps: number of energy increments
            gaussian_smearing: smearing in eV when using gaussian broadening
            bwdf: enables the bond-weighted distribution function (BWDF). Value is the binning interval
            advanced_options: dict with additional advanced options. See lobster user guide for further details
        """

        if basis_set and basis_set.lower() not in self.accepted_basis_sets:
            raise ValueError("Wrong basis set {}".format(basis_set))
        self.basis_set = basis_set
        self.basis_functions = basis_functions or []
        self.include_orbitals = include_orbitals
        self.atom_pairs = atom_pairs or []
        self.dist_range = dist_range or []
        self.orbitalwise = orbitalwise
        self.start_en = start_en
        self.end_en = end_en
        self.en_steps = en_steps
        self.gaussian_smearing = gaussian_smearing
        self.bwdf = bwdf
        self.advanced_options = advanced_options or {}

        if not all(opt in self.available_advanced_options for opt in self.advanced_options.keys()):
            raise ValueError("Unknown adavanced options")

    @classmethod
    def _get_basis_functions_from_abinit_pseudos(cls, pseudos):
        """
        Extracts the basis function used from the PAW abinit pseudopotentials

        Args:
            pseudos: a list of Pseudos objects.
        """

        basis_functions = []
        for p in pseudos:
            if not hasattr(p, "valence_states"):
                raise RuntimeError("Only PAW pseudos in PAWXML format are supported by Lobster interface.")
            el = p.symbol + " ".join(str(vs['n'] + OrbitalType(int(vs['l'])).name)
                                     for vs in p.valence_states.values() if 'n' in vs)
            basis_functions.append(el)
        return basis_functions

    def set_basis_functions_from_abinit_pseudos(self, pseudos):
        """
        Sets the basis function used from the PAW abinit pseudopotentials

        Args:
            pseudos: a list of Pseudos objects.
        """
        basis_functions = self._get_basis_functions_from_abinit_pseudos(pseudos)

        self.basis_functions = basis_functions

    @classmethod
    def _get_basis_functions_from_potcar(cls, potcar):
        """
        Extracts the basis function used from a POTCAR.

        Args:
            potcar: a pymatgen.io.vasp.inputs.Potcar object
        """
        basis_functions = []
        for p in potcar:
            basis_functions.append(p.element +" "+ " ".join(str(vs[0]) + vs[1] for vs in p.electron_configuration))
        return basis_functions

    def set_basis_functions_from_potcar(self, potcar):
        """
        Sets the basis function used from a POTCAR.

        Args:
            potcar: a pymatgen.io.vasp.inputs.Potcar object
        """
        basis_functions = self._get_basis_functions_from_potcar(potcar)

        self.basis_functions = basis_functions

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """
        String representation.
        """
        lines = []

        if self.basis_set:
            lines.append("basisSet {}".format(self.basis_set))

        for bf in self.basis_functions:
            lines.append("basisFunctions {}".format(bf))

        for ap in self.atom_pairs:
            line = "cohpBetween atom {} atom {}".format(*ap)
            if self.orbitalwise:
                line += " orbitalwise"
            lines.append(line)

        for dr in self.dist_range:
            line = "cohpGenerator from {} to {}".format(dr[0], dr[1])
            if len(dr) > 2:
                line += " type {} type {}".format(dr[2], dr[3])
            if self.orbitalwise:
                line += " orbitalwise"
            lines.append(line)

        if self.start_en:
            lines.append("COHPStartEnergy {}".format(self.start_en))

        if self.end_en:
            lines.append("COHPEndEnergy {}".format(self.end_en))

        if self.en_steps:
            lines.append("COHPSteps {}".format(self.en_steps))

        if self.gaussian_smearing:
            lines.append("gaussianSmearingWidth {}".format(self.gaussing_smearing))

        if self.bwdf:
            lines.append("BWDF {}".format(self.bwdf))

        for k, v in self.advanced_options.items():
            lines.append("{} {}".format(k, v))

        return "\n".join(lines)

    @classmethod
    def from_run_dir(cls, dirpath, dE=0.01, **kwargs):
        """
        Generates an instance of the class based on the output folder of a DFT calculation.
        Reads the information from the pseudopotentials in order to determine the
        basis functions.

        Args:
            dirpath: the path to the calculation directory. For abinit it should contain the
                "files" file, for vasp it should contain the vasprun.xml and the POTCAR.
            dE: The spacing of the energy sampling in eV.
            kwargs: the inputs for the init method, except for basis_functions, start_en,
                end_en and en_steps.

        Returns:
            A LobsterInput.
        """

        # Try to determine the code used for the calculation
        dft_code = None
        if os.path.isfile(os.path.join(dirpath, 'vasprun.xml')):
            dft_code = "vasp"
            vr = Vasprun(os.path.join(dirpath, 'vasprun.xml'))

            en_min = np.min([bands_spin for bands_spin in vr.eigenvalues.values()])
            en_max = np.max([bands_spin for bands_spin in vr.eigenvalues.values()])
            fermie = vr.efermi

            potcar = Potcar.from_file(os.path.join(dirpath, 'POTCAR'))
            basis_functions = cls._get_basis_functions_from_potcar(potcar)

        elif glob.glob(os.path.join(dirpath, '*.files')):
            dft_code = "abinit"
            ff = glob.glob(os.path.join(dirpath, '*.files'))[0]
            with open(ff, "rt") as files_file:
                ff_lines = files_file.readlines()
            out_path = ff_lines[3].strip()
            if not os.path.isabs(out_path):
                out_path = os.path.join(dirpath, out_path)

            with GsrFile.from_file(out_path+'_GSR.nc') as gsr:
                en_min = gsr.ebands.eigens.min()
                en_max = gsr.ebands.eigens.max()
                fermie = gsr.ebands.fermie

            pseudo_paths = []
            for l in ff_lines[5:]:
                l = l.strip()
                if l:
                    if not os.path.isabs(l):
                        l = os.path.join(dirpath, l)
                    pseudo_paths.append(l)

            pseudos = [Pseudo.from_file(p) for p in pseudo_paths]

            basis_functions = cls._get_basis_functions_from_abinit_pseudos(pseudos)
        else:
            raise ValueError('Unable to determine the code used in dir {}'.format(dirpath))

        start_en = en_min + fermie
        end_en = en_max - fermie

        # shift the energies so that are divisible by dE and the value for the fermi level (0 eV) is included
        start_en = np.floor(start_en/dE)*dE
        end_en= np.ceil(end_en/dE)*dE

        en_steps = int((end_en-start_en)/dE)

        return cls(basis_functions=basis_functions, start_en=start_en, end_en=end_en, en_steps=en_steps, **kwargs)


    def write(self, dirpath='.'):
        """
        Write the input file 'lobsterin' in dirpath.
        """

        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        # Write the input file.
        with open(os.path.join(dirpath, 'lobsterin'), "wt") as f:
            f.write(str(self))
