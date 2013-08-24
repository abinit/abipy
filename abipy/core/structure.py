"""This module defines the structure object that store information on the crystalline structure and its symmetries."""
from __future__ import division, print_function

import numpy as np

from .constants import Bohr2Ang, Ang2Bohr
from .symmetries import SpaceGroup
from abipy.iotools import as_etsfreader, Visualizer
from abipy.iotools import xsf

import pymatgen

__all__ = [
    "Structure",
    "Lattice",
]


class Lattice(pymatgen.Lattice):

    @classmethod
    def from_abivars(cls, d):
        rprim = d.get("rprim", None)
        angdeg = d.get("angdeg", None)
        acell = d["acell"]

        # Call pymatgen constructors (note that pymatgen uses Angstrom instead of Bohr).
        if rprim is not None:
            rprim = np.reshape(rprim, (3,3))
            rprimd = [acell[i] * rprim[i] for i in range(3)]
            return cls(Bohr2Ang(rprimd))

        elif angdeg is not None:
            # angdeg(0) is the angle between the 2nd and 3rd vectors,
            # angdeg(1) is the angle between the 1st and 3rd vectors,
            # angdeg(2) is the angle between the 1st and 2nd vectors,
            raise NotImplementedError("angdeg convention should be tested")
            angles = angdeg
            angles[1] = -angles[1]
            new = cls.from_lengths_and_angles(Bohr2Ang(acell), angdeg)
            new.__class__ = cls
            return new

        else:
            raise ValueError("Don't know how to construct a Lattice from dict: %s" % d)


class Structure(pymatgen.Structure):

    @classmethod
    def from_file(cls, filepath):
        """
        Return a new instance from a NetCDF file containing crystallographic data in the ETSF-IO format.

        Args:
            ncdata:
                filename or NetcdfReader instance.
        """
        if filepath.endswith(".nc"):
            file, closeit = as_etsfreader(filepath)

            new = file.read_structure()
            # Change the class of new.
            new.__class__ = cls

            new.set_spacegroup(SpaceGroup.from_file(file))

            if closeit:
                file.close()
        else:
            # TODO: Spacegroup is missing here.
            from pymatgen.io.smartio import read_structure
            new = read_structure(filepath)
            # Change the class of new.
            new.__class__ = cls

        return new

    @property
    def spacegroup(self):
        """`SpaceGroup` instance."""
        try:
            return self._spacegroup
        except AttributeError:
            return None

    def set_spacegroup(self, spacegroup):
        """`SpaceGroup` setter."""
        self._spacegroup = spacegroup

    @property
    def has_spacegroup(self):
        """True is self contains info on the spacegroup."""
        return self.spacegroup is not None

    @property
    def is_symmorphic(self):
        """True if at least one fractional translation is non-zero."""
        return self.spacegroup.is_symmorphic

    @property
    def fm_symmops(self):
        """Tuple with ferromagnetic symmetries (time-reversal is included, if present)."""
        return self.spacegroup.symmops(afm_sign=+1)

    @property
    def afm_symmops(self):
        """Tuple with Anti-ferromagnetic symmetries (time-reversal is included, if present)."""
        return self.spacegroup.symmops(afm_sign=-1)

    @property
    def hsym_kpath(self):
        """
        Returns an instance of the pymatgen class `HighSymmKpath`
        (Database of high symmetry k-points and high symmetry lines).
        """
        try:
            return self._hsym_kpath

        except AttributeError:
            from pymatgen.symmetry.bandstructure import HighSymmKpath
            self._hsym_kpath = HighSymmKpath(self)
            return self._hsym_kpath

    @property
    def hsym_kpoints(self):
        """`KpointList` object with the high-symmetry K-points."""
        try:
            return self._hsym_kpoints

        except AttributeError:
            # Get mapping name --> frac_coords for the special k-points in the database.
            name2frac_coords = self.hsym_kpath.kpath["kpoints"]

            frac_coords, names = [], []
            for (name, fc) in name2frac_coords.items():
                frac_coords.append(fc)
                names.append(name)

            # Build KpointList instance.
            from .kpoints import KpointList
            self._hsym_kpoints = KpointList(self.reciprocal_lattice, frac_coords, weights=None, names=names) 
            return self._hsym_kpoints

    @property
    def hsym_stars(self):
        """
        List of `Star` objects. Each start is associated to one of the special k-points 
        present in the pymatgen database.
        """
        try:
            return self._hsym_stars

        except AttributeError:
            # Construct the stars.
            self._hsym_stars = [kpoint.compute_star(self.fm_symmops) for kpoint in self.hsym_kpoints]
            return self._hsym_stars

    def findname_in_hsym_stars(self, kpoint):
        """Returns the name of the special k-point, None if kpoint is unknown.""" 
        for star in self.hsym_stars:
            if star.find(kpoint) != -1:
                return star.name
        else:
            return None

    def show_bz(self, **kwargs):
        """
        Gives the plot (as a matplotlib object) of the symmetry line path in the Brillouin Zone.

        Returns:
            `matplotlib` figure.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        show              True to show the figure (Default).
        savefig           'abc.png' or 'abc.eps'* to save the figure to a file.
        ================  ==============================================================
        """
        return self.hsym_kpath.get_kpath_plot(**kwargs)

    def export(self, filename):
        """
        Export the crystalline structure on file filename.

        Returns:
            Instance of :class:`Visualizer`

        The format is defined by the extension in filename:
        See :class:`Visualizer` for the list of applications and formats supported.

            #. "prefix.xsf" for XcrysDen files.

        An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
        """
        if "." not in filename:
            raise ValueError("Cannot detect extension in filename %s: " % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: 
            # filename == ".ext" ==> Create temporary file.
            import tempfile
            filename = tempfile.mkstemp(suffix="."+ext, text=True)[1]

        with open(filename, mode="w") as fh:
            if ext == "xsf": # xcrysden
                xsf.xsf_write_structure(fh, structures=[self])
            else:
                raise Visualizer.Error("extension %s is not supported." % ext)

        return Visualizer.from_file(filename)

    def visualize(self, visualizer):
        """
        Visualize the crystalline structure with visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        extensions = Visualizer.exts_from_appname(visualizer)

        for ext in extensions:
            ext = "." + ext
            try:
                return self.export(ext)
            except Visualizer.Error:
                pass
        else:
            raise Visualizer.Error("Don't know how to export data for %s" % visualizer)

    def to_abivars(self):
        """Returns a dictionary with the abinit variables."""
        types_of_specie = self.types_of_specie
        natom = self.num_sites

        znucl_type = [specie.number for specie in types_of_specie]

        znucl_atoms = self.atomic_numbers

        typat = np.zeros(natom, np.int)
        for (atm_idx, site) in enumerate(self):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

        #significant_figures = 12
        #format_str = "{{:.{0}f}}".format(significant_figures)
        #fmt = format_str.format

        #lines = []
        #for vec in Ang2Bohr(self.lattice.matrix):
        #    lines.append(" ".join([fmt(c) for c in vec]))
        #rprim = "\n" + "\n".join(lines)

        #lines = []
        #for (i, site) in enumerate(self):
        #    coords = site.frac_coords
        #    lines.append( " ".join([fmt(c) for c in coords]) + " # " + site.species_string )
        #xred = '\n' + "\n".join(lines)

        rprim = Ang2Bohr(self.lattice.matrix)
        xred = np.reshape([site.frac_coords for site in self], (-1,3))

        return {
            "acell" : 3 * [1.0],
            "rprim" : rprim,
            "natom" : natom,
            "ntypat": len(types_of_specie),
            "typat" : typat,
            "znucl" : znucl_type,
            "xred"  : xred,
        }

    @classmethod
    def from_abivars(cls, d):
        lattice = Lattice.from_abivars(d)

        coords, coords_are_cartesian = d.get("xred", None), False

        if coords is None:
            coords = d.get("xcart", None)
            if coords is not None:
                coords = Bohr2Ang(coords)
            else:
                coords = d.get("xangst", None)
            coords_are_cartesian = True
        
        if coords is None:
            raise ValueError("Cannot extract atomic coordinates from dict %s" % str(d))

        coords = np.reshape(coords, (-1,3))

        znucl_type, typat = d["znucl"], d["typat"]
        assert len(typat) == len(coords)

        # Note Fortan --> C indexing 
        species = [znucl_type[typ-1] for typ in typat]

        return cls(lattice, species, coords, validate_proximity=False,
                   to_unit_cell=False, coords_are_cartesian=coords_are_cartesian)

    def write_structure(self, filename):
        """See `pymatgen.io.smartio.write_structure`"""
        from pymatgen.io.smartio import write_structure
        write_structure(self, filename)

