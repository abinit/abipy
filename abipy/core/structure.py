# coding: utf-8
"""This module defines basic objects representing the crystalline structure."""
from __future__ import print_function, division, unicode_literals

import collections
import pymatgen
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.core.units import ArrayWithUnit
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from abipy.core.symmetries import SpaceGroup
from abipy.iotools import as_etsfreader, Visualizer
from abipy.iotools import xsf

__all__ = [
    "Lattice",
    "Structure",
]


class Lattice(pymatgen.Lattice):
    """
    Extends class:`pymatgen.Lattice` with methods that allows one
    to construct a Lattice object from ABINIT variables.
    """
    @classmethod
    def from_abivars(cls, d):
        """
        Returns a new instance from a dictionary with the variables 
        used in ABINIT to define the unit cell.
        """
        rprim = d.get("rprim", None)
        angdeg = d.get("angdeg", None)
        acell = d["acell"]

        # Call pymatgen constructors (note that pymatgen uses Angstrom instead of Bohr).
        if rprim is not None:
            assert angdeg is None
            rprim = np.reshape(rprim, (3,3))
            rprimd = [float(acell[i]) * rprim[i] for i in range(3)]
            return cls(ArrayWithUnit(rprimd, "bohr").to("ang"))

        elif angdeg is not None:
            # angdeg(0) is the angle between the 2nd and 3rd vectors,
            # angdeg(1) is the angle between the 1st and 3rd vectors,
            # angdeg(2) is the angle between the 1st and 2nd vectors,
            raise NotImplementedError("angdeg convention should be tested")
            angles = angdeg
            angles[1] = -angles[1]
            l = ArrayWithUnit(acell, "bohr").to("ang")
            return cls.from_lengths_and_angles(l, angdeg)

        else:
            raise ValueError("Don't know how to construct a Lattice from dict: %s" % str(d))

    #def to_abivars(self, **kwargs):
    #    # Should we use (rprim, acell) or (angdeg, acell) to specify the lattice?
    #    geomode = kwargs.pop("geomode", "rprim")
    #    if geomode == "rprim":
    #        return dict(acell=3 * [1.0], rprim=rprim))
    #    elif geomode == "angdeg":
    #        return dict(acell=3 * [1.0], angdeg=angdeg))
    #    else:
    #        raise ValueError("Wrong value for geomode: %s" % geomode)


class Structure(pymatgen.Structure):
    """
    Extends :class:`pymatgen.Structure` with methods that allows one
    to construct a Structure object from ABINIT variables.
    """
    @classmethod
    def from_file(cls, filepath):
        """
        Return a new Structure instance from a NetCDF file 

        Args:
            filename: netcdf file with crystallographic data in the ETSF-IO format.
                      or any other file format supported by `pymatgen.io.smartio`.
        """
        if filepath.endswith(".nc"):
            file, closeit = as_etsfreader(filepath)

            new = file.read_structure()
            # Change the class of new.
            new.__class__ = cls

            new.set_spacegroup(SpaceGroup.from_file(file))
            if closeit: file.close()

        else:
            # TODO: Spacegroup is missing here.
            new = super(Structure, cls).from_file(filepath)
            # Change the class of new.
            if new.__class__ != cls: new.__class__ = cls

        return new

    @classmethod
    def boxed_molecule(cls, pseudos, cart_coords, acell=3*(10,)):
        """
        Creates a molecule in a periodic box of lengths acell [Bohr]

        Args:
            pseudos: List of pseudopotentials
            cart_coords: Cartesian coordinates
            acell: Lengths of the box in *Bohr*
        """
        cart_coords = np.atleast_2d(cart_coords)

        molecule = Molecule([p.symbol for p in pseudos], cart_coords)

        l = ArrayWithUnit(acell, "bohr").to("ang")

        structure = molecule.get_boxed_structure(l[0], l[1], l[2])

        return cls(structure)

    @classmethod
    def boxed_atom(cls, pseudo, cart_coords=3*(0,), acell=3*(10,)):
        """
        Creates an atom in a periodic box of lengths acell [Bohr]

        Args:
            pseudo: Pseudopotential object.
            cart_coords: Cartesian coordinates
            acell: Lengths of the box in *Bohr*
        """
        return cls.boxed_molecule([pseudo], cart_coords, acell=acell)

    @classmethod
    def bcc(cls, a, species, **kwargs):
        """
        Build a primitive bcc crystal structure.

        Args:
            a: Lattice parameter in Angstrom.
            species: Chemical species. See __init__ method of :class:`pymatgen.Structue`
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structue`
        """
        lattice = 0.5 * float(a) * np.array([
            -1,  1,  1,
             1, -1,  1,
             1,  1, -1])

        return cls(lattice, species, coords=[[0, 0, 0]],  **kwargs)

    @classmethod
    def fcc(cls, a, species, **kwargs):
        """
        Build a primitive fcc crystal structure.

        Args:
            a: Lattice parameter in Angstrom.
            species: Chemical species. See __init__ method of :class:`pymatgen.Structure`
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structure`
        """
        # This is problematic
        lattice = 0.5 * float(a) * np.array([
            0,  1,  1,
            1,  0,  1,
            1,  1,  0])

        return cls(lattice, species, coords=[[0, 0, 0]], **kwargs)

    #@classmethod
    #def rocksalt(cls, a, sites, **kwargs):
    #    lattice = 0.5 * float(a) * np.array([
    #        0,  1,  1,
    #        1,  0,  1,
    #        1,  1,  0])
    #    coords = np.reshape([0, 0, 0, 0.5, 0.5, 0.5], (2,3))
    #    return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    #@classmethod
    #def ABO3(cls, a, sites, **kwargs)
    #   """Peroviskite structures."""
    #    return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    #@classmethod
    #def hH(cls, a, sites, **kwargs)
    #    return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    @property
    def spacegroup(self):
        """:class:`SpaceGroup` instance."""
        try:
            return self._spacegroup
        except AttributeError:
            return None

    def set_spacegroup(self, spacegroup):
        """`SpaceGroup` setter."""
        self._spacegroup = spacegroup

    @property
    def has_spacegroup(self):
        """True is the structure contains info on the spacegroup."""
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

    @lazy_property
    def hsym_kpath(self):
        """
        Returns an instance of :class:`HighSymmKpath`. (Database of high symmetry k-points and high symmetry lines).
        """
        from pymatgen.symmetry.bandstructure import HighSymmKpath
        return HighSymmKpath(self)

    @lazy_property
    def hsym_kpoints(self):
        """:class:`KpointList` object with the high-symmetry K-points."""
        # Get mapping name --> frac_coords for the special k-points in the database.
        name2frac_coords = self.hsym_kpath.kpath["kpoints"]
        kpath = self.hsym_kpath.kpath["path"]

        frac_coords, names = [], []
        for segment in kpath:
            for name in segment:
                fc = name2frac_coords[name]
                frac_coords.append(fc)
                names.append(name)

        # Build KpointList instance.
        from .kpoints import KpointList
        return KpointList(self.reciprocal_lattice, frac_coords, weights=None, names=names) 

    @lazy_property
    def hsym_stars(self):
        """
        List of :class:`Star` objects. Each star is associated to one of the special k-points
        present in the `pymatgen` database.
        """
        # Construct the stars.
        return [kpoint.compute_star(self.fm_symmops) for kpoint in self.hsym_kpoints]

    def get_sorted_structure_z(self):
        """Orders the structure according to increasing Z of the elements"""
        sites = sorted(self.sites, key=lambda site: site.specie.Z)
        return self.__class__.from_sites(sites)

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

        Returns: `matplotlib` figure.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        show              True to show the figure (Default).
        savefig           'abc.png' or 'abc.eps'* to save the figure to a file.
        ================  ==============================================================
        """
        return self.hsym_kpath.get_kpath_plot(**kwargs)

    def export(self, filename, visu=None):
        """
        Export the crystalline structure on file filename. 

        Args:
            filename: String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example,
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            visu: `Visualizer` subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See :class:`Visualizer` for the list of applications and formats supported.

        Returns: Instance of :class:`Visualizer`
        """
        if "." not in filename:
            raise ValueError("Cannot detect extension in filename %s: " % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: 
            # filename == ".ext" ==> Create temporary file.
            import tempfile
            filename = tempfile.mkstemp(suffix="." + ext, text=True)[1]

        with open(filename, mode="w") as fh:
            if ext == "xsf":  
                # xcrysden
                xsf.xsf_write_structure(fh, structures=[self])
            else:
                raise Visualizer.Error("extension %s is not supported." % ext)

        if visu is None:
            return Visualizer.from_file(filename)
        else:
            return visu(filename)

    def visualize(self, visu_name):
        """
        Visualize the crystalline structure with visualizer.
        See :class:`Visualizer` for the list of applications and formats supported.
        """
        # Get the Visualizer subclass from the string.
        visu = Visualizer.from_name(visu_name)

        # Try to export data to one of the formats supported by the visualizer
        # Use a temporary file (note "." + ext)
        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export(ext, visu=visu)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for %s" % visu_name)

    def to_abivars(self, **kwargs):
        """Returns a dictionary with the ABINIT variables."""
        types_of_specie = self.types_of_specie
        natom = self.num_sites

        znucl_type = [specie.number for specie in types_of_specie]

        znucl_atoms = self.atomic_numbers

        typat = np.zeros(natom, np.int)
        for (atm_idx, site) in enumerate(self):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

        rprim = ArrayWithUnit(self.lattice.matrix, "ang").to("bohr")
        xred = np.reshape([site.frac_coords for site in self], (-1,3))

        # Set small values to zero. This usually happens when the CIF file
        # does not give structure parameters with enough digits.
        #rprim = np.where(np.abs(rprim) > 1e-8, rprim, 0.0)
        #xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

        # Info on atoms.
        d = dict(
            natom=natom,
            ntypat=len(types_of_specie),
            typat=typat,
            znucl=znucl_type,
            xred=xred,
        )

        # Add info on the lattice.
        # Should we use (rprim, acell) or (angdeg, acell) to specify the lattice?
        geomode = kwargs.pop("geomode", "rprim")
        #latt_dict = self.lattice.to_abivars(geomode=geomode)

        if geomode == "rprim":
            d.update(dict(
                acell=3 * [1.0],
                rprim=rprim))

        elif geomode == "angdeg":
            d.update(dict(
                acell=3 * [1.0],
                angdeg=angdeg))
        else:
            raise ValueError("Wrong value for geomode: %s" % geomode)

        return d

    @classmethod
    def from_abivars(cls, d):
        """Build a `Structure` object from a dictionary containing ABINIT variables."""
        lattice = Lattice.from_abivars(d)

        coords, coords_are_cartesian = d.get("xred", None), False

        if coords is None:
            coords = d.get("xcart", None)
            if coords is not None:
                coords = ArrayWithUnit(coords, "bohr").to("ang")
            else:
                coords = d.get("xangst", None)
            coords_are_cartesian = True
        
        if coords is None:
            raise ValueError("Cannot extract atomic coordinates from dict %s" % str(d))

        coords = np.reshape(coords, (-1,3))

        znucl_type, typat = d["znucl"], d["typat"]

        if not isinstance(znucl_type, collections.Iterable):
            znucl_type = [znucl_type]

        if not isinstance(typat, collections.Iterable):
            typat = [typat]

        assert len(typat) == len(coords)

        # Note Fortan --> C indexing 
        species = [znucl_type[typ-1] for typ in typat]

        return cls(lattice, species, coords, validate_proximity=False,
                   to_unit_cell=False, coords_are_cartesian=coords_are_cartesian)

    def write_structure(self, filename):
        """See :ref:`pymatgen.io.smartio.write_structure`"""
        if filename.endswith(".nc"):
            raise NotImplementedError("Cannot write a structure to a netcdfile file yet")

        else:
            self.to(filename=filename)

    def convert(self, format="cif"):
        """
        Convert the Abinit structure to CIF, POSCAR, CSSR  and pymatgen's JSON serialized structures (json, mson)
        """
        prefix_dict = {
            "POSCAR": "POSCAR",
        }

        # FIXME:
        # Do we need symmetry operations here?
        # perhaps if the CIF file is used.
        suffix_dict = { 
            "cif": ".cif",
            "cssr": ".cssr",
            "json": ".json",
            "mson": ".mson",
        }

        if format not in prefix_dict.keys() and format not in suffix_dict.keys():
            raise ValueError("Unknown format %s" % format)

        prefix = prefix_dict.get(format, "tmp")
        suffix = suffix_dict.get(format, "")

        import tempfile
        tmp_file = tempfile.NamedTemporaryFile(suffix=suffix, prefix=prefix, mode="rw")

        self.write_structure(tmp_file.name)

        tmp_file.seek(0)
        return tmp_file.read()

    #def max_overlap_and_sites(self, pseudos):
    #    # For each site in self:
    #    # 1) Get the radius of the pseudopotential sphere 
    #    # 2) Get the neighbors of the site (considering the periodic images).

    #    max_overlap, ovlp_sites = 0.0, None

    #    for site in self:
    #        #site.specie
    #        #r = Length(pseudo.r_cut, "Bohr").to("ang")
    #        sitedist_list = self.get_neighbors(site, r, include_index=False)

    #        if sitedist_list:
    #            # Spheres are overlapping: compute overlap and update the return values 
    #            # if the new overlap is larger than the previous one.
    #            for (other_site, dist) in sitedist_list:
    #                # Eq 16 of http://mathworld.wolfram.com/Sphere-SphereIntersection.html
    #                #overlap = sphere_overlap(site.coords, r1, other_site.coords, r2)

    #                if overlap > max_overlap:
    #                    max_overlap = overlap
    #                    ovlp_sites = (site, other_site)

    #    return max_overlap, ovlp_sites

    def displace(self, displ, eta, frac_coords=True):
        """
        Displace the sites of the structure along the displacement vector displ.

        The displacement vector is first rescaled so that the maxium atomic displacement 
        is one Angstrom, and then multiplied by eta. Hence passing eta=0.001, will move 
        all the atoms so that the maximum atomic displacement is 0.001 Angstrom.

        Args:
            displ: Displacement vector with 3*len(self) entries (fractional coordinates).
            eta: Scaling factor.
            frac_coords: Boolean stating whether the vector corresponds to fractional or cartesian coordinates.
        """
        # Get a copy since we are going to modify displ.
        displ = np.reshape(displ, (-1,3)).copy()

        if len(displ) != len(self):
            raise ValueError("Displ must contains 3 * natom entries")

        if np.iscomplexobj(displ):
            raise TypeError("Displacement cannot be complex")

        if not frac_coords:
            # Convert to fractional coordinates.
            displ = np.reshape([self.lattice.get_fractional_coords(vec) for vec in displ], (-1,3))

        # Normalize the displacement so that the maximum atomic displacement is 1 Angstrom.
        dnorm = self.norm(displ, space="r")
        displ /= np.max(np.abs(dnorm))

        # Displace the sites.
        for i in range(len(self)):
            self.translate_sites(indices=i, vector=eta * displ[i, :], frac_coords=True)

    def get_smallest_supercell(self, qpoint, max_supercell):
        """
        Args:
            qpoint: q vector in reduced coordinate in reciprocal space
            max_supercell: vector with the maximum supercell size

        Returns: the scaling matrix of the supercell
        """
        if np.allclose(qpoint, 0):
            scale_matrix = np.eye(3, 3)
            return scale_matrix

        l = max_supercell

        # Inspired from Exciting Fortran code phcell.F90
        # It should be possible to improve this code taking advantage of python !
        scale_matrix = np.zeros((3,3),dtype=np.int)
        dmin = np.inf
        found = False

        # Try to reduce the matrix
        rprimd = self.lattice.matrix
        for l1 in np.arange(-l[0], l[0]+1):
            for l2 in np.arange(-l[1], l[1]+1):
                for l3 in np.arange(-l[2], l[2]+1):
                    lnew = np.array([l1, l2, l3])
                    ql = np.dot(lnew, qpoint)
                    # Check if integer and non zero !
                    if np.abs(ql - np.round(ql)) < 1e-6:
                        Rl = np.dot(lnew, rprimd)
                        # Normalize the displacement so that the maximum atomic displacement is 1 Angstrom.
                        dnorm = np.sqrt(np.dot(Rl,Rl))
                        if dnorm < dmin and dnorm > 1e-6:
                            found = True
                            scale_matrix[:, 0] = lnew
                            dmin = dnorm
        if not found:
            raise ValueError('max_supercell is not large enough for this q-point')

        found = False
        dmin = np.inf
        for l1 in np.arange(-l[0], l[0]+1):
            for l2 in np.arange(-l[1], l[1]+1):
                for l3 in np.arange(-l[2], l[2]+1):
                    lnew = np.array([l1, l2, l3])
                    # Check if not parallel !
                    cp = np.cross(lnew, scale_matrix[:,0])
                    if np.dot(cp,cp) > 1e-6:
                        ql = np.dot(lnew, qpoint)
                        # Check if integer and non zero !
                        if np.abs(ql - np.round(ql)) < 1e-6:
                            Rl = np.dot(lnew, rprimd)
                            dnorm = np.sqrt(np.dot(Rl, Rl))
                            if dnorm < dmin and dnorm > 1e-6:
                                found = True
                                scale_matrix[:, 1] = lnew
                                dmin = dnorm
        if not found:
            raise ValueError('max_supercell is not large enough for this q-point')

        dmin = np.inf
        found = False
        for l1 in np.arange(-l[0], l[0]+1):
            for l2 in np.arange(-l[1], l[1]+1):
                for l3 in np.arange(-l[2], l[2]+1):
                    lnew = np.array([l1, l2, l3])
                    # Check if not parallel !
                    cp = np.dot(np.cross(lnew, scale_matrix[:, 0]), scale_matrix[:, 1])
                    if cp > 1e-6:
                        # Should be positive as (R3 X R1).R2 > 0 for abinit !
                        ql = np.dot(lnew, qpoint)
                        # Check if integer and non zero !
                        if np.abs(ql - np.round(ql)) < 1e-6:
                            Rl = np.dot(lnew, rprimd)
                            dnorm = np.sqrt(np.dot(Rl,Rl))
                            if dnorm < dmin and dnorm > 1e-6:
                                found = True
                                scale_matrix[:, 2] = lnew
                                dmin = dnorm
        if not found:
            raise ValueError('max_supercell is not large enough for this q-point')

        # Fortran 2 python!!!
        return scale_matrix.T

    def get_trans_vect(self, scale_matrix):
        """
        Returns the translation vectors for a given scale matrix.

        Args:
            scale_matrix: Scale matrix defining the new lattice vectors in term of the old ones

        Return: the translation vectors
        """
        scale_matrix = np.array(scale_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)

        def range_vec(i):
            low = 0
            high = 0
            for z in scale_matrix[:, i]:
                if z > 0:
                    high += z
                else:
                    low += z
            return np.arange(low, high+1)
        arange = range_vec(0)[:, None] * np.array([1, 0, 0])[None, :]
        brange = range_vec(1)[:, None] * np.array([0, 1, 0])[None, :]
        crange = range_vec(2)[:, None] * np.array([0, 0, 1])[None, :]
        all_points = arange[:, None, None] + brange[None, :, None] +\
            crange[None, None, :]
        all_points = all_points.reshape((-1, 3))

        #find the translation vectors (in terms of the initial lattice vectors)
        #that are inside the unit cell defined by the scale matrix
        #we're using a slightly offset interval from 0 to 1 to avoid numerical
        #precision issues
        inv_matrix = np.linalg.inv(scale_matrix)


        frac_points = np.dot(all_points, inv_matrix)
        tvects = all_points[np.where(np.all(frac_points < 1-1e-10, axis=1)
                                     & np.all(frac_points >= -1e-10, axis=1))]
        assert len(tvects) == np.round(abs(np.linalg.det(scale_matrix)))

        return tvects

    def write_vib_file(self, xyz_file, qpoint, displ, do_real=True, frac_coords=True, scale_matrix=None, max_supercell=None):
        """
        write into the file descriptor xyz_file the positions and displacements of the atoms

        :param xyz_file: file_descriptor
        :param qpoint: qpoint to be analyzed
        :param displ: eigendisplacements to be analyzed
        :param do_real: True if you want to get only real part, False means imaginary part
        :param frac_coords: True if the eigendisplacements are given in fractional coordinates
        :param scale_matrix: Scale matrix for supercell
        :param max_supercell: Maximum size of supercell vectors with respect to primitive cell

        :return: nothing
        """
        if scale_matrix is None:
            if max_supercell is None:
                raise ValueError("If scale_matrix is not provided, please provide max_supercell !")

            scale_matrix = self.get_smallest_supercell(qpoint, max_supercell=max_supercell)

        old_lattice = self._lattice
        new_lattice = Lattice(np.dot(scale_matrix, old_lattice.matrix))

        tvects = self.get_trans_vect(scale_matrix)

        new_displ = np.zeros(3, dtype=np.float)

        fmtstr = "{{}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}}\n".format(6)

        for at,site in enumerate(self):
            for t in tvects:
                if do_real:
                    new_displ[:] = np.real(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ[at,:])
                else:
                    new_displ[:] = np.imag(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ[at,:])
                if frac_coords:
                    # Convert to fractional coordinates.
                    new_displ = self.lattice.get_cartesian_coords(new_displ)

                # We don't normalize here !!!
                fcoords = site.frac_coords + t

                coords = old_lattice.get_cartesian_coords(fcoords)

                new_fcoords = new_lattice.get_fractional_coords(coords)

                # New_fcoords -> map into 0 - 1
                new_fcoords = np.mod(new_fcoords, 1)
                coords = new_lattice.get_cartesian_coords(new_fcoords)

                xyz_file.write(fmtstr.format(site.specie, coords[0], coords[1], coords[2], new_displ[0], new_displ[1], new_displ[2]))

    def frozen_phonon(self, qpoint, displ, do_real=True, frac_coords=True, scale_matrix=None, max_supercell=None):
        """
        Compute the supercell needed for a given qpoint and add the displacement.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space.
            displ: displacement in real space of the atoms, will be normalized to 1 Angstrom.
            eta: pre-factor multiplying the displacement.
            do_real: true if we want only the real part of the displacement.
        """
        # I've copied code from make_supercell since the loop over supercell images
        # is inside make_supercell and I don't want to create a mapping

        if scale_matrix is None:
            if max_supercell is None:
                raise ValueError("If scale_matrix is not provided, please provide max_supercell !")

            scale_matrix = self.get_smallest_supercell(qpoint, max_supercell=max_supercell)

        scale_matrix = np.array(scale_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)

        old_lattice = self._lattice
        new_lattice = Lattice(np.dot(scale_matrix, old_lattice.matrix))

        tvects = self.get_trans_vect(scale_matrix)

        new_displ = np.zeros(3, dtype=np.float)
        new_sites = []
        for at,site in enumerate(self):
            for t in tvects:
                if(do_real):
                    new_displ[:] = np.real(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ[at,:])
                else:
                    new_displ[:] = np.imag(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ[at,:])
                if not frac_coords:
                    # Convert to fractional coordinates.
                    new_displ = self.lattice.get_fractional_coords(new_displ)

                # We don't normalize here !!!
                fcoords = site.frac_coords + t + new_displ
                coords = old_lattice.get_cartesian_coords(fcoords)
                new_site = PeriodicSite(
                    site.species_and_occu, coords, new_lattice,
                    coords_are_cartesian=True, properties=site.properties,
                    to_unit_cell=True)
                new_sites.append(new_site)

        self._sites = new_sites
        self._lattice = new_lattice

    def calc_kptbounds(self):
        """Returns the suggested value for the ABINIT variable `kptbounds`."""
        kptbounds = [k.frac_coords for k in self.hsym_kpoints]
        return np.reshape(kptbounds, (-1, 3))

    def calc_ksampling(self, nksmall, symprec=0.01, angle_tolerance=5):
        """
        Return the k-point sampling from the number of divisions to be used for
        the smallest lattive vectors of the reciprocal lattice.
        """
        ngkpt = self.calc_ngkpt(nksmall)
        shiftk = self.calc_shiftk(symprec=symprec, angle_tolerance=angle_tolerance)

        return AttrDict(
            ngkpt=ngkpt,
            shiftk=shiftk)

    def calc_ngkpt(self, nksmall): 
        """
        Compute the ABINIT variable `ngkpt` from the number of divisions used for the smallest lattice vector.

        Args:
            nksmall: Number of division for the smallest lattice vector.
        """
        lengths = self.lattice.reciprocal_lattice.abc
        lmin = np.min(lengths)

        ngkpt = np.ones(3, dtype=np.int)
        for i in range(3):
            ngkpt[i] = int(round(nksmall * lengths[i] / lmin))
            if ngkpt[i] == 0:
                ngkpt[i] = 1

        return ngkpt

    def calc_shiftk(self, symprec=0.01, angle_tolerance=5):
        """
        Find the values of `shiftk` and `nshiftk` appropriated for the sampling of the Brillouin zone.

        Returns
            Suggested value of shiftk

        .. note::

            When the primitive vectors of the lattice do NOT form a FCC or a BCC lattice, 
            the usual (shifted) Monkhorst-Pack grids are formed by using nshiftk=1 and shiftk 0.5 0.5 0.5 . 
            This is often the preferred k point sampling. For a non-shifted Monkhorst-Pack grid, 
            use nshiftk=1 and shiftk 0.0 0.0 0.0 , but there is little reason to do that.

            2) When the primitive vectors of the lattice form a FCC lattice, with rprim

                    0.0 0.5 0.5
                    0.5 0.0 0.5
                    0.5 0.5 0.0

            the (very efficient) usual Monkhorst-Pack sampling will be generated by using nshiftk= 4 and shiftk

                    0.5 0.5 0.5
                    0.5 0.0 0.0
                    0.0 0.5 0.0
                    0.0 0.0 0.5

            3) When the primitive vectors of the lattice form a BCC lattice, with rprim

                    -0.5  0.5  0.5
                    0.5 -0.5  0.5
                    0.5  0.5 -0.5

            the usual Monkhorst-Pack sampling will be generated by using nshiftk= 2 and shiftk

                    0.25  0.25  0.25
                    -0.25 -0.25 -0.25

            However, the simple sampling nshiftk=1 and shiftk 0.5 0.5 0.5 is excellent.

            4) For hexagonal lattices with hexagonal axes, e.g. rprim

                    1.0  0.0       0.0
                   -0.5  sqrt(3)/2 0.0
                    0.0  0.0       1.0

            one can use nshiftk= 1 and shiftk 0.0 0.0 0.5

            In rhombohedral axes, e.g. using angdeg 3*60., this corresponds to shiftk 0.5 0.5 0.5, 
            to keep the shift along the symmetry axis. 
        """
        # Find lattice type.
        sym = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        lattice_type = sym.get_lattice_type() 
        spg_symbol = sym.get_spacegroup_symbol()

        # Generate the appropriate set of shifts.
        shiftk = None

        if lattice_type == "cubic":
            if "F" in spg_symbol:  
                # FCC
                shiftk = [0.5, 0.5, 0.5,
                          0.5, 0.0, 0.0,
                          0.0, 0.5, 0.0,
                          0.0, 0.0, 0.5]

            elif "I" in spg_symbol:  
                # BCC
                shiftk = [0.25,  0.25,  0.25,
                         -0.25, -0.25, -0.25]

                #shiftk = [0.5, 0.5, 05])

        elif lattice_type == "hexagonal":
            # Find the hexagonal axis and set the shift along it.
            for i, angle in enumerate(self.lattice.angles):
                if abs(angle - 120) < 1.0:
                    j = (i + 1) % 3
                    k = (i + 2) % 3
                    hex_ax = [ax for ax in range(3) if ax not in [j,k]][0] 
                    break
            else:
                raise ValueError("Cannot find hexagonal axis")

            shiftk = [0.0, 0.0, 0.0]
            shiftk[hex_ax] = 0.5 

        if shiftk is None:
            # Use default value.
            shiftk = [0.5, 0.5, 0.5]

        return np.reshape(shiftk, (-1,3))

    def num_valence_electrons(self, pseudos):
        """
        Returns the number of valence electrons.

        Args:
            pseudos: List of :class:`Pseudo` objects or list of filenames.
        """
        table = PseudoTable.as_table(pseudos)

        nval = 0
        for site in self:
            symbol = site.species_string
            pseudos = table.pseudos_with_symbol(symbol)
            assert len(pseudos) == 1
            nval += pseudos[0].Z_val 

        return nval


class StructureModifier(object):
    """
    This object provides an easy-to-use interface for 
    generating new structures according to some algorithm.

    The main advantages of this approach are:
        
        *) Client code does not have to worry about the fact
           that many methods of Structure modify the object in place.

        *) One can render the interface more user-friendly. For example 
           some arguments might have a unit that can be specified in input.
           For example one can pass a length in Bohr that will be automatically 
           converted into Angstrom before calling the pymatgen methods
    """
    def __init__(self, structure):
        """
        Args:
            structure: Structure object.
        """
        # Get a copy to avoid any modification of the input. 
        self._original_structure = structure.copy()

    def copy_structure(self):
        """Returns a copy of the original structure."""
        return self._original_structure.copy()

    def scale_lattice(self, vol_ratios):
        """
        Scale the lattice vectors so that length proportions and angles are preserved.

        Args:
            vol_ratios: List with the ratios v/v0 where v0 is the volume of the original structure.

        Return: List of new structures with desired volume.
        """
        vol_ratios = np.array(vol_ratios)
        new_volumes = self._original_structure.volume * vol_ratios

        news = []
        for vol in new_volumes:
            new_structure = self.copy_structure()
            new_structure.scale_lattice(vol)
            news.append(new_structure)

        return news

    def make_supercell(self, scaling_matrix):
        """
        Create a supercell.

        Args:
            scaling_matrix:
                A scaling matrix for transforming the lattice vectors.
                Has to be all integers. Several options are possible:

                a. A full 3x3 scaling matrix defining the linear combination
                   the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
                   1]] generates a new structure with lattice vectors a' =
                   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
                   vectors of the original structure.
                b. An sequence of three scaling factors. E.g., [2, 1, 1]
                   specifies that the supercell should have dimensions 2a x b x
                   c.
                c. A number, which simply scales all lattice vectors by the
                   same factor.

        Returns:
            New structure.
        """
        new_structure = self.copy_structure()
        new_structure.make_supercell(scaling_matrix)
        return new_structure

    def displace(self, displ, etas, frac_coords=True):
        """
        Displace the sites of the structure along the displacement vector displ.

        The displacement vector is first rescaled so that the maxium atomic displacement 
        is one Angstrom, and then multiplied by eta. Hence passing eta=0.001, will move 
        all the atoms so that the maximum atomic displacement is 0.001 Angstrom.

        Args:
            displ:
                Displacement vector with 3*len(self) entries (fractional coordinates).
            eta:
                Scaling factor. 
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            List of new structures with displaced atoms.
        """
        if not isinstance(etas, collections.Iterable):
            etas = [etas]

        news = []
        for eta in etas:
            new_structure = self.copy_structure()
            new_structure.displace(displ, eta, frac_coords=frac_coords)
            news.append(new_structure)

        return news

    def frozen_phonon(self, qpoint, displ, do_real=True, frac_coords=True, scale_matrix=None, max_supercell=None):

        new_structure = self.copy_structure()
        new_structure.frozen_phonon(qpoint, displ, do_real, frac_coords, scale_matrix, max_supercell)

        return new_structure
