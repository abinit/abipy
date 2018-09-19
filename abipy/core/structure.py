# coding: utf-8
"""
This module defines basic objects representing the crystalline structure.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import collections
import tempfile
import numpy as np
import pickle
import pymatgen
import pymatgen.core.units as pmg_units

from pprint import pprint, pformat
from warnings import warn
from collections import OrderedDict
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.string import is_string, marquee, list_strings
from monty.termcolor import cprint
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.flowtk import PseudoTable
from abipy.core.mixins import NotebookWriter
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.iotools import as_etsfreader, Visualizer, xsf
from abipy.flowtk.abiobjects import structure_from_abivars, structure_to_abivars


__all__ = [
    "mp_match_structure",
    "mp_search",
    "cod_search",
    "Structure",
    "dataframes_from_structures",
]


def mp_match_structure(obj, api_key=None, endpoint=None, final=True):
    """
    Finds matching structures on the Materials Project database.

    Args:
        obj: filename or |Structure| object.
        api_key (str): A String API key for accessing the MaterialsProject REST interface.
        endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
        final (bool): Whether to get the final structure, or the initial
            (pre-relaxation) structure. Defaults to True.

    Returns:
        :class:`MpStructures` object with
            structures: List of matching structures and list of Materials Project identifier.
    """
    structure = Structure.as_structure(obj)
    # Must use pymatgen structure else server does not know how to handle the JSON doc.
    structure.__class__ = pymatgen.Structure

    from abipy.core import restapi
    structures = []
    with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
        try:
            mpids = rest.find_structure(structure)
            if mpids:
                structures = [Structure.from_mpid(mid, final=final, api_key=api_key, endpoint=endpoint)
                        for mid in mpids]

        except rest.Error as exc:
            cprint(str(exc), "red")

        finally:
            # Back to abipy structure
            structure = Structure.as_structure(structure)
            return restapi.MpStructures(structures=structures, ids=mpids)


def mp_search(chemsys_formula_id, api_key=None, endpoint=None):
    """
    Connect to the materials project database.
    Get a list of structures corresponding to a chemical system, formula, or materials_id.

    Args:
        chemsys_formula_id (str): A chemical system (e.g., Li-Fe-O),
            or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).
        api_key (str): A String API key for accessing the MaterialsProject REST interface.
            If this is None, the code will check if there is a `PMG_MAPI_KEY` in your .pmgrc.yaml.
        endpoint (str): Url of endpoint to access the MaterialsProject REST interface.

    Returns:
        :class:`MpStructures` object with
            List of Structure objects, Materials project ids associated to structures.
            and List of dictionaries with MP data (same order as structures).

        Note that the attributes evalute to False if no match is found
    """
    chemsys_formula_id = chemsys_formula_id.replace(" ", "")

    structures, mpids, data = [], [], None
    from abipy.core import restapi
    with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
        try:
            data = rest.get_data(chemsys_formula_id, prop="")
            if data:
                structures = [Structure.from_str(d["cif"], fmt="cif", primitive=False, sort=False)
                                for d in data]
                mpids = [d["material_id"] for d in data]
                # Want AbiPy structure.
                structures = list(map(Structure.as_structure, structures))

        except rest.Error as exc:
            cprint(str(exc), "magenta")

        return restapi.MpStructures(structures, mpids, data=data)


def cod_search(formula, primitive=False):
    """
    Connect to the COD_ database. Get a list of structures corresponding to a chemical formula

    Args:
        formula (str): Chemical formula (e.g., Fe2O3)
        primitive (bool): True if primitive structures are wanted. Note that many COD structures are not primitive.

    Returns:
        :class:`CodStructures` object with
            List of Structure objects, COD ids associated to structures.
            and List of dictionaries with COD data (same order as structures).

        Note that the attributes evalute to False if no match is found
    """
    from pymatgen.ext.cod import COD
    data = COD().get_structure_by_formula(formula)

    cod_ids = [e.pop("cod_id") for e in data]
    # Want AbiPy structure.
    structures = list(map(Structure.as_structure, [e.pop("structure") for e in data]))
    if primitive:
        structures = [s.get_primitive_structure() for s in structures]

    from abipy.core import restapi
    return restapi.CodStructures(structures, cod_ids, data=data)


class Structure(pymatgen.Structure, NotebookWriter):
    """
    Extends :class:`pymatgen.core.structure.Structure` with Abinit-specific methods.

    A jupyter_ notebook documenting the usage of this object is available at
    <https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/structure.ipynb>

    For the pymatgen project see :cite:`Ong2013`.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Structure
    """
    @classmethod
    def as_structure(cls, obj):
        """
        Convert obj into a |Structure|. Accepts:

            - Structure object.
            - Filename
            - Dictionaries (JSON_ format or dictionaries with abinit variables).
            - Objects with a ``structure`` attribute.
        """
        if isinstance(obj, cls): return obj
        if isinstance(obj, pymatgen.Structure):
            obj.__class__ = cls
            return obj

        if is_string(obj):
            return cls.from_file(obj)

        if isinstance(obj, collections.Mapping):
            try:
                return Structure.from_abivars(obj)
            except:
                try:
                    return Structure.from_dict(obj)
                except:
                    raise TypeError("Don't know how to convert dict %s into a structure" % str(obj))

        if hasattr(obj, "structure"):
            return cls.as_structure(obj.structure)
        elif hasattr(obj, "final_structure"):
            # This for HIST.nc file
            return cls.as_structure(obj.final_structure)

        raise TypeError("Don't know how to convert %s into a structure" % type(obj))

    @classmethod
    def from_file(cls, filepath, primitive=False, sort=False):
        """
        Reads a structure from a file. For example, anything ending in
        a "cif" is assumed to be a Crystallographic Information Format file.
        Supported formats include CIF_, POSCAR/CONTCAR, CHGCAR, LOCPOT,
        vasprun.xml, CSSR, Netcdf and pymatgen's JSON serialized structures.

        Netcdf files supported:
            All files produced by ABINIT with info of the crystalline geometry
            HIST.nc, in this case the last structure of the history is returned.

        Args:
            filename (str): The filename to read from.
            primitive (bool): Whether to convert to a primitive cell
                Only available for cifs, POSCAR, CSSR, JSON, YAML
                Defaults to True.
            sort (bool): Whether to sort sites. Default to False.

        Returns: |Structure| object
        """
        if filepath.endswith("_HIST.nc"):
            # Abinit history file. In this case we return the last structure!
            # Note that HIST does not follow the etsf-io conventions.
            from abipy.dynamics.hist import HistFile
            with HistFile(filepath) as hist:
                return hist.structures[-1]

        elif filepath.endswith(".nc"):
            # Generic netcdf file.
            ncfile, closeit = as_etsfreader(filepath)

            new = ncfile.read_structure(cls=cls)
            new.set_abi_spacegroup(AbinitSpaceGroup.from_file(ncfile))
            if closeit: ncfile.close()

        elif filepath.endswith(".abi") or filepath.endswith(".in"):
            # Abinit input file.
            # Here I assume that the input file contains a single structure.
            from abipy.abio.abivars import AbinitInputFile
            return AbinitInputFile.from_file(filepath).structure

        elif filepath.endswith(".abo") or filepath.endswith(".out"):
            # Abinit output file. We can have multi-datasets and multiple initial/final structures!
            # By desing, we return the last structure if out is completed else the initial one.
            # None is returned if the structures are different.
            from abipy.abio.outputs import AbinitOutputFile
            with AbinitOutputFile(filepath) as out:
                #print("initial_structures:\n", out.initial_structures, "\nfinal_structures:\n", out.final_structures)
                if out.final_structures: return out.final_structure
                if out.initial_structures: return out.initial_structure
            raise ValueError("Cannot find structure in Abinit output file `%s`" % filepath)

        elif filepath.endswith("_DDB"):
            # DDB file.
            from abipy.abilab import abiopen
            with abiopen(filepath) as abifile:
                return abifile.structure

        elif filepath.endswith(".pickle"):
            # From pickle.
            with open(filepath, "rb") as fh:
                new = pickle.load(fh)
                if not isinstance(new, pymatgen.Structure):
                    # Is it a object with a structure property?
                    if hasattr(new, "structure"): new = new.structure

                if not isinstance(new, pymatgen.Structure):
                    raise TypeError("Don't know how to extract a Structure from file %s, received type %s" %
                        (filepath, type(new)))

                if new.__class__ != cls: new.__class__ = cls

        else:
            # Invoke pymatgen and change class
            # Note that AbinitSpacegroup is missing here.
            new = super(Structure, cls).from_file(filepath, primitive=primitive, sort=sort)
            if new.__class__ != cls: new.__class__ = cls

        return new

    @classmethod
    def from_mpid(cls, material_id, final=True, api_key=None, endpoint=None):
        """
        Get a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string, e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.
            api_key (str): A String API key for accessing the MaterialsProject
                REST interface. Please apply on the Materials Project website for one.
                If this is None, the code will check if there is a ``PMG_MAPI_KEY`` in your .pmgrc.yaml.
                If so, it will use that environment
                This makes easier for heavy users to simply add this environment variable
                to their setups and MPRester can then be called without any arguments.
            endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
                Defaults to the standard Materials Project REST address, but
                can be changed to other urls implementing a similar interface.

        Returns: |Structure| object.
        """
        # Get pytmatgen structure and convert it to abipy structure
        from abipy.core import restapi
        with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
            new = rest.get_structure_by_material_id(material_id, final=final)
            return cls.as_structure(new)

    @classmethod
    def from_cod_id(cls, cod_id, primitive=False, **kwargs):
        """
        Queries the COD_ for a structure by id. Returns |Structure| object.

        Args:
            cod_id (int): COD id.
            primitive (bool): True if primitive structures are wanted. Note that many COD structures are not primitive.
            kwargs: Arguments passed to ``get_structure_by_id``

        Returns: |Structure| object.
        """
        from pymatgen.ext.cod import COD
        new = COD().get_structure_by_id(cod_id, **kwargs)
        if primitive: new = new.get_primitive_structure()
        return cls.as_structure(new)

    @classmethod
    def from_ase_atoms(cls, atoms):
        """
        Returns structure from ASE Atoms.

        Args:
            atoms: ASE Atoms object

        Returns:
            Equivalent Structure
        """
        import pymatgen.io.ase as aio
        return aio.AseAtomsAdaptor.get_structure(atoms, cls=cls)

    def to_ase_atoms(self):
        """
        Returns ASE_ Atoms object from structure.
        """
        import pymatgen.io.ase as aio
        return aio.AseAtomsAdaptor.get_atoms(self)

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
        molecule = pymatgen.Molecule([p.symbol for p in pseudos], cart_coords)
        l = pmg_units.ArrayWithUnit(acell, "bohr").to("ang")

        new = molecule.get_boxed_structure(l[0], l[1], l[2])
        return cls.as_structure(new)

    @classmethod
    def boxed_atom(cls, pseudo, cart_coords=3*(0,), acell=3*(10,)):
        """
        Creates an atom in a periodic box of lengths acell [Bohr]

        Args:
            pseudo: Pseudopotential object.
            cart_coords: Cartesian coordinates in Angstrom
            acell: Lengths of the box in *Bohr* (Abinit input variable)
        """
        return cls.boxed_molecule([pseudo], cart_coords, acell=acell)

    @classmethod
    def bcc(cls, a, species, primitive=True, units="ang", **kwargs):
        """
        Build a primitive or a conventional bcc crystal structure.

        Args:
            a: Lattice parameter (Angstrom if units is not given)
            species: Chemical species. See __init__ method of |pymatgen-Structure|
            primitive: if True a primitive cell will be produced, otherwise a conventional one
            units: Units of input lattice parameters e.g. "bohr", "pm"
            kwargs: All keyword arguments accepted by |pymatgen-Structure|.
        """
        a = pmg_units.Length(a, units).to("ang")
        if primitive:
            lattice = 0.5 * float(a) * np.array([
                -1,  1,  1,
                 1, -1,  1,
                 1,  1, -1])

            coords = [[0, 0, 0]]

        else:
            lattice = float(a) * np.eye(3)
            coords = [[0, 0, 0],
                      [0.5, 0.5, 0.5]]
            species = np.repeat(species, 2)

        return cls(lattice, species, coords=coords,  **kwargs)

    @classmethod
    def fcc(cls, a, species, primitive=True, units="ang", **kwargs):
        """
        Build a primitive or a conventional fcc crystal structure.

        Args:
            a: Lattice parameter (Angstrom if units is not given)
            species: Chemical species. See __init__ method of :class:`pymatgen.Structure`
            primitive: if True a primitive cell will be produced, otherwise a conventional one
            units: Units of input lattice parameters e.g. "bohr", "pm"
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structure`
        """
        a = pmg_units.Length(a, units).to("ang")
        if primitive:
            lattice = 0.5 * float(a) * np.array([
                0,  1,  1,
                1,  0,  1,
                1,  1,  0])
            coords = [[0, 0, 0]]
        else:
            lattice = float(a) * np.eye(3)
            species = np.repeat(species, 4)
            coords = [[0, 0, 0],
                      [0.5, 0.5, 0],
                      [0.5, 0, 0.5],
                      [0, 0.5, 0.5]]

        return cls(lattice, species, coords=coords, **kwargs)

    @classmethod
    def zincblende(cls, a, species, units="ang", **kwargs):
        """
        Build a primitive zincblende crystal structure.

        Args:
            a: Lattice parameter (Angstrom if units is not given)
            species: Chemical species. See __init__ method of :class:`pymatgen.Structure`
            units: Units of input lattice parameters e.g. "bohr", "pm"
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structure`

        Example::

            Structure.zincblende(a, ["Zn", "S"])

        """
        a = pmg_units.Length(a, units).to("ang")
        lattice = 0.5 * float(a) * np.array([
            0,  1,  1,
            1,  0,  1,
            1,  1,  0])

        frac_coords = np.reshape([0, 0, 0, 0.25, 0.25, 0.25], (2, 3))
        return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    @classmethod
    def rocksalt(cls, a, species, units="ang", **kwargs):
        """
        Build a primitive fcc crystal structure.

        Args:
            a: Lattice parameter (Angstrom if units is not given)
            units: Units of input lattice parameters e.g. "bohr", "pm"
            species: Chemical species. See __init__ method of :class:`pymatgen.Structure`
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structure`

        Example::

            Structure.rocksalt(a, ["Na", "Cl"])

        """
        a = pmg_units.Length(a, units).to("ang")
        lattice = 0.5 * float(a) * np.array([
            0,  1,  1,
            1,  0,  1,
            1,  1,  0])

        frac_coords = np.reshape([0, 0, 0, 0.5, 0.5, 0.5], (2, 3))
        return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    @classmethod
    def ABO3(cls, a, species, units="ang", **kwargs):
       """
       Peroviskite structures.

       Args:
            a: Lattice parameter (Angstrom if units is not given)
            species: Chemical species. See __init__ method of :class:`pymatgen.Structure`
            units: Units of input lattice parameters e.g. "bohr", "pm"
            kwargs: All keyword arguments accepted by :class:`pymatgen.Structure`
       """
       a = pmg_units.Length(a, units).to("ang")
       lattice = float(a) * np.eye(3)
       frac_coords = np.reshape([
          0,     0,   0,  # A (2a)
          0.5, 0.5, 0.5,  # B (2a)
          0.5, 0.5, 0.0,  # O (6b)
          0.5, 0.0, 0.5,  # O (6b)
          0.0, 0.5, 0.5,  # O (6b)
         ], (5, 3))

       return cls(lattice, species, frac_coords, coords_are_cartesian=False, **kwargs)

    @classmethod
    def from_abistring(cls, string):
        """Initialize Structure from string with Abinit input variables."""
        from abipy.abio.abivars import AbinitInputFile
        return AbinitInputFile.from_string(string).structure

    @classmethod
    def from_abivars(cls, *args, **kwargs):
        """
        Build a |Structure| object from a dictionary with ABINIT variables.

        Example::

            al_structure = Structure.from_abivars(
                acell=3*[7.5],
                rprim=[0.0, 0.5, 0.5,
                       0.5, 0.0, 0.5,
                       0.5, 0.5, 0.0],
                typat=1,
                xred=[0.0, 0.0, 0.0],
                ntypat=1,
                znucl=13,
            )

        ``xred`` can be replaced with ``xcart`` or ``xangst``.
        """
        return structure_from_abivars(cls, *args, **kwargs)

    def __str__(self):
        return self.to_string()

    def to_string(self, title=None, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))
        if verbose:
            app(self.spget_summary(verbose=verbose))
        else:
            app(super(Structure, self).__str__())

        if self.abi_spacegroup is not None:
            app("\nAbinit Spacegroup: %s" % self.abi_spacegroup.to_string(verbose=verbose))

        return "\n".join(lines)

    def to(self, fmt=None, filename=None, **kwargs):
        __doc__ = pymatgen.Structure.to.__doc__ + \
        "\n Accepts also fmt='abivars' and `.abi` as Abinit input file extension"

        filename = filename or ""
        fmt = "" if fmt is None else fmt.lower()
        fname = os.path.basename(filename)

        if fmt in ("abi", "abivars") or fname.endswith(".abi"):
            if filename:
                with open(filename, "wt") as f:
                    f.write(self.abi_string)
            else:
                return self.abi_string
        else:
            return super(Structure, self).to(fmt=fmt, filename=filename, **kwargs)

    def __mul__(self, scaling_matrix):
        """
        Makes a supercell. Allowing to have sites outside the unit cell
        See pymatgen for docs.

        Wraps __mul__ operator of pymatgen structure to return abipy structure
        """
        new = super(Structure, self).__mul__(scaling_matrix)
        return self.__class__.as_structure(new)

    __rmul__ = __mul__

    def to_abivars(self, **kwargs):
        """Returns a dictionary with the ABINIT variables."""
        return structure_to_abivars(self, **kwargs)

    @property
    def latex_formula(self):
        """LaTeX formatted formula. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$."""
        from pymatgen.util.string import latexify
        return latexify(self.formula)

    @property
    def abi_string(self):
        """Return a string with the ABINIT input associated to this structure."""
        from abipy.abio.variable import InputVariable
        lines = []
        app = lines.append
        abivars = self.to_abivars()
        for varname, value in abivars.items():
            app(str(InputVariable(varname, value)))

        return("\n".join(lines))

    def get_conventional_standard_structure(self, international_monoclinic=True,
                                            symprec=1e-3, angle_tolerance=5):
        """
        Gives a structure with a conventional cell according to certain
	standards. The standards are defined in :cite:`Setyawan2010`
        They basically enforce as much as possible norm(a1) < norm(a2) < norm(a3)

        Returns:
            The structure in a conventional standardized cell
        """
        spga = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        new = spga.get_conventional_standard_structure(international_monoclinic=international_monoclinic)
        return self.__class__.as_structure(new)

    def abi_primitive(self, symprec=1e-3, angle_tolerance=5, no_idealize=0):
        #TODO: this should be moved to pymatgen in the get_refined_structure or so ...
        # to be considered in February 2016
        import spglib
        from pymatgen.io.ase import AseAtomsAdaptor
        try:
            from ase.atoms import Atoms
        except ImportError:
            raise ImportError('Could not import Atoms from ase. Install it with `conda install ase` or pip')

        s = self.get_sorted_structure()
        ase_adaptor = AseAtomsAdaptor()
        ase_atoms = ase_adaptor.get_atoms(structure=s)
        standardized = spglib.standardize_cell(ase_atoms, to_primitive=1, no_idealize=no_idealize,
                                               symprec=symprec, angle_tolerance=angle_tolerance)
        standardized_ase_atoms = Atoms(scaled_positions=standardized[1], numbers=standardized[2], cell=standardized[0])
        standardized_structure = ase_adaptor.get_structure(standardized_ase_atoms)

        return self.__class__.as_structure(standardized_structure)

    def abi_sanitize(self, symprec=1e-3, angle_tolerance=5, primitive=True, primitive_standard=False):
        """
        Returns a new structure in which:

            * Structure is refined.
            * Reduced to primitive settings.
            * Lattice vectors are exchanged if the triple product is negative

        Args:
            symprec (float): Symmetry precision used to refine the structure.
            angle_tolerance (float): Tolerance on angles.
                if ``symprec`` is None and `angle_tolerance` is None, no structure refinement is peformed.
            primitive (bool): Whether to convert to a primitive cell following :cite:`Setyawan2010`
            primitive_standard (bool): Returns most primitive structure found.
        """
        from pymatgen.transformations.standard_transformations import PrimitiveCellTransformation, SupercellTransformation
        structure = self.__class__.from_sites(self)

        # Refine structure
        if symprec is not None and angle_tolerance is not None:
            sym_finder = SpacegroupAnalyzer(structure=structure, symprec=symprec, angle_tolerance=angle_tolerance)
            structure = sym_finder.get_refined_structure()

        # Convert to primitive structure.
        if primitive:
            if primitive_standard:
                # Setyawan, W., & Curtarolo, S.
                sym_finder_prim = SpacegroupAnalyzer(structure=structure, symprec=symprec, angle_tolerance=angle_tolerance)
                structure = sym_finder_prim.get_primitive_standard_structure(international_monoclinic=False)
            else:
                # Find most primitive structure.
                get_prim = PrimitiveCellTransformation()
                structure = get_prim.apply_transformation(structure)

        # Exchange last two lattice vectors if triple product is negative.
        m = structure.lattice.matrix
        x_prod = np.dot(np.cross(m[0], m[1]), m[2])
        if x_prod < 0:
            trans = SupercellTransformation(((1, 0, 0), (0, 0, 1), (0, 1, 0)))
            structure = trans.apply_transformation(structure)
            m = structure.lattice.matrix
            x_prod = np.dot(np.cross(m[0], m[1]), m[2])
            if x_prod < 0: raise RuntimeError("x_prod is still negative!")

        return self.__class__.as_structure(structure)

    def get_oxi_state_decorated(self, **kwargs):
        """
        Use :class:`pymatgen.analysis.bond_valence.BVAnalyzer` to estimate oxidation states
        Return oxidation state decorated structure.
        This currently works only for ordered structures only.

        Args:
            kwargs: Arguments passed to BVAnalyzer

        Returns:
            A modified structure that is oxidation state decorated.
        """
        from pymatgen.analysis.bond_valence import BVAnalyzer
        new = BVAnalyzer(**kwargs).get_oxi_state_decorated_structure(self)
        return self.__class__.as_structure(new)

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        try:
            from nbjsmol import nbjsmol_display
            return nbjsmol_display(self.to(fmt="cif"), ext=".cif", html=True)
        except ImportError as exc:
            # Write warning only once.
            cls = self.__class__
            if not hasattr(cls, "_repr_html_num_warns"): cls._repr_html_num_warns = 0
            if cls._repr_html_num_warns == 0:
                cls._repr_html_num_warns += 1
                warn(str(exc) +
                     "\n_repr_html_ requires nbjsmol package\n."
                     "Install it with `pip install nbjsmol.`\n"
                     "See also https://github.com/gmatteo/nbjsmol\n"
                     "Returning `str(self)` in HTML form.")
            return str(self).replace("\n", "<br>")

    @property
    def reciprocal_lattice(self):
        """
        Reciprocal lattice of the structure.
        """
        return self._lattice.reciprocal_lattice

    def lattice_vectors(self, space="r"):
        """
        Returns the vectors of the unit cell in Angstrom.

        Args:
            space: "r" for real space vectors, "g" for reciprocal space basis vectors.
        """
        if space.lower() == "r":
            return self.lattice.matrix
        if space.lower() == "g":
            return self.lattice.reciprocal_lattice.matrix
        raise ValueError("Wrong value for space: %s " % str(space))

    def spget_lattice_type(self, symprec=1e-3, angle_tolerance=5):
        """
        Call spglib to get the lattice for the structure, e.g., (triclinic,
        orthorhombic, cubic, etc.).This is the same than the
        crystal system with the exception of the hexagonal/rhombohedral lattice

        Args:
            symprec (float): Symmetry precision for distance
            angle_tolerance (float): Tolerance on angles.

        Returns:
            (str): Lattice type for structure or None if type cannot be detected.
        """
        spgan = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        return spgan.get_lattice_type()

    def spget_equivalent_atoms(self, symprec=1e-3, angle_tolerance=5, printout=False):
        """
        Call spglib_ to find the inequivalent atoms and build symmetry tables.

        Args:
            symprec (float): Symmetry precision for distance.
            angle_tolerance (float): Tolerance on angles.
            printout (bool): True to print symmetry tables.

        Returns:
            ``namedtuple`` (irred_pos, eqmap, spgdata) with the following attributes::

                * irred_pos: array giving the position of the i-th irred atom in the structure.
                    The number of irred atoms is len(irred_pos)
                *   eqmap: Mapping irred atom position --> list with positions of symmetrical atoms
                *   spgdata: spglib dataset with additional data reported by spglib_.

         :Example:

            for irr_pos in irred_pos:
                eqmap[irr_pos]   # List of symmetrical positions associated to the irr_pos atom.
        """
        spgan = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        spgdata = spgan.get_symmetry_dataset()
        equivalent_atoms = spgdata["equivalent_atoms"]
        irred_pos = []
        eqmap = collections.defaultdict(list)
        for pos, eqpos in enumerate(equivalent_atoms):
            eqmap[eqpos].append(pos)
            # Add it to irred_pos if it's irreducible.
            if pos == eqpos: irred_pos.append(pos)

        # Convert to numpy arrays
        irred_pos = np.array(irred_pos)
        for eqpos in eqmap:
            eqmap[eqpos] = np.array(eqmap[eqpos], dtype=np.int)

        if printout:
            print("Found %d inequivalent position(s)." % len(irred_pos))
            for i, irr_pos in enumerate(sorted(eqmap.keys())):
                print("Irred_Site: %s" % str(self[irr_pos]))
                for eqind in eqmap[irr_pos]:
                    if eqind == irr_pos: continue
                    print("\tSymEq: %s" % str(self[eqind]))

        return dict2namedtuple(irred_pos=irred_pos, eqmap=eqmap, spgdata=spgdata)

    def spget_summary(self, symprec=1e-3, angle_tolerance=5, verbose=0):
        """
        Return string with full information about crystalline structure i.e.
        space group, point group, wyckoff positions, equivalent sites.

        Args:
            symprec (float): Symmetry precision for distance.
            angle_tolerance (float): Tolerance on angles.
            verbose (int): Verbosity level.
        """
        spgan = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        spgdata = spgan.get_symmetry_dataset()
        # Get spacegroup number computed by Abinit if available.
        abispg_number = None if self.abi_spacegroup is None else self.abi_spacegroup.spgid

        # Print lattice info
        outs = ["Full Formula ({s})".format(s=self.composition.formula),
                "Reduced Formula: {}".format(self.composition.reduced_formula)]
        app = outs.append
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        app("")
        app("Spglib space group info (magnetic symmetries are not taken into account).")
        app("Spacegroup: %s (%s), Hall: %s, Abinit spg_number: %s" % (
             spgan.get_space_group_symbol(), spgan.get_space_group_number(), spgan.get_hall(), str(abispg_number)))
        app("Crystal_system: %s, Lattice_type: %s, Point_group: %s" % (
            spgan.get_crystal_system(), spgan.get_lattice_type(), spgan.get_point_group_symbol()))
        app("")

        wickoffs, equivalent_atoms = spgdata["wyckoffs"], spgdata["equivalent_atoms"]
        table = [["Idx", "Symbol", "Reduced_Coords", "Wyck", "EqIdx"]]
        for i, site in enumerate(self):
            table.append([
                i,
                site.species_string,
                "%+.5f %+.5f %+.5f" % tuple(site.frac_coords),
                "%s" % wickoffs[i],
                "%d" % equivalent_atoms[i],
            ])

        from tabulate import tabulate
        app(tabulate(table, headers="firstrow"))

        # Print entire dataset.
        if verbose > 1:
            app("\nSpglib dataset:")
            app(pformat(spgdata, indent=4))

        return "\n".join(outs)

    @property
    def abi_spacegroup(self):
        """
        :class:`AbinitSpaceGroup` instance with Abinit symmetries read from the netcd file.
        None if abinit symmetries are not available e.g. if the structure has been created
        from a CIF file.
        """
        try:
            return self._abi_spacegroup
        except AttributeError:
            return None

    def set_abi_spacegroup(self, spacegroup):
        """``AbinitSpaceGroup`` setter."""
        self._abi_spacegroup = spacegroup

    @property
    def has_abi_spacegroup(self):
        """True is the structure contains info on the spacegroup."""
        return self.abi_spacegroup is not None

    def spgset_abi_spacegroup(self, has_timerev, overwrite=False):
        """
        Call spglib to find the spacegroup of the crystal, create new
	:class:`AbinitSpaceGroup` object and store it in ``self.abi_spacegroup``.

        Args:
            has_timerev (bool): True if time-reversal can be used.
            overwrite (bool): By default, the method raises `ValueError` if the object
                already has the list of symmetries found by Abinit.

        Returns:
	    :class:`AbinitSpaceGroup`

        .. warning:

            This method should be called only if the Abipy structure does not have
            spacegroup symmetries e.g. if we are reading a CIF file or if the structure
            is initialized from an output file produced by another code.
        """
        if self.has_abi_spacegroup and not overwrite:
            raise ValueError(("Structure object already has an Abinit spacegroup object.\n"
                              "Use `overwrite=True` to allow modification."))

        msg = ("Structure object does not have symmetry operations computed from Abinit.\n"
               "Will call spglib to get symmetry operations.")
        cprint(msg, "magenta")

        spglib_data = SpacegroupAnalyzer(self).get_symmetry_dataset()
        spgid = spglib_data["number"]
        symrel, tnons = spglib_data["rotations"], spglib_data["translations"]
        symafm = [1] * len(symrel)  # TODO: Anti-ferromagnetic symmetries are not supported by spglib

        abispg = AbinitSpaceGroup(spgid, symrel, tnons, symafm, has_timerev, inord="C")
        self.set_abi_spacegroup(abispg)

        return abispg

    def abiget_spginfo(self, tolsym=None, pre=None):
        """
	Call Abinit to get spacegroup information.
	Return dictionary with e.g. {'bravais': 'Bravais cF (face-center cubic)', 'spg_number': 227, 'spg_symbol': 'Fd-3m'}.

	Args:
            tolsym: Abinit tolsym input variable. None correspondes to the default value.
	    pre: Keywords in dictionary are prepended with this string
        """
        from abipy.data.hgh_pseudos import HGH_TABLE
        from abipy.abio import factories
        gsinp = factories.gs_input(self, HGH_TABLE, spin_mode="unpolarized")
        gsinp["chkprim"] = 0
        d = gsinp.abiget_spacegroup(tolsym=tolsym, retdict=True)
        if pre: d = {pre + k: v for k, v in d.items()}
        return d

    def print_neighbors(self, radius=2.0):
        """
        Get neighbors for each atom in the unit cell, out to a distance ``radius`` in Angstrom
        Print results.
        """
        print(" ")
        print("Finding neighbors for each atom in the unit cell, out to a distance %s [Angstrom]" % radius)
        print(" ")

        ns = self.get_all_neighbors(radius, include_index=False)
        for i, (site, sited_list) in enumerate(zip(self, ns)):
            print("[%s] site %s has %s neighbors:" % (i, repr(site), len(sited_list)))
            for s, dist in sorted(sited_list, key=lambda t: t[1]):
                print("\t", repr(s), " at distance", dist)
            print("")

    @lazy_property
    def hsym_kpath(self):
        """
        Returns an instance of :class:`pymatgen.symmetry.bandstructure.HighSymmKpath`.
        (Database of high symmetry k-points and high symmetry lines).
        """
        from pymatgen.symmetry.bandstructure import HighSymmKpath
        return HighSymmKpath(self)

    @lazy_property
    def hsym_kpoints(self):
        """|KpointList| object with the high-symmetry K-points."""
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

    def get_kcoords_from_names(self, knames, cart_coords=False):
        """
        Return numpy array with the fractional coordinates of the high-symmetry k-points listed in `knames`.

        Args:
            knames: List of strings with the k-point labels.
            cart_coords: True if the ``coords`` dataframe should contain Cartesian cordinates
                instead of Reduced coordinates.
        """
        kname2frac = {k.name: k.frac_coords for k in self.hsym_kpoints}
        # Add aliases for Gamma.
        if r"$\Gamma$" in kname2frac:
            kname2frac["G"] = kname2frac[r"$\Gamma$"]
            kname2frac["Gamma"] = kname2frac[r"$\Gamma$"]

        try:
            kcoords = np.reshape([kname2frac[name] for name in list_strings(knames)], (-1, 3))
        except KeyError:
            cprint("Internal list of high-symmetry k-points:\n" % str(self.hsym_kpoints))
            raise

        if cart_coords:
            kcoords = self.reciprocal_lattice.get_cartesian_coords(kcoords)

        return kcoords

    @lazy_property
    def hsym_stars(self):
        """
        List of |KpointStar| objects. Each star is associated to one of the special k-points
        present in the pymatgen database.
        """
        # Construct the stars.
        return [kpoint.compute_star(self.abi_spacegroup.fm_symmops) for kpoint in self.hsym_kpoints]

    # TODO
    #def get_star_kpoint(self, kpoint):

    #    # Call spglib to get spacegroup if Abinit spacegroup is not available.
    #    if self.abi_spacegroup is None:
    #        self.spgset_abi_spacegroup(has_timerev=not options.no_time_reversal)

    #    kpoint = Kpoint(options.kpoint, self.reciprocal_lattice)
    #    kstar = kpoint.compute_star(self.abi_spacegroup, wrap_tows=True)
    #    return kstar
    #    #print("Found %s points in the star of %s\n" % (len(kstar), repr(kpoint)))
    #    #for k in kstar:
    #    #    print(4 * " ", repr(k))

    def get_sorted_structure_z(self):
        """Order the structure according to increasing Z of the elements"""
        return self.__class__.from_sites(sorted(self.sites, key=lambda site: site.specie.Z))

    def findname_in_hsym_stars(self, kpoint):
        """
        Returns the name of the special k-point, None if kpoint is unknown.
        """
        if self.abi_spacegroup is None: return None

        from .kpoints import Kpoint
        kpoint = Kpoint.as_kpoint(kpoint, self.reciprocal_lattice)

        # Try to find kpoint in hsym_stars without taking into accout symmetry operation (compare with base_point)
        # Important if there are symmetry equivalent k-points in hsym_kpoints e.g. K and U in FCC lattice
        # as U should not be mapped onto K as done in the second loop below.
        from .kpoints import issamek
        for star in self.hsym_stars:
            if issamek(kpoint.frac_coords, star.base_point.frac_coords):
                return star.name

	# Now check if kpoint is in one of the stars.
        for star in self.hsym_stars:
            i = star.find(kpoint)
            if i != -1:
                #print("input kpt:", kpoint, "star image", star[i], star[i].name)
                return star.name
        else:
            return None

    def get_symbol2indices(self):
        """
        Return a dictionary mapping chemical symbols to numpy array with the position of the atoms.

        Example:

            MgB2 --> {Mg: [0], B: [1, 2]}
        """
        return {symbol: np.array(self.indices_from_symbol(symbol)) for symbol in self.symbol_set}

    def get_symbol2coords(self):
        """
        Return a dictionary mapping chemical symbols to a [ntype_symbol, 3] numpy array
        with the fractional coordinates.
        """
        # TODO:
        #use structure.frac_coords but add reshape in pymatgen.
        #fcoords = np.reshape([s.frac_coords for s in self], (-1, 3))
        coords = {}
        for symbol in self.symbol_set:
            coords[symbol] = np.reshape(
                [site.frac_coords for site in self if site.specie.symbol == symbol], (-1, 3))

        return coords

    def dot(self, coords_a, coords_b, space="r", frac_coords=False):
        """
        Compute the scalar product of vector(s) either in real space or
        reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        lattice = {"r": self.lattice,
                   "g": self.reciprocal_lattice}[space.lower()]
        return lattice.dot(coords_a, coords_b, frac_coords=frac_coords)

    def norm(self, coords, space="r", frac_coords=True):
        """
        Compute the norm of vector(s) either in real space or reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        return np.sqrt(self.dot(coords, coords, space=space, frac_coords=frac_coords))

    def get_dict4pandas(self, symprec=1e-2, angle_tolerance=5.0, with_spglib=True):
        """
        Return a :class:`OrderedDict` with the most important structural parameters:

            - Chemical formula and number of atoms.
            - Lattice lengths, angles and volume.
            - The spacegroup number computed by Abinit (set to None if not available).
            - The spacegroup number and symbol computed by spglib (if `with_spglib`).

        Useful to construct pandas DataFrames

        Args:
            with_spglib (bool): If True, spglib is invoked to get the spacegroup symbol and number
            symprec (float): Symmetry precision used to refine the structure.
            angle_tolerance (float): Tolerance on angles.
        """
        abc, angles = self.lattice.abc, self.lattice.angles

        # Get spacegroup info from spglib.
        spglib_symbol, spglib_number, spglib_lattice_type = None, None, None
        if with_spglib:
            try:
                spglib_symbol, spglib_number = self.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
                spglib_lattice_type = self.spget_lattice_type(symprec=symprec, angle_tolerance=angle_tolerance)
            except Exception as exc:
                cprint("Spglib couldn't find space group symbol and number for composition %s" % str(self.composition), "red")
                print("Exception:\n", exc)

        # Get spacegroup number computed by Abinit if available.
        abispg_number = None if self.abi_spacegroup is None else self.abi_spacegroup.spgid

        od = OrderedDict([
            ("formula", self.formula), ("natom", self.num_sites),
            ("alpha", angles[0]), ("beta", angles[1]), ("gamma", angles[2]),
            ("a", abc[0]), ("b", abc[1]), ("c", abc[2]), ("volume", self.volume),
            ("abispg_num", abispg_number),
        ])
        if with_spglib:
            od["spglib_symb"] = spglib_symbol
            od["spglib_num"] = spglib_number
            od["spglib_lattice_type"] = spglib_lattice_type

        return od

    @add_fig_kwargs
    def plot(self, **kwargs):
        """
        Plot structure with matplotlib. Return matplotlib Figure
        See plot_structure for kwargs
        """
        from abipy.tools.plotting import plot_structure
        return plot_structure(self, **kwargs)

    @add_fig_kwargs
    def plot_bz(self, ax=None, pmg_path=True, with_labels=True, **kwargs):
        """
        Gives the plot (as a matplotlib object) of the symmetry line path in the Brillouin Zone.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            pmg_path (bool): True if the default path used in pymatgen should be show.
            with_labels (bool): True to plot k-point labels.

        Returns: |matplotlib-Figure|.
        """
        from pymatgen.electronic_structure.plotter import plot_brillouin_zone, plot_brillouin_zone_from_kpath
        labels = None if not with_labels else self.hsym_kpath.kpath["kpoints"]
        #pprint(labels)
        if pmg_path:
            return plot_brillouin_zone_from_kpath(self.hsym_kpath, ax=ax, show=False, **kwargs)
        else:
            return plot_brillouin_zone(self.reciprocal_lattice, ax=ax, labels=labels, show=False, **kwargs)

    @add_fig_kwargs
    def plot_xrd(self, wavelength="CuKa", symprec=0, debye_waller_factors=None,
                 two_theta_range=(0, 90), annotate_peaks=True, ax=None, **kwargs):
        """
        Use pymatgen :class:`XRDCalculator` to show the XRD plot.

        Args:
            wavelength (str/float): The wavelength can be specified as either a
                float or a string. If it is a string, it must be one of the
                supported definitions in the AVAILABLE_RADIATION class
                variable, which provides useful commonly used wavelengths.
                If it is a float, it is interpreted as a wavelength in
                angstroms. Defaults to "CuKa", i.e, Cu K_alpha radiation.
            symprec (float): Symmetry precision for structure refinement. If
                set to 0, no refinement is done. Otherwise, refinement is
                performed using spglib_ with provided precision.
            debye_waller_factors ({element symbol: float}): Allows the
                specification of Debye-Waller factors. Note that these
                factors are temperature dependent.
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (bool): Whether to annotate the peaks with plane information.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        xrd = XRDCalculator(wavelength=wavelength, symprec=symprec, debye_waller_factors=debye_waller_factors)
        xrd.get_plot(self, two_theta_range=two_theta_range, annotate_peaks=annotate_peaks, ax=ax)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.plot_bz(show=False)

    def export(self, filename, visu=None, verbose=1):
        """
        Export the crystalline structure to file ``filename``.

        Args:
            filename (str): String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example,
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            visu: |Visualizer| subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See |Visualizer| for the list of applications and formats supported.
            verbose: Verbosity level

        Returns: ``Visulizer`` instance.
        """
        if "." not in filename:
            raise ValueError("Cannot detect extension in filename %s:" % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]
        #print("tokens", tokens, "ext", ext)
        #if ext == "POSCAR":

        if not tokens[0]:
            # filename == ".ext" ==> Create temporary file.
            # nbworkdir in cwd is needed when we invoke the method from a notebook.
            from abipy.core.globals import abinb_mkstemp
            _, rpath = abinb_mkstemp(force_abinb_workdir=False, use_relpath=False,
                                     suffix="." + ext, text=True)
            #if abilab.in_notebook():
            #    _, filename = tempfile.mkstemp(suffix="." + ext, dir=abilab.get_abipy_nbworkdir(), text=True)
            #else:
            #    _, filename = tempfile.mkstemp(suffix="." + ext, text=True)

        if ext.lower() in ("xsf", "poscar", "cif"):
            if verbose:
                print("Writing data to:", filename, "with fmt:", ext.lower())
            s = self.to(fmt=ext)
            with open(filename, "wt") as fh:
                fh.write(s)

        if visu is None:
            return Visualizer.from_file(filename)
        else:
            return visu(filename)

    def chemview(self, **kwargs):
        """
        Visualize structure in jupyter_ notebook using chemview package.
        """
        from pymatgen.vis.structure_chemview import quick_view
        return quick_view(self, **kwargs)

    def vtkview(self, show=True, **kwargs):
        """
        Visualize structure with VTK. Requires vtk python bindings.

        Args:
            show: True to show structure immediately.
            kwargs: keyword arguments passed to :class:`StructureVis`.

        Return:
            StructureVis object.
        """
        from pymatgen.vis.structure_vtk import StructureVis
        vis = StructureVis(**kwargs)
        vis.set_structure(self, to_unit_cell=True)
        if show: vis.show()
        return vis

    def mayaview(self, figure=None, show=True, **kwargs):
        """Visualize the crystalline structure with mayaview"""
        from abipy.display import mvtk
        return mvtk.plot_structure(self, figure=figure, show=show, **kwargs)

    def visualize(self, appname="vesta"):
        """
        Visualize the crystalline structure with visualizer.
        See |Visualizer| for the list of applications and formats supported.
        """
        if appname in ("mpl", "matplotlib"): return self.plot()
        if appname == "vtk": return self.vtkview()
        if appname == "mayavi": return self.mayaview()

        # Get the Visualizer subclass from the string.
        visu = Visualizer.from_name(appname)

        # Try to export data to one of the formats supported by the visualizer
        # Use a temporary file (note "." + ext)
        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export(ext, visu=visu)()
            except visu.Error as exc:
                print(exc)
                pass
        else:
            raise visu.Error("Don't know how to export data for %s" % appname)

    def convert(self, fmt="cif", **kwargs):
        """
        Return string with the structure in the given format `fmt`
        Options include "abivars", "cif", "xsf", "poscar", "siesta", "wannier90", "cssr", "json".
        """
        if fmt in ("abivars", "abinit"):
            return self.abi_string
        elif fmt == "abipython":
            return pformat(self.to_abivars(), indent=4)
        elif fmt == "qe":
            from pymatgen.io.pwscf import PWInput
            return str(PWInput(self, pseudo={s: s + ".pseudo" for s in self.symbol_set}))
        elif fmt == "siesta":
            return structure2siesta(self)
        elif fmt in ("wannier90", "w90"):
            from abipy.wannier90.win import structure2wannier90
            return structure2wannier90(self)
        else:
            return super(Structure, self).to(fmt=fmt, **kwargs)

    #def max_overlap_and_sites(self, pseudos):
    #    # For each site in self:
    #    # 1) Get the radius of the pseudopotential sphere
    #    # 2) Get the neighbors of the site (considering the periodic images).

    #    pseudos = PseudoTable.as_table(pseudos)
    #    max_overlap, ovlp_sites = 0.0, None
    #    for site in self:
    #        symbol = site.specie.symbol
    #        pseudo = pseudos[symbol]
    #        r1 = Length(pseudo.r_cut, "Bohr").to("ang")
    #        sitedist_list = self.get_neighbors(site, r1, include_index=False)

    #        if sitedist_list:
    #            # Spheres are overlapping: compute overlap and update the return values
    #            # if the new overlap is larger than the previous one.
    #            for other_site, dist in sitedist_list:
    #                other_symbol = other_site.specie.symbol
    #                other_pseudo = pseudos[other_symbol]
    #                r2 = Length(other_pseudo.r_cut, "Bohr").to("ang")
    #                # Eq 16 of http://mathworld.wolfram.com/Sphere-SphereIntersection.html
    #                overlap = sphere_overlap(site.coords, r1, other_site.coords, r2)
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
        displ = np.reshape(displ, (-1, 3)).copy()

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
        scale_matrix = np.zeros((3, 3), dtype=np.int)
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
                        if dnorm < (dmin-1e-6) and dnorm > 1e-6:
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
                            if dnorm < (dmin-1e-6) and dnorm > 1e-6:
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
                            if dnorm < (dmin-1e-6) and dnorm > 1e-6:
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

        # find the translation vectors (in terms of the initial lattice vectors)
        # that are inside the unit cell defined by the scale matrix
        # we're using a slightly offset interval from 0 to 1 to avoid numerical
        # precision issues
        inv_matrix = np.linalg.inv(scale_matrix)

        frac_points = np.dot(all_points, inv_matrix)
        tvects = all_points[np.where(np.all(frac_points < 1-1e-10, axis=1)
                                     & np.all(frac_points >= -1e-10, axis=1))]
        assert len(tvects) == np.round(abs(np.linalg.det(scale_matrix)))

        return tvects

    def write_vib_file(self, xyz_file, qpoint, displ, do_real=True, frac_coords=True,
                       scale_matrix=None, max_supercell=None):
        """
        write into the file descriptor xyz_file the positions and displacements of the atoms

        Args:
            xyz_file: file_descriptor
            qpoint: qpoint to be analyzed
            displ: eigendisplacements to be analyzed
            do_real: True if you want to get only real part, False means imaginary part
            frac_coords: True if the eigendisplacements are given in fractional coordinates
            scale_matrix: Scale matrix for supercell
            max_supercell: Maximum size of supercell vectors with respect to primitive cell
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

        for at, site in enumerate(self):
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

                xyz_file.write(fmtstr.format(site.specie, coords[0], coords[1], coords[2],
                               new_displ[0], new_displ[1], new_displ[2]))

    def frozen_2phonon(self, qpoint, displ1, displ2, eta=1, frac_coords=False, scale_matrix=None, max_supercell=None):
        """
        Creates the supercell needed for a given qpoint and adds the displacements.
        The displacements are normalized so that the largest atomic displacement will correspond to the
        value of eta in Angstrom.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space.
            displ1: first displacement in real space of the atoms.
            displ2: second displacement in real space of the atoms.
            eta: pre-factor multiplying the displacement. Gives the value in Angstrom of the
                largest displacement.
            frac_coords: whether the displacements are given in fractional or cartesian coordinates
            scale_matrix: the scaling matrix of the supercell. If None a scaling matrix suitable for
                the qpoint will be determined.
            max_supercell: mandatory if scale_matrix is None, ignored otherwise. Defines the largest
                supercell in the search for a scaling matrix suitable for the q point.

        Returns:
            A namedtuple with a Structure with the displaced atoms, a numpy array containing the
            displacements applied to each atom and the scale matrix used to generate the supercell.
        """

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

        if frac_coords:
            displ1 = np.array((old_lattice.get_cartesian_coords(d) for d in displ1))
            displ2 = np.array((old_lattice.get_cartesian_coords(d) for d in displ2))
        else:
            displ1 = np.array(displ1)
            displ2 = np.array(displ2)
        # from here displ are in cartesian coordinates

        norm_factor = np.linalg.norm(displ1+displ2, axis=1).max()

        displ1 = eta * displ1 / norm_factor
        displ2 = eta * displ2 / norm_factor

        new_displ1 = np.zeros(3, dtype=np.float)
        new_displ2 = np.zeros(3, dtype=np.float)
        new_sites = []
        displ_list = []
        for at,site in enumerate(self):
            for t in tvects:
                new_displ1[:] = np.real(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ1[at,:])

                new_displ2[:] = np.real(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ2[at,:])

                displ_list.append(new_displ1+new_displ2)
                coords = site.coords + old_lattice.get_cartesian_coords(t) + new_displ1 + new_displ2
                new_site = PeriodicSite(
                    site.species_and_occu, coords, new_lattice,
                    coords_are_cartesian=True, properties=site.properties,
                    to_unit_cell=True)
                new_sites.append(new_site)

        new_structure = self.__class__.from_sites(new_sites)

        return dict2namedtuple(structure=new_structure, displ=np.array(displ_list), scale_matrix=scale_matrix)

    def frozen_phonon(self, qpoint, displ, eta=1, frac_coords=False, scale_matrix=None, max_supercell=None):
        """
        Creates a supercell with displaced atoms for the specified q-point.
        The displacements are normalized so that the largest atomic displacement will correspond to the
        value of eta in Angstrom.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space.
            displ: displacement in real space of the atoms.
            eta: pre-factor multiplying the displacement. Gives the value in Angstrom of the
                largest displacement.
            frac_coords: whether the displacements are given in fractional or cartesian coordinates
            scale_matrix: the scaling matrix of the supercell. If None a scaling matrix suitable for
                the qpoint will be determined.
            max_supercell: mandatory if scale_matrix is None, ignored otherwise. Defines the largest
                supercell in the search for a scaling matrix suitable for the q point.

        Returns:
            A namedtuple with a Structure with the displaced atoms, a numpy array containing the
            displacements applied to each atom and the scale matrix used to generate the supercell.
        """

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

        if frac_coords:
            displ = np.array((old_lattice.get_cartesian_coords(d) for d in displ))
        else:
            displ = np.array(displ)
        # from here displ are in cartesian coordinates

        displ = eta * displ / np.linalg.norm(displ, axis=1).max()

        new_displ = np.zeros(3, dtype=np.float)
        new_sites = []
        displ_list = []
        for at, site in enumerate(self):
            for t in tvects:
                new_displ[:] = np.real(np.exp(2*1j*np.pi*(np.dot(qpoint,t)))*displ[at,:])
                displ_list.append(list(new_displ))

                coords = site.coords + old_lattice.get_cartesian_coords(t) + new_displ
                new_site = PeriodicSite(
                    site.species_and_occu, coords, new_lattice,
                    coords_are_cartesian=True, properties=site.properties,
                    to_unit_cell=True)
                new_sites.append(new_site)

        new_structure = self.__class__.from_sites(new_sites)

        return dict2namedtuple(structure=new_structure, displ=np.array(displ_list), scale_matrix=scale_matrix)

    def calc_kptbounds(self):
        """Returns the suggested value for the ABINIT variable ``kptbounds``."""
        kptbounds = [k.frac_coords for k in self.hsym_kpoints]
        return np.reshape(kptbounds, (-1, 3))

    def get_kpath_input_string(self, fmt="abinit", line_density=10):
        """
        Return string with input variables for band-structure calculations
        in the format used by code `fmt`.
        Use `line_density` points for the smallest segment (if supported by code).
        """
        lines = []; app = lines.append
        if fmt in ("abinit", "abivars"):
            app("# Abinit Structure")
            app(self.convert(fmt=fmt))
            app("\n# K-path in reduced coordinates:")
            app("# tolwfr 1e-20 iscf -2 getden ??")
            app(" ndivsm %d" % line_density)
            app(" kptopt %d" % -(len(self.hsym_kpoints) - 1))
            app(" kptbounds")
            for k in self.hsym_kpoints:
                app("    {:+.5f}  {:+.5f}  {:+.5f}  # {kname}".format(*k.frac_coords, kname=k.name))

        elif fmt in ("wannier90", "w90"):
            app("# Wannier90 structure")
            from abipy.wannier90.win import Wannier90Input
            win = Wannier90Input(self)
            win.set_kpath()
            app(win.to_string())

        elif fmt == "siesta":
            app("# Siesta structure")
            app(self.convert(fmt=fmt))
            # Build normalized k-path.
            from .kpoints import Kpath
            vertices_names = [(k.frac_coords, k.name) for k in self.hsym_kpoints]
            kpath = Kpath.from_vertices_and_names(self, vertices_names, line_density=line_density)
            app("%block BandLines")
            prev_ik = 0
            for ik, k in enumerate(kpath):
                if not k.name: continue
                n = ik - prev_ik
                app("{}  {:+.5f}  {:+.5f}  {:+.5f}  # {kname}".format(n if n else 1, *k.frac_coords, kname=k.name))
                prev_ik = ik
            app("%endblock BandLines")

        else:
            raise ValueError("Don't know how to generate string for code: `%s`" % str(fmt))

        return "\n".join(lines)

    #def ksampling_from_jhudb(self, **kwargs):
    #    from pymatgen.ext.jhu import get_kpoints
    #    __doc__ = get_kpoints.__doc__
    #    kpoints = get_kpoints(self, **kwargs)
    #    print(kpoints)
    #    print(kpoints.style)
    #    print("num_kpts", kpoints.num_kpts)

    #    d = {"kptopt": 0,
    #         "kpt": kpoints.kpts,
    #         "nkpt": kpoints.num_kpts,
    #         #"kptnrm": kptnrm,
    #         "wtk": kpoints.kpts_weights,
    #         "shiftk": kpoints.kpts_shift,
    #         "chksymbreak": 0,
    #    }
    #    print(d)

    #    #from pymatgen.io.abinit.abiobjects import KSampling
    #    #return KSampling(mode=KSamplingModes.automatic,
    #    #         num_kpts= 0,
    #    #         kpts=((1, 1, 1),),
    #    #         kpt_shifts=(0.5, 0.5, 0.5),
    #    #         kpts_weights=None, use_symmetries=True, use_time_reversal=True, chksymbreak=None,
    #    #         comment=None)

    def calc_ksampling(self, nksmall, symprec=0.01, angle_tolerance=5):
        """
        Return the k-point sampling from the number of divisions ``nksmall`` to be used for
        the smallest reciprocal lattice vector.
        """
        ngkpt = self.calc_ngkpt(nksmall)
        shiftk = self.calc_shiftk(symprec=symprec, angle_tolerance=angle_tolerance)

        return AttrDict(ngkpt=ngkpt, shiftk=shiftk)

    def calc_ngkpt(self, nksmall):
        """
        Compute the ABINIT variable ``ngkpt`` from the number of divisions used
        for the smallest lattice vector.

        Args:
            nksmall (int): Number of division for the smallest lattice vector.
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
        Find the values of ``shiftk`` and ``nshiftk`` appropriated for the sampling of the Brillouin zone.

        When the primitive vectors of the lattice do NOT form a FCC or a BCC lattice,
        the usual (shifted) Monkhorst-Pack grids are formed by using nshiftk=1 and shiftk 0.5 0.5 0.5 .
        This is often the preferred k point sampling. For a non-shifted Monkhorst-Pack grid,
        use `nshiftk=1` and `shiftk 0.0 0.0 0.0`, but there is little reason to do that.

        When the primitive vectors of the lattice form a FCC lattice, with rprim::

                0.0 0.5 0.5
                0.5 0.0 0.5
                0.5 0.5 0.0

        the (very efficient) usual Monkhorst-Pack sampling will be generated by using nshiftk= 4 and shiftk::

            0.5 0.5 0.5
            0.5 0.0 0.0
            0.0 0.5 0.0
            0.0 0.0 0.5

        When the primitive vectors of the lattice form a BCC lattice, with rprim::

               -0.5  0.5  0.5
                0.5 -0.5  0.5
                0.5  0.5 -0.5

        the usual Monkhorst-Pack sampling will be generated by using nshiftk= 2 and shiftk::

                0.25  0.25  0.25
               -0.25 -0.25 -0.25

        However, the simple sampling nshiftk=1 and shiftk 0.5 0.5 0.5 is excellent.

        For hexagonal lattices with hexagonal axes, e.g. rprim::

                1.0  0.0       0.0
               -0.5  sqrt(3)/2 0.0
                0.0  0.0       1.0

        one can use nshiftk= 1 and shiftk 0.0 0.0 0.5
        In rhombohedral axes, e.g. using angdeg 3*60., this corresponds to shiftk 0.5 0.5 0.5,
        to keep the shift along the symmetry axis.

        Returns:
            Suggested value of shiftk.
        """
        # Find lattice type.
        sym = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        lattice_type, spg_symbol = sym.get_lattice_type(), sym.get_space_group_symbol()

        # Check if the cell is primitive
        is_primitve = len(sym.find_primitive()) == len(self)

        # Generate the appropriate set of shifts.
        shiftk = None

        if is_primitve:
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

            elif lattice_type == "tetragonal":
                if "I" in spg_symbol:
                    # BCT
                    shiftk = [0.25,  0.25,  0.25,
                             -0.25, -0.25, -0.25]

        if shiftk is None:
            # Use default value.
            shiftk = [0.5, 0.5, 0.5]

        return np.reshape(shiftk, (-1, 3))

    def num_valence_electrons(self, pseudos):
        """
        Returns the number of valence electrons.

        Args:
            pseudos: List of |Pseudo| objects or list of filenames.
        """
        nval, table = 0, PseudoTable.as_table(pseudos)
        for site in self:
            pseudo = table.pseudo_with_symbol(site.specie.symbol)
            nval += pseudo.Z_val

        return int(nval) if int(nval) == nval else nval

    def valence_electrons_per_atom(self, pseudos):
        """
        Returns the number of valence electrons for each atom in the structure.

        Args:
            pseudos: List of |Pseudo| objects or list of filenames.
        """
        table = PseudoTable.as_table(pseudos)
        psp_valences = []
        for site in self:
            pseudo = table.pseudo_with_symbol(site.specie.symbol)
            psp_valences.append(pseudo.Z_val)

        return psp_valences

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        # The notebook will reconstruct the object from this file
        _, tmpfile = tempfile.mkstemp(suffix='.pickle')
        with open(tmpfile, "wb") as fh:
            pickle.dump(self, fh)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("structure = abilab.Structure.from_file('%s')" % tmpfile),
            nbv.new_code_cell("print(structure)"),
            nbv.new_code_cell("print(structure.abi_string)"),
            nbv.new_code_cell("structure"),
            nbv.new_code_cell("print(structure.spglib_summary())"),
            nbv.new_code_cell("if structure.abi_spacegroup is not None: print(structure.abi_spacegroup)"),
            nbv.new_code_cell("print(structure.hsym_kpoints)"),
            nbv.new_code_cell("structure.plot_bz();"),
            nbv.new_code_cell("structure.plot_xrd();"),
            nbv.new_code_cell("# sanitized = structure.abi_sanitize(); print(sanitized)"),
            nbv.new_code_cell("# ase_atoms = structure.to_ase_atoms()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


def dataframes_from_structures(struct_objects, index=None, symprec=1e-2, angle_tolerance=5,
	                       with_spglib=True, cart_coords=False):
    """
    Build two pandas Dataframes_ with the most important geometrical parameters associated to
    a list of structures or a list of objects that can be converted into structures.

    Args:
        struct_objects: List of objects that can be converted to structure.
            Support filenames, structure objects, Abinit input files, dicts and many more types.
            See ``Structure.as_structure`` for the complete list.
        index: Index of the |pandas-DataFrame|.
        symprec (float): Symmetry precision used to refine the structure.
        angle_tolerance (float): Tolerance on angles.
        with_spglib (bool): If True, spglib_ is invoked to get the spacegroup symbol and number.
        cart_coords: True if the ``coords`` dataframe should contain Cartesian cordinates
            instead of Reduced coordinates.

    Return:
        namedtuple with two |pandas-DataFrames| named ``lattice`` and ``coords``
        ``lattice`` contains the lattice parameters. ``coords`` the atomic positions..
        The list of structures is available in the ``structures`` entry.

    .. code-block:: python

        dfs = dataframes_from_structures(files)
        dfs.lattice
        dfs.coords
        for structure in dfs.structures:
            print(structure)
    """
    structures = [Structure.as_structure(obj) for obj in struct_objects]
    # Build Frame with lattice parameters.
    # Use OrderedDict to have columns ordered nicely.
    odict_list = [(structure.get_dict4pandas(with_spglib=with_spglib, symprec=symprec, angle_tolerance=angle_tolerance))
	          for structure in structures]

    import pandas as pd
    lattice_frame = pd.DataFrame(odict_list, index=index,
                                 columns=list(odict_list[0].keys()) if odict_list else None)

    # Build Frame with atomic positions.
    vtos = lambda v: "%+0.6f %+0.6f %+0.6f" % (v[0], v[1], v[2])
    max_numsite = max(len(s) for s in structures)
    odict_list = []
    for structure in structures:
        if cart_coords:
            odict_list.append({i: (site.species_string, vtos(site.coords)) for i, site in enumerate(structure)})
        else:
            odict_list.append({i: (site.species_string, vtos(site.frac_coords)) for i, site in enumerate(structure)})

    coords_frame = pd.DataFrame(odict_list, index=index,
                                columns=list(range(max_numsite)) if odict_list else None)

    return dict2namedtuple(lattice=lattice_frame, coords=coords_frame, structures=structures)


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
            scaling_matrix: A scaling matrix for transforming the lattice vectors.
                Has to be all integers. Several options are possible:

                a. A full 3x3 scaling matrix defining the linear combination of the old lattice vectors.
                    E.g., [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with lattice vectors
                    a' = 2a + b, b' = 3b, c' = c
                    where a, b, and c are the lattice vectors of the original structure.
                b. A sequence of three scaling factors. e.g., [2, 1, 1]
                   specifies that the supercell should have dimensions 2a x b x c.
                c. A number, which simply scales all lattice vectors by the same factor.

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
            displ: Displacement vector with 3*len(self) entries (fractional coordinates).
            eta: Scaling factor.
            frac_coords: Boolean stating whether the vector corresponds to fractional or cartesian coordinates.

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

    def frozen_phonon(self, qpoint, displ, eta=1, frac_coords=False, scale_matrix=None, max_supercell=None):

        return self._original_structure.frozen_phonon(qpoint, displ, eta, frac_coords, scale_matrix, max_supercell)

    def frozen_2phonon(self, qpoint, displ1, displ2, eta=1, frac_coords=False, scale_matrix=None, max_supercell=None):

        return self._original_structure.frozen_2phonon(qpoint, displ1, displ2, eta, frac_coords, scale_matrix,
                                                       max_supercell)


def diff_structures(structures, fmt="cif", mode="table", headers=(), file=sys.stdout):
    """
    Convert list of structure to string using format `fmt`, print diff to file `file`.

    Args:
        structures: List of structures or list of objects that can be converted into structure e.g. filepaths
        fmt: Any output format supported by `structure.to` method. Non-case sensitive.
        mode: `table` to show results in tabular form or `diff` to show differences with unified diff.
        headers: can be an explicit list of column headers Otherwise a headerless table is produced
        file: Output Stream
    """
    outs = [s.convert(fmt=fmt).splitlines() for s in map(Structure.as_structure, structures)]

    if mode == "table":
        try:
            from itertools import izip_longest as zip_longest  # Py2
        except ImportError:
            from itertools import zip_longest  # Py3k

        table = [r for r in zip_longest(*outs, fillvalue=" ")]
        from tabulate import tabulate
        print(tabulate(table, headers=headers), file=file)

    elif mode == "diff":
        import difflib
        fromfile, tofile = "", ""
        for i in range(1, len(outs)):
            if headers: fromfile, tofile = headers[0], headers[i]
            diff = "\n".join(difflib.unified_diff(outs[0], outs[i], fromfile=fromfile, tofile=tofile))
            print(diff, file=file)

    else:
        raise ValueError("Unsupported mode: `%s`" % str(mode))


def structure2siesta(structure, verbose=0):
    """
    Return string with structural information in Siesta format from pymatgen structure

    Args:
        structure: pymatgen structure.
        verbose: Verbosity level.
    """

    if not structure.is_ordered:
        raise NotImplementedError("""\
Received disordered structure with partial occupancies that cannot be converted into a Siesta input
Please use OrderDisorderedStructureTransformation or EnumerateStructureTransformation
to build an appropriate supercell from partial occupancies or alternatively use the Virtual Crystal Approximation.""")

    types_of_specie = structure.types_of_specie

    lines = []
    app = lines.append
    app("NumberOfAtoms %d" % len(structure))
    app("NumberOfSpecies %d" % structure.ntypesp)

    if verbose:
        app("# The species number followed by the atomic number, and then by the desired label")
    app("%block ChemicalSpeciesLabel")
    for itype, specie in enumerate(types_of_specie):
        app("    %d %d %s" % (itype + 1, specie.number, specie.symbol))
    app("%endblock ChemicalSpeciesLabel")

    # Write lattice vectors.
    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    lvectors = np.where(np.abs(structure.lattice.matrix) > 1e-8, structure.lattice.matrix, 0.0)
    app("LatticeConstant 1.0 Ang")
    app("%block LatticeVectors")
    for r in lvectors:
        app("    %.10f %.10f %.10f" % (r[0], r[1], r[2]))
    app("%endblock LatticeVectors")

    # Write atomic coordinates
    #% block AtomicCoordinatesAndAtomicSpecies
    #4.5000 5.0000 5.0000 1
    #5.5000 5.0000 5.0000 1
    #% endblock AtomicCoordinatesAndAtomicSpecies
    app("AtomicCoordinatesFormat Fractional")
    app("%block AtomicCoordinatesAndAtomicSpecies")
    for i, site in enumerate(structure):
        itype = types_of_specie.index(site.specie)
        fc = np.where(np.abs(site.frac_coords) > 1e-8, site.frac_coords, 0.0)
        app("    %.10f %.10f %.10f %d" % (fc[0], fc[1], fc[2], itype + 1))
    app("%endblock AtomicCoordinatesAndAtomicSpecies")

    return "\n".join(lines)