"""This module defines the structure object that store information on the crystalline structure and its symmetries."""
from __future__ import division, print_function

from .constants import Bohr2Ang
from .symmetries import SpaceGroup
from abipy.iotools import as_etsfreader, Visualizer
from abipy.iotools import xsf

import pymatgen

__all__ = [
    "Structure",
]


class Structure(pymatgen.Structure):

    @classmethod
    def from_file(cls, file):
        """
        Return a new instance from a NetCDF file containing crystallographic data in the ETSF-IO format.

        Args:
            ncdata:
                filename or NetcdfReader instance.
        """
        file, closeit = as_etsfreader(file)

        new = file.read_structure()
        # Change the class of new.
        new.__class__ = cls

        new.set_spacegroup(SpaceGroup.from_file(file))

        if closeit:
            file.close()

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
    def is_symmorphic(self):
        """True if at least one fractional translation is non-zero."""
        return self.spacegroup.issymmorphic

    @property
    def fm_symops(self):
        """Tuple with ferromagnetic symmetries (time-reversal is included, if present)."""
        return self.spacegroup.symmops(afm_sgn=+1)

    @property
    def afm_symops(self):
        """Tuple with Anti-ferromagnetic symmetries (time-reversal is included, if present)."""
        return self.spacegroup.symmops(afm_sgn=-1)

    @property
    def high_symm_kpath(self):
        """
        Returns an instance of the pymatgen class `HighSymmKpath`
        (Database of special k-point and high symmetry lines).
        """
        try:
            return self._high_symm_kpath
        except AttributeError:
            from pymatgen.symmetry.bandstructure import HighSymmKpath
            self._high_symm_kpath = HighSymmKpath(self)
            return self._high_symm_kpath

    @property
    def special_kpoints(self):
        try:
            return self._special_kpoints

        except AttributeError:
            name2frac_coords = self.high_symm_kpath.kpath["kpoints"]
            import pprint
            pprint.pprint(name2frac_coords)
            self._special_kpoints = name2frac_coords
            for symop in self.spacegroup:
                print(symop)
                symop.rotate_k( wrap_tows=True)

            return self._special_kpoints

    def show_bz(self):
        """
        Gives the plot (as a matplotlib object) of the symmetry line path in the Brillouin Zone.
        """
        return self.high_symm_kpath.get_kpath_plot()

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

        if not tokens[0]: # filename == ".ext" ==> Create temporary file.
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
        #from pymatgen.vis.structure_vtk import StructureVis
        #vis = StructureVis()
        #vis.set_structure(self)
        #vis.show()
        extensions = Visualizer.exts_from_appname(visualizer)

        for ext in extensions:
            ext = "." + ext
            try:
                return self.export(ext)
            except Visualizer.Error:
                pass
        else:
            raise Visualizer.Error("Don't know how to export data for %s" % visualizer)

