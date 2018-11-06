# coding: utf-8
"""The interatomic force constants calculated by anaddb."""
from __future__ import print_function, division, absolute_import # unicode_literals,

import numpy as np

from monty.functools import lazy_property
from abipy.core.mixins import Has_Structure
from abipy.iotools import ETSF_Reader
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt


class InteratomicForceConstants(Has_Structure):
    """
    The interatomic force constants calculated by anaddb.
    Read from anaddb.nc
    """

    def __init__(self, structure, atoms_indices, neighbours_indices, ifc_cart_coord,
                 ifc_cart_coord_short_range, local_vectors, distances):
        """
        Args:
            structure: |Structure| object.
            atoms_index: List of integers representing the indices in the structure of the analyzed atoms.
            neighbours_index: List of integers representing the indices in the structure of the neighbour atoms.
            ifc_cart_coord: ifc in Cartesian coordinates
            ifc_cart_coord_short_range: short range part of the ifc in Cartesian coordinates
            local_vectors: local basis used to determine the ifc_local_coord
        """
        self._structure = structure
        self.atoms_indices = atoms_indices
        self.neighbours_indices = neighbours_indices
        self.ifc_cart_coord = ifc_cart_coord
        self.ifc_cart_coord_short_range = ifc_cart_coord_short_range
        self.local_vectors = local_vectors
        self.distances = distances

    @property
    def number_of_atoms(self):
        """Number of atoms is structure."""
        return len(self.structure)

    @classmethod
    def from_file(cls, filepath):
        """Create the object from a netcdf_ file."""
        with ETSF_Reader(filepath) as r:
            try:
                structure = r.read_structure()
                atoms_indices = r.read_value("ifc_atoms_indices") - 1
                neighbours_indices = r.read_value("ifc_neighbours_indices") - 1
                distances = r.read_value("ifc_distances")
                ifc_cart_coord = r.read_value("ifc_matrix_cart_coord")
                ifc_cart_coord_short_range = r.read_value("ifc_matrix_cart_coord_short_range", default=None)
                local_vectors = r.read_value("ifc_local_vectors", default=None)
            except:
                import traceback
                msg = traceback.format_exc()
                msg += ("Error while trying to read IFCs from file.\n"
                       "Verify that the required variables are used in anaddb: ifcflag, natifc, atifc, ifcout\n")
                raise ValueError(msg)

            return cls(structure=structure, atoms_indices=atoms_indices, neighbours_indices=neighbours_indices,
                       ifc_cart_coord=ifc_cart_coord,
                       ifc_cart_coord_short_range=ifc_cart_coord_short_range, local_vectors=local_vectors,
                       distances=distances)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")

        return "\n".join(lines)

    @property
    def number_of_neighbours(self):
        """Number of neighbouring atoms for which the ifc are present. ifcout in anaddb."""
        return np.shape(self.neighbours_indices)[1]

    @lazy_property
    def ifc_cart_coord_ewald(self):
        """Ewald part of the ifcs in cartesian coordinates"""
        if self.ifc_cart_coord_short_range is None:
            return None
        else:
            return self.ifc_cart_coord-self.ifc_cart_coord_short_range

    @lazy_property
    def ifc_local_coord(self):
        """Ifcs in local coordinates"""
        if self.local_vectors is None:
            return None
        else:
            return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord, self.local_vectors)

    @lazy_property
    def ifc_local_coord_short_range(self):
        """Short range part of the ifcs in cartesian coordinates"""
        if self.local_vectors is None:
            return None
        else:
            return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord_short_range, self.local_vectors)

    @lazy_property
    def ifc_local_coord_ewald(self):
        """Ewald part of the ifcs in local coordinates"""
        return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord_ewald, self.local_vectors)

    def _filter_ifc_indices(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Internal method that provides the indices of the neighouring atoms in self.neighbours_indices that satisfy
        the required conditions. All the arguments are optional. If None the filter will not be applied.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """

        if atom_indices is not None and atom_element is not None:
            raise ValueError("atom_index and atom_element cannot be specified simultaneously")

        if atom_indices is not None and not isinstance(atom_indices, (list, tuple)):
            atom_indices = [atom_indices]

        if atom_element:
            atom_indices = self.structure.indices_from_symbol(atom_element)

        if atom_indices is None:
            atom_indices = range(len(self.structure))

        # apply the filter: construct matrices of num_atoms*num_neighbours size, all conditions should be satisfied.
        ind = np.where(
            (np.tile(np.in1d(self.atoms_indices, atom_indices), [self.number_of_neighbours, 1])).T &
            (self.distances > min_dist if min_dist is not None else True) &
            (self.distances < max_dist if max_dist is not None else True) &
            (np.in1d(self.neighbours_indices, self.structure.indices_from_symbol(neighbour_element))
             .reshape(self.number_of_atoms, self.number_of_neighbours) if neighbour_element is not None else True)
        )

        return ind

    def get_ifc_cartesian(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Filters the IFCs in cartesian coordinates
        All the arguments are optional. If None the filter will not be applied.
        Returns two arrays containing the distances and the corresponding filtered ifcs.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """
        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        return self.distances[ind], self.ifc_cart_coord[ind]

    def get_ifc_local(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Filters the IFCs in local coordinates
        All the arguments are optional. If None the filter will not be applied.
        Returns two arrays containing the distances and the corresponding filtered ifcs.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        return self.distances[ind], self.ifc_local_coord[ind]

    def get_plot_ifc(self, ifc, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None,
                     max_dist=None, ax=None, **kwargs):
        """
        Plots the specified ifcs, filtered according to the optional arguments.
        An array with shape number_of_atoms*number_of_neighbours, so only one of the components of the ifc matrix can
        be plotted at a time.

        Args:
            ifc: an array with shape number_of_atoms * number_of_neighbours of the ifc that should be plotted
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        dist, filtered_ifc =  self.distances[ind], ifc[ind]

        if 'color' not in kwargs: kwargs['color'] = 'blue'
        if 'marker' not in kwargs: kwargs['marker'] = 'o'
        if 'linewidth' not in kwargs and 'lw' not in kwargs: kwargs['lw'] = 0

        ax.set_xlabel('Distance (Bohr)')
        ax.set_ylabel(r'IFC (Ha/Bohr$^2$)')
        ax.grid(True)
        ax.plot(dist, filtered_ifc, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_longitudinal_ifc(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None,
                              max_dist=None, ax=None, **kwargs):
        """
        Plots the total longitudinal ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        return self.get_plot_ifc(self.ifc_local_coord[:, :, 0, 0], atom_indices=atom_indices, atom_element=atom_element,
                                 neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist, ax=ax, **kwargs)

    @add_fig_kwargs
    def plot_longitudinal_ifc_short_range(self, atom_indices=None, atom_element=None, neighbour_element=None,
                                          min_dist=None, max_dist=None, ax=None, **kwargs):
        """
        Plots the short range longitudinal ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        if self.ifc_local_coord_short_range is None:
            raise ValueError("Ewald contribution is missing, Run anaddb with dipdip=1")

        return self.get_plot_ifc(self.ifc_local_coord_short_range[:, :, 0, 0], atom_indices=atom_indices,
                                 atom_element=atom_element, neighbour_element=neighbour_element, min_dist=min_dist,
                                 max_dist=max_dist, ax=ax, **kwargs)

    @add_fig_kwargs
    def plot_longitudinal_ifc_ewald(self, atom_indices=None, atom_element=None, neighbour_element=None,
                                    min_dist=None, max_dist=None, ax=None, **kwargs):
        """
        Plots the Ewald part of the ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        if self.ifc_local_coord_ewald is None:
            raise ValueError("Ewald contribution is missing, Run anaddb with dipdip=1")

        return self.get_plot_ifc(self.ifc_local_coord_ewald[:, :, 0, 0], atom_indices=atom_indices,
                                 atom_element=atom_element, neighbour_element=neighbour_element, min_dist=min_dist,
                                 max_dist=max_dist, ax=ax, **kwargs)