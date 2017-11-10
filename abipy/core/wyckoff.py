"""
This module provides functions to generate crystalline structures from Wyckoff positions
or to retrieve Wyckoff parameters from a given structure.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

#from monty.functools import lazy_property
from collections import Sequence, OrderedDict


def _validate_params(params, wpos):
    """
    Check keys in params. Raises ValueError.
    """
    check = params.copy()
    for k in params:
        if k not in wpos.names:
            raise ValueError("Cannot find key: %s in symbol_names: %s" % (k, wpos.names))
        check.pop(k)
    if check:
        raise ValueError("Found unknown symbol names in %s" % check)


class WyckoffPositions(Sequence):
    """
    A numpy array with sympy expression.
    """
    @classmethod
    def from_string(cls, s):
        """
        Initialize the object from string `s`.
        """
        # Build a matrix of strings with expressions.
        s = s.replace("\n", "")
        s = s.replace("(", "").replace(")", "")
        s = s.replace("[", "").replace("]", "")

        try:
            matrix_str = np.reshape(s.split(","), (-1, 3))
        except ValueError:
            raise ValueError("Wrong string in input, perhaps missing comma. Input string:\n %s" % str(s))

        return cls(matrix_str)

    def __init__(self, matrix_str, letter=None, site_symmetry=None):
        """

        Args:
            matrix_str:
            letter:
            site_symmetry:
        """
        self.letter, self.site_symmetry = None, None

        # Build numpy matrix of sympy expressions from the matrix of strings.
        from sympy import Symbol
        from sympy.parsing.sympy_parser import parse_expr
        exprs = []
        for vec in np.reshape(matrix_str, (-1, 3)):
            exprs.extend([parse_expr(e) for e in vec])

        self._expr_mat = np.reshape(exprs, (-1, 3))

        # Get set of sympy symbols.
        symbols = []
        for expr in self.ravel():
            for arg in expr.args:
                if isinstance(arg, Symbol): symbols.append(arg)

        # SymPy symbols and names.
        self.symbols = list(set(symbols))
        self.names = [s.name for s in self.symbols]

    def __len__(self):
        return len(self._expr_mat)

    def __getitem__(self, key):
        return self._expr_mat.__getitem__(key)

    def __str__(self):
        """String representation."""
        lines = []
        for vec in self:
            lines.append("(" + ",".join(str(e) for e in vec) + ")")
        return ", ".join(lines)

    def ravel(self):
        """Flat list, similar to np.ravel"""
        return self._expr_mat.ravel()

    @property
    def mult(self):
        """Multiplicity"""
        return len(self)

    @property
    def num_params(self):
        """The number of parameters."""
        return len(self.symbols)

    def params_from_frac_coords(self, frac_coords):
        """
        Compute the Wyckoff parameters from the fractional coordinates.
        """
        frac_coords = np.reshape(frac_coords, (-1, 3))
        assert len(self) == len(frac_coords)
        equations = []
        for i in range(len(self)):
            eqs = self[i, :] - frac_coords[i, :]
            equations.extend(eqs)

        # Solve the problem.
        from sympy.solvers import solve
        d = solve(equations, self.symbols)
        #print("solution: ", d)

        # Return results in a dict with symbol names.
        if not d:
            return None

        return {symbol.name: d[symbol] for symbol in self.symbols}

    def error_of_params(self, params, frac_coords):
        """
        Return a numpy array with the difference between the
        wyckoff positions computed from the sympy expressions and the input frac_coords
        """
        _validate_params(params, self)
        frac_coords = np.reshape(frac_coords, (-1, 3))
        assert len(self) == len(frac_coords)

        frac_errors = np.empty((self.mult, 3))
        for i in range(len(self)):
            for j in range(3):
                frac_errors[i,j] = self[i,j].subs(params) - frac_coords[i,j]

        # Modulo lattice vector (with tolerance?)
        return frac_errors % 1


class Wyckoff(object):
    """
    """
    #@classmethod
    #def from_site2wpos_ordict(cls, site2wpos):
    #    return cls(site2wpos)

    def __init__(self, site2wpos):
        # TODO: Ordered dict.
        d = {}
        for k, s in site2wpos.items():
            d[k] = WyckoffPositions.from_string(s)

        self.site2wpos = d

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("Wyckoff positions:")
        for site, wpos in self.site2wpos.items():
            app("%s site:\n\t%s" % (site, wpos))

        return "\n".join(lines)

    def generate_structure(self, lattice, site2params, to_unit_cell=False, struct_cls=None):
        """
        Generate structure object from `lattice` and dictionary `site2params`
        """
        if not isinstance(site2params, OrderedDict):
            raise ValueError("Please use a OrderedDict for site2params so that\n" +
                             "the algorithm used to build the structure is deterministic.")

        # Build species and coords arrays.
        species, coords = [], []
        for elsym, wpos in self.site2wpos.items():
            params = site2params[elsym]
            _validate_params(params, wpos)
            for vec_expr in wpos:
                species.append(elsym)
                # Evaluate sympy expression with params.
                coords.append([float(expr.subs(params)) for expr in vec_expr])

        #for sp, c in zip(species, coords):
        #    print(sp, c, type(c), type(c[0]))

        from abipy.core.structure import Structure
        cls = struct_cls if struct_cls is not None else Structure
        return cls(lattice, species, coords, validate_proximity=True,
                   to_unit_cell=to_unit_cell, coords_are_cartesian=False, site_properties=None)

    def find_params(self, structure):
        """
        Compute the value of the parameter given a Structure object.
        """
        frac_coords = structure.get_symbol2coords()

        # Solve the problem.
        params = {}
        for elsym, wpos in self.site2wpos.items():
            params[elsym] = wpos.params_from_frac_coords(frac_coords[elsym])

        # Return results in a dict with symbol names.
        return params

    def error_of_params(self, params, structure):
        frac_coords = structure.get_symbol2coords()

        for elsym, wpos in self.site2wpos.items():
            frac_errors = wpos.error_of_params(params[elsym], frac_coords[elsym])
            print(frac_errors)
