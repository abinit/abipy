#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import numpy as np

from monty.functools import lazy_property
from collections import Sequence, OrderedDict


class WyckoffPositions(Sequence):
    """
    A numpy array with sympy expression.
    """
    @classmethod
    def from_string(cls, s):
        """Initialize the object from a string."""
        s = s.replace("\n", "")
        s = s.replace("(", "").replace(")", "")
        s = s.replace("[", "").replace("]", "")

        try:
            matrix = np.reshape(s.split(","), (-1, 3))
        except ValueError:
            raise ValueError("Wrong string in input, perhaps missing comma")

        from sympy.parsing.sympy_parser import parse_expr 
        exprs, symbols = [], []
        for vec in matrix:
            exprs.extend([parse_expr(e) for e in vec])

        return cls(exprs)

    def __init__(self, exprs, letter=None, site_symmetry=None):
        """

            Args:
                exprs:
                letter:
                site_symmetry:
        """
        self._expr_mat = np.reshape(exprs, (-1, 3))
        self.letter, self.site_symmetry = None, None

        # Get set of sympy symbols.
        from sympy import Symbol 
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

    def solve_frac_coords(self, frac_coords):
        frac_coords = np.reshape(frac_coords, (-1, 3))
        assert len(self) == len(frac_coords)
        equations = []
        for i in range(len(self)):
            eqs = self[i,:] - frac_coords[i,:]
            equations.extend(eqs)

        # Solve the problem.
        from sympy.solvers import solve
        d = solve(equations, self.symbols)
        #print("solution: ", d)

        # Return results in a dict with symbol names.
        if not d:
            return None

        return {symbol.name: d[symbol] for symbol in self.symbols}


class Wyckoff(object):
    """
    """
    def __init__(self, site2wpos):
        # Todo: Ordered dict.
        d = {}
        for k, s in site2wpos.items():
            d[k] = WyckoffPositions.from_string(s)

        self.site2wpos = d

    def __str__(self):
        lines = []
        app = lines.append
        app("Wyckoff positions:")
        for site, wpos in self.site2wpos.items():
            app("%s site:\n\t%s" % (site, wpos))
        return "\n".join(lines)

    def gen_structure(self, lattice, site2params, to_unit_cell=False, struct_cls=None):
        """
        """
        if not isinstance(site2params, OrderedDict):
            raise ValueError("Please use a OrderedDict for site2params so that\n" + 
                             "the algorithm used to build the structure is deterministic.")

        def check_params(params, wpos):
            """Check keys in params"""
            check = params.copy()
            for k in params:
                if k not in wpos.names:
                    raise ValueError("Cannot find %s in symbol_names: %s" % (k, wpos.names))
                check.pop(k)
            if check:
                raise ValueError("Found unknown symbol names in %s" % check)

        # Build species and coords arrays.
        species, coords = [], []
        for elsym, wpos in self.site2wpos.items():
            params = site2params[elsym]
            check_params(params, wpos)
            for vec_expr in wpos:
                species.append(elsym)
                # Evaluate sympy expression with params.
                coords.append([expr.subs(params) for expr in vec_expr])

        #print(list(zip(species, coords)))
        from abipy.core.structure import Structure
        cls = struct_cls if struct_cls is not None else Structure
        return cls(lattice, species, coords, validate_proximity=True,
                   to_unit_cell=to_unit_cell, coords_are_cartesian=False, site_properties=None)

    def find_params(self, structure):
        # TODO: 
        #use structure.frac_coords but add reshape in pymatgen.
        #fcoords = np.reshape([s.frac_coords for s in structure.sites], (-1, 3))

        fcoords = {}
        for elsym in structure.symbol_set:
            print(elsym)
            #inds = structure.indices_form_symbol(elsym)
            #print([site.specie.symbol for site in structure])
            fcoords[elsym] = np.reshape(
                [site.frac_coords for site in structure if site.specie.symbol == elsym], (-1, 3))
        #print(fcoords)

        # Solve the problem.
        params = {}
        for elsym, wpos in self.site2wpos.items():
            params[elsym] = wpos.solve_frac_coords(fcoords[elsym])

        # Return results in a dict with symbol names.
        return params


def test_sio2():
    # http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?gnum=152
    #6	c	1	
    site2wpos = {}

    #3	b	.2.	
    #site2wpos["Si"] = "(x,0,5/6)  (0,x,1/6)   (-x,-x,1/2)"

    #3	a	.2.	
    site2wpos["Si"] = "[x,0,1/3], (0,x,2/3), (-x,-x,0)"

    site2wpos["O"] = """
        (x,y,z), (-y,x-y,z+1/3), (-x+y,-x,z+2/3), (y,x,-z),
        (x-y,-y,-z+2/3), (-x,-x+y,-z+1/3)"""

    wyck = Wyckoff(site2wpos)
    print(wyck)

    site2params = OrderedDict([
        ("Si", dict(x=0.4763)),
        ("O", dict(x=0.1588, y=0.7439, z=0.4612)),
        #("O", dict(x=0.1588, y=0.7439, z=0.4612, yyy=3)),
    ])

    from abipy.core import Lattice
    lattice = abilab.Lattice.hexagonal(a=4.971, c=5.473)

    for to_unit_cell in [False]:
    #for to_unit_cell in [False, True]:
        structure = wyck.gen_structure(lattice, site2params, to_unit_cell=to_unit_cell)
        print("to_unit_cell %s\n" % to_unit_cell, structure)

        #from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        #spga = SpacegroupAnalyzer(structure)
        #structure = spga.get_refined_structure()
        #print("Refined:", structure)

        wyck_params = wyck.find_params(structure)

        for site, params in site2params.items():
            print(wyck_params[site])
            assert wyck_params[site] == params


if __name__ == "__main__":
    test_sio2()
