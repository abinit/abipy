"""
This module provides functions to generate crystalline structures from Wyckoff positions
or to retrieve Wyckoff parameters from a given structure.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import sympy as sp

#from monty.functools import lazy_property
from collections import Sequence, OrderedDict
from monty.termcolor import cprint
from abipy.core.mixins import Has_Structure


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


class SiteSymmetries(Has_Structure):

    def __init__(self, structure):
        # Get Spacegroup from spglib if not available from Abinit.
        if not structure.has_abi_spacegroup:
            structure.spgset_abi_spacegroup(has_timerev=True, overwrite=False)

        self._structure = structure
        abispg = structure.abi_spacegroup
        nsym = len(abispg.symrel)

        #self.eq_atoms = structure.spget_equivalent_atoms()

        # Precompute sympy objects.
        self.sp_symrel, self.sp_symrec = [], []
        self.sp_tnons, self.sp_inv_symrel = [], []
        from abipy.core.symmetries import mati3inv

        for symr, symc, tau in zip(abispg.symrel, abispg.symrec, abispg.tnons):
            self.sp_symrel.append(sp.Matrix((symr)))
            inv_symr = mati3inv(symr, trans=False)
            self.sp_inv_symrel.append(sp.Matrix((inv_symr)))
            self.sp_symrec.append(sp.Matrix((symc)))
            # FIXME: Should convert to rational
            # Permissible translations are unit cell translations or fractions thereof
            # that are consistent with the rotational symmetry (e.g. 1/2, 1/3, 1/4, and 1/6), plus combinations.
            tau = np.around(tau, decimals=5)
            self.sp_tnons.append(sp.Matrix(tau))

        #from abipy.core.symmetries import indsym_from_symrel
        #other_indsym = indsym_from_symrel(abispg.symrel, abispg.tnons, self.structure, tolsym=1e-8)
        #assert np.all(self.structure.indsym == other_indsym)

        # Compute symmetry operations in Cartesian coordinates (numpy arrays)
        a = self.structure.lattice.matrix.T
        self.symcart = np.matmul(a, np.matmul(abispg.symrel, np.linalg.inv(a)))

        indsym = self.structure.indsym

        import spglib
        self.sitesym_labels = []
        for iatom, site in enumerate(self.structure):
            rotations = [abispg.symrel[isym] for isym in range(nsym) if
                         indsym[iatom, isym, 3] == iatom and abispg.symafm[isym] == 1]
            # Passing a 0-length rotations list to spglib can segfault.
            herm_symbol, ptg_num = "1", 1
            if len(rotations) != 0:
                herm_symbol, ptg_num, trans_mat = spglib.get_pointgroup(rotations)

            self.sitesym_labels.append("%s (#%d) ord: %d" % (herm_symbol.strip(), ptg_num, len(rotations)))

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def get_wyckoff_dataframe(self, view="all", select_symbols=None, verbose=0):
        """
        Find Wyckoff positions.

        Args:
            view:
            select_symbols:
            verbose: Verbosity level.

        Return |pandas-DataFrame| with cartesian tensor components as columns and (inequivalent) sites along the rows.
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols, verbose=verbose)

        sitesym_labels = self.structure.spget_site_symmetries()

        frac_symbols = sp.symbols("xfrac, yfrac, zfrac")
        vector = sp.Matrix((frac_symbols))

        indsym = self.structure.indsym
        abispg = self.structure.abi_spacegroup
        rows = []
        for (iatom, wlabel) in zip(aview.iatom_list, aview.wyck_labels):
            site = self.structure[iatom]
            system = []
            for isym, (rm1, tau) in enumerate(zip(self.sp_inv_symrel, self.sp_tnons)):
                if indsym[iatom, isym, 3] != iatom: continue
                l0 = indsym[iatom, isym, :3]
                l0 = sp.Matrix(l0)
                m = rm1 * (vector - tau) - (l0 + vector)
                if verbose: sp.pprint(m)
                system.append(m)

            # Solve system of linear equations.
            solutions = sp.solve(system, dict=True)
            #print(solutions)
            if verbose and not solutions:
                cprint("No solution for iatom %d" % iatom, "yellow")

            if solutions:
                d = OrderedDict()
                d["element"] = site.specie.symbol
                d["site_index"] = iatom
                d["cart_coords"] = np.round(site.coords, decimals=5)
                d["frac_coords"] = np.round(site.frac_coords, decimals=5)
                d["wyckoff"] = wlabel
                d["site_symmetry"] = sitesym_labels[iatom]
                for s in frac_symbols:
                    d[s] = s
                if len(solutions) > 1:
                    cprint("Found multiple solutions for iatom %d" % iatom, "red")
                d.update(**solutions[0])
                rows.append(d)

        import pandas as pd
        df = pd.DataFrame(rows, index=None, columns=list(rows[0].keys()) if rows else None)
        return df

    def get_tensor_rank2_dataframe(self, view="all", select_symbols=None, verbose=0):
        """

        Args:
            view:
            select_symbols:
            verbose: Verbosity level

        Return |pandas-DataFrame| with cartesian tensor components as columns and (inequivalent) sites along the rows.
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols, verbose=verbose)

        sitesym_labels = self.structure.spget_site_symmetries()

        # Symmetric tensor.
        axx, ayy, azz, axy, axz, ayz = symbols = sp.symbols('axx, ayy, azz, axy, axz, ayz')
        tensor = sp.Matrix(([axx, axy, axz], [axy, ayy, ayz], [axz, ayz, azz]))

        indsym = self.structure.indsym
        rows = []
        for (iatom, wlabel) in zip(aview.iatom_list, aview.wyck_labels):
            site = self.structure[iatom]
            system = []
            for isym, rotf in enumerate(self.sp_symrel):
                if indsym[iatom, isym, 3] != iatom: continue
                m = rotf * tensor * rotf.T - tensor
                if verbose: sp.pprint(m)
                system.append(m)

            # Solve system of linear equations.
            solutions = sp.solve(system, dict=True)
            #print(solutions)
            if verbose and not solutions:
                cprint("No solution for iatom %d" % iatom, "yellow")

            if solutions:
                d = OrderedDict()
                d["element"] = site.specie.symbol
                d["site_index"] = iatom
                d["frac_coords"] = np.round(site.frac_coords, decimals=5)
                d["cart_coords"] = np.round(site.coords, decimals=5)
                d["wyckoff"] = wlabel
                d["site_symmetry"] = sitesym_labels[iatom]
                for s in symbols:
                    d[s] = s
                if len(solutions) > 1:
                    cprint("Found multiple solutions for iatom %d" % iatom, "red")
                d.update(**solutions[0])
                rows.append(d)

        import pandas as pd
        df = pd.DataFrame(rows, index=None, columns=list(rows[0].keys()) if rows else None)
        return df

    def check_site_symmetries(self, values_cart, verbose=0):
        """
        Test whether a set of tensors associated to the crystalline sites are compatible
        with the space group symmetries.

        Args:
            values_cart:
            verbose: Verbosity level.

        Return: max_err
        """
        natom = len(self.structure)
        indsym = self.structure.indsym
        nsym = len(self.symcart)

        max_err = 0.0
        for iatom in range(natom):
            ref_mat = values_cart[iatom]
            sym_mat = np.zeros_like(ref_mat)
            count = 0
            for isym, scart in enumerate(self.symcart):
                if indsym[iatom, isym, 3] != iatom: continue
                count += 1
                sym_mat += np.matmul(scart, np.matmul(ref_mat, scart.T))
                #sym_mat += np.matmul(scart.T, np.matmul(ref_mat, scart))

            if (nsym // count) * count != nsym: max_err = 1e+23
            sym_mat /= count
            diff_mat = sym_mat - ref_mat
            max_err = max(max_err, np.abs(diff_mat).sum())
            if count != 1 and verbose:
                print("For iatom", iatom, "count:", count, "ref_mat, sym_mat, diff_mat")
                print(np.hstack((ref_mat, sym_mat, diff_mat)))

        for iatom in range(natom):
            ref_mat = values_cart[iatom]
            for isym, scart in enumerate(self.symcart):
                jatom = indsym[iatom, isym, 3]
                if jatom == iatom: continue
                #sym_mat = np.matmul(scart, np.matmul(ref_mat, scart.T))
                sym_mat = np.matmul(scart.T, np.matmul(ref_mat, scart))
                diff_mat = sym_mat - values_cart[jatom]
                max_err = max(max_err, np.abs(diff_mat).sum())
                if verbose:
                    print("For iatom", iatom, "ref_mat, sym_mat, diff_mat")
                    print(np.hstack((values_cart[jatom], sym_mat, diff_mat)))

        print("Max error:", max_err)

        return max_err
