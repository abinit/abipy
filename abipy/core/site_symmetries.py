"""
This module provides functions to generate crystalline structures from Wyckoff positions
or to retrieve Wyckoff parameters from a given structure.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import sympy as sp

#from monty.functools import lazy_property
from collections import OrderedDict
from monty.termcolor import cprint
from abipy.core.mixins import Has_Structure


class SiteSymmetries(Has_Structure):

    def __init__(self, structure):
        # Get Spacegroup from spglib if not available from Abinit.
        if not structure.has_abi_spacegroup:
            structure.spgset_abi_spacegroup(has_timerev=True, overwrite=False)

        self._structure = structure

        abispg = structure.abi_spacegroup
        nsym = len(abispg.symrel)
        indsym = self.structure.indsym

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

        import spglib
        self.sitesym_labels = []
        for iatom, site in enumerate(self.structure):
            rotations = [abispg.symrel[isym] for isym in range(nsym) if
                         indsym[iatom, isym, 3] == iatom and abispg.symafm[isym] == 1]
            # Passing a 0-length rotations list to spglib can segfault.
            herm_symbol, ptg_num = "1", 1
            if len(rotations) != 0:
                herm_symbol, ptg_num, trans_mat = spglib.get_pointgroup(rotations)

            self.sitesym_labels.append("%s (#%d,nsym:%d)" % (herm_symbol.strip(), ptg_num, len(rotations)))

        #for irred_isite in self.irred_isites:
        #    for isite_eq, rm1, tau, l0 in self.eq_sites[irred_isite]:
        #        # isite_eq = rm1(irred_site - tau) + l0

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level verbose."""
        lines = []; app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        return "\n".join(lines)

    def get_wyckoff_dataframe(self, view="all", select_symbols=None, decimals=5, verbose=0):
        """
        Find Wyckoff positions.

        Args:
            view:
            select_symbols:
            decimals: Number of decimal places to round to.
                If decimals is negative, it specifies the number of positions to the left of the decimal point.
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
                if verbose:
                    print(92 * "=")
                    print("System of linear equations for iatom %d, %s" % (iatom, repr(site)))
                    sp.pprint(m)
                    print("Rotation:")
                    sp.pprint(rm1)

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
                d["cart_coords"] = np.round(site.coords, decimals=decimals)
                d["frac_coords"] = np.round(site.frac_coords, decimals=decimals)
                d["wyckoff"] = wlabel
                d["site_symmetry"] = sitesym_labels[iatom]
                for s in frac_symbols:
                    d[str(s)] = str(s)
                if len(solutions) > 1:
                    cprint("Found multiple solutions for iatom %d" % iatom, "red")
                d.update({str(k): str(v) for k, v in solutions[0].items()})
                rows.append(d)

        import pandas as pd
        df = pd.DataFrame(rows, index=None, columns=list(rows[0].keys()) if rows else None)
        return df

    def get_tensor_rank2_dataframe(self, view="all", select_symbols=None, decimals=5, verbose=0):
        """

        Args:
            view:
            select_symbols:
            decimals: Number of decimal places to round to.
                If decimals is negative, it specifies the number of positions to the left of the decimal point.
            verbose: Verbosity level

        Return |pandas-DataFrame| with cartesian tensor components as columns and (inequivalent) sites along the rows.
        """
        # Select atoms.
        aview = self._get_atomview(view, select_symbols, verbose=verbose)

        sitesym_labels = self.structure.spget_site_symmetries()

        # Symmetric tensor in reduced coordinates (direct lattice)
        # Operations in reduced coords are given by integer matrices and this facilitates the solution
        # of the system of equations with sympy.
        Txx, Tyy, Tzz, Txy, Txz, Tyz = symbols = sp.symbols('Txx, Tyy, Tzz, Txy, Txz, Tyz')
        tensor = sp.Matrix(([Txx, Txy, Txz], [Txy, Tyy, Tyz], [Txz, Tyz, Tzz]))

        indsym = self.structure.indsym
        rows = []
        for (iatom, wlabel) in zip(aview.iatom_list, aview.wyck_labels):
            site = self.structure[iatom]
            system = []
            for isym, rotf in enumerate(self.sp_symrel):
                if indsym[iatom, isym, 3] != iatom: continue
                m = rotf * tensor * rotf.T - tensor
                if verbose:
                    print(92 * "=")
                    print("System of linear equations for iatom %d, %s" % (iatom, repr(site)))
                    sp.pprint(m)
                    print("Rotation:")
                    sp.pprint(rotf)
                    print(92 * "=")

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
                d["frac_coords"] = np.round(site.frac_coords, decimals=decimals)
                d["cart_coords"] = np.round(site.coords, decimals=decimals)
                d["wyckoff"] = wlabel
                d["site_symmetry"] = sitesym_labels[iatom]
                for s in symbols:
                    d[str(s)] = str(s)
                if len(solutions) > 1:
                    cprint("Found multiple solutions for iatom %d" % iatom, "red")
                d.update({str(k): str(v) for k, v in solutions[0].items()})
                rows.append(d)

        import pandas as pd
        df = pd.DataFrame(rows, index=None, columns=list(rows[0].keys()) if rows else None)
        return df

    def check_site_symmetries(self, tcart, verbose=0):
        """
        Test whether a set of tensors associated to the crystalline sites are compatible
        with the space group symmetries.

        Args:
            tcart: (natom, 3, 3) array
            verbose: Verbosity level.

        Return: max_err
        """
        natom = len(self.structure)
        indsym = self.structure.indsym
        nsym = len(self.symcart)
        tcart = np.reshape(tcart, (natom, 3, 3))

        max_err = 0.0
        for iatom in range(natom):
            ref_mat = tcart[iatom]
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
            ref_mat = tcart[iatom]
            for isym, scart in enumerate(self.symcart):
                jatom = indsym[iatom, isym, 3]
                if jatom == iatom: continue
                #sym_mat = np.matmul(scart, np.matmul(ref_mat, scart.T))
                sym_mat = np.matmul(scart.T, np.matmul(ref_mat, scart))
                diff_mat = sym_mat - tcart[jatom]
                max_err = max(max_err, np.abs(diff_mat).sum())
                if verbose:
                    print("For iatom", iatom, "ref_mat, sym_mat, diff_mat")
                    print(np.hstack((tcart[jatom], sym_mat, diff_mat)))

        print("Max error:", max_err)

        return max_err
