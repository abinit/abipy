# coding: utf-8
"""Classes for the analysis of electronic fatbands and projected DOSes."""
from __future__ import annotations

import traceback
import numpy as np

#from tabulate import tabulate
from numpy.linalg import inv, det, eigvals
from monty.termcolor import cprint
from monty.functools import lazy_property
from monty.string import marquee
#from pymatgen.core.periodic_table import Element
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.structure import Structure
from abipy.electrons.ebands import ElectronBands, ElectronsReader
#from abipy.tools.numtools import gaussian
from abipy.tools.typing import Figure
from abipy.tools.plotting import set_axlims, get_ax_fig_plt, get_axarray_fig_plt, add_fig_kwargs
#from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
#    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels, get_figs_plotly,
#    get_fig_plotly, add_plotly_fig_kwargs, PlotlyRowColDesc, plotly_klabels, plotly_set_lims)

def print_options_decorator(**kwargs):
    """
    A decorator to temporarily set np.printoptions inside a function.
    Use the decorator to apply print settings inside the function.

    @print_options_decorator(precision=2, suppress=True)
    def print_array():
        print(np.array([1.234567, 7.891234, 3.456789]))
    """
    def decorator(func):
        import functools
        @functools.wraps(func)
        def wrapper(*args, **kwargs_inner):
            with np.printoptions(**kwargs):
                return func(*args, **kwargs_inner)
        return wrapper
    return decorator



class OrbmagAnalyzer:
    """
    """

    def __init__(self, filepaths):
        self.orb_files = [OrbmagFile(path) for path in filepaths]

        # Consistency check
        structure0 = self.orb_files[0].structure
        if any(orb_file.structure != structure0 for orb_file in self.orb_files[1:]):
            raise RuntimeError("ORBMAG.nc files have different structures")

        # This piece of code is taken from merge_orbmag_mesh
        # the main difference is that here ncroots[0] is replaced by
        # the reader instance of the first OrbmagFile.
        r0 = self.orb_files[0].r

        mband = r0.read_dimvalue('mband')
        nkpt = r0.read_dimvalue('nkpt')
        nsppol = r0.read_dimvalue('nsppol')
        ndir = r0.read_dimvalue('ndir')
        orbmag_nterms = r0.read_dimvalue('orbmag_nterms')
        rprimd = r0.read_value('primitive_vectors')
        gprimd = inv(rprimd)
        ucvol = det(rprimd)

        # merging here means combine the 3 files: each delivers a 3-vector (ndir = 3),
        # output is a 3x3 matrix (ndir x ndir) for each term, sppol, kpt, band
        self.orbmag_merge_mesh = np.zeros((orbmag_nterms, ndir, ndir, nsppol, nkpt, mband))

        # here ncroots have been replaced by a list of reader instances.
        readers = [orb.r for orb in self.orb_files]

        for iband in range(mband):
            for isppol in range(nsppol):
                for ikpt in range(nkpt):
                    for iterm in range(orbmag_nterms):
                        for idir, r in enumerate(readers):
                            # convert terms to Cart coords; formulae differ depending on term. First
                            # four were in k space, remaining two already in real space
                            if iterm < 4:
                                omtmp = ucvol*np.matmul(gprimd, r.read_variable('orbmag_mesh')[iterm,0:ndir,isppol,ikpt,iband])
                            else:
                                omtmp = np.matmul(rprimd, r.read_variable('orbmag_mesh')[iterm,0:ndir,isppol,ikpt,iband])
                            self.orbmag_merge_mesh[iterm,idir,0:ndir,isppol,ikpt,iband] = omtmp

        # This piace of code has been taken from orbmag_sigij_mesh
        wtk = r0.read_value('kpoint_weights')
        occ = r0.read_value('occupations')
        orbmag_nterms = r0.read_dimvalue('orbmag_nterms')

        self.orbmag_merge_sigij_mesh = np.zeros((orbmag_nterms, nsppol, nkpt, mband, ndir, ndir))

        for iband in range(mband):
            for isppol in range(nsppol):
                for ikpt in range(nkpt):
                    # weight factor for each band and k point
                    trnrm=occ[0,ikpt,iband]*wtk[ikpt]/ucvol
                    for iterm in range(orbmag_nterms):
                        # sigij = \sigma_ij the 3x3 shielding tensor term for each sppol, kpt, and band
                        # additional ucvol factor converts to induced dipole moment (was dipole moment density,
                        # that is, magnetization)
                        self.orbmag_merge_sigij_mesh[iterm,isppol,ikpt,iband,0:ndir,0:ndir]= \
                        ucvol*trnrm*self.orbmag_merge_mesh[iterm,0:ndir,0:ndir,isppol,ikpt,iband]

    #def __str__(self) -> str:
    #    """String representation"""
    #    return self.to_string()

    #def to_string(self, verbose: int = 0) -> str:
    #    lines = []; app = lines.append
    #    return "\n".join(lines)

    #@print_options_decorator(precision=2, suppress=True)
    def report_eigvals(self, report_type):
        """
        """
        np.set_printoptions(precision=2)

        terms = ['CC   ','VV1  ','VV2  ','NL   ','LR   ','A0An ']
        orbmag_merge_sigij_mesh = self.orbmag_merge_sigij_mesh

        total_sigij = orbmag_merge_sigij_mesh.sum(axis=(0,1,2,3))
        eigenvalues = -1.0E6 * np.real(eigvals(total_sigij))
        isotropic = eigenvalues.sum() / 3.0
        span = eigenvalues.max() - eigenvalues.min()
        skew = 3.0 * (eigenvalues.sum() - eigenvalues.max() -eigenvalues.min() -isotropic) / span

        print('\nShielding tensor eigenvalues, ppm : ', eigenvalues)
        print('Shielding tensor iso, span, skew, ppm : %6.2f %6.2f %6.2f \n' % (isotropic,span,skew))

        if report_type == 'T':
            print('Term totals')
            term_sigij = orbmag_merge_sigij_mesh.sum(axis=(1, 2, 3))
            for iterm in range(np.size(orbmag_merge_sigij_mesh, axis=0)):
                eigenvalues = -1.0E6 * np.real(eigvals(term_sigij[iterm]))
                print(terms[iterm] + ': ', eigenvalues)
            print('\n')

        elif report_type == 'B':
            print('Band totals')
            band_sigij = orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2))
            for iband in range(np.size(orbmag_merge_sigij_mesh, axis=3)):
                eigenvalues = -1.0E6 * np.real(eigvals(band_sigij[iband]))
                print('band ' + str(iband) + ' : ', eigenvalues)
            print('\n')

        elif report_type == 'TB':
            print('Terms in each band')
            tband_sigij = orbmag_merge_sigij_mesh.sum(axis=(1,2))
            for iband in range(np.size(orbmag_merge_sigij_mesh, axis=3)):
                print('band ' + str(iband) + ' : ')
                for iterm in range(np.size(orbmag_merge_sigij_mesh, axis=0)):
                    eigenvalues = -1.0E6 * np.real(eigvals(tband_sigij[iterm, iband]))
                    print('   ' + terms[iterm] + ': ', eigenvalues)
                print('\n')
            print('\n')


    @add_fig_kwargs
    def plot_fatbands(self, ebands_kpath, ax=None, **kwargs) -> Figure:
        """
        Plot fatbands.

        Args:
            ebands_kpath
        """
        ebands_kpath = ElectronBands.as_ebands(ebands_kpath)

        # TODO: Use SKW to interpolate ...
        # I'd start by weighting each band and kpt by trace(sigij)/3.0, the isotropic part of sigij,
        # both as a term-by-term plot and as a single weighting produced by summing over all 6 terms.
        ebands_kpath.plot(ax=ax, show=False)

        #from abipy.core.skw import SkwInterpolator
        #skw = SkwInterpolator(lpratio, kpts, eigens, fermie, nelect, cell, symrel, has_timrev,
        #                      filter_params=None, verbose=1)
        #print(skw)
        #skw.eval_sk(spin, kpt, der1=None, der2=None) -> np.ndarray:

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        return fig


class OrbmagFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    """

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram::
    """

    @classmethod
    def from_file(cls, filepath: str) -> OrbmagFile:
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = ElectronsReader(filepath)

        # Initialize the electron bands from file
        self.natom = len(self.structure)

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self) -> dict:
        """dict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def close(self) -> None:
        """Called at the end of the ``with`` context manager."""
        return self.r.close()

    def __str__(self) -> str:
        """String representation"""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=True, title="Electronic Bands"))
        app("")
        #app(marquee("Fatbands Info", mark="="))
        #app("prtdos: %d, prtdosm: %d, mbesslang: %d, pawprtdos: %d, usepaw: %d" % (
        #    self.prtdos, self.prtdosm, self.mbesslang, self.pawprtdos, self.usepaw))
        app("nsppol: %d, nkpt: %d, mband: %d" % (self.nsppol, self.nkpt, self.mband))
        app("")

        if verbose > 1:
            app("")
            app(self.hdr.to_str(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)
