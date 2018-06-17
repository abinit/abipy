# coding: utf-8
"""Interface to the ABIWAN netcdf file produced by abinit by calling wannier90 in library mode."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

#from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt #, get_axarray_fig_plt
from abipy.electrons.ebands import ElectronBands, ElectronsReader


class AbiwanFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File produced by Abinit with the unitary matrices obtained by
    calling wannier90 in library mode.

    Usage example:

    .. code-block:: python

        with abilab.abiopen("foo_ABIWAN.nc") as abiwan:
            print(abiwan)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbiwanFile
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(AbiwanFile, self).__init__(filepath)

        self.reader = AbiwanReader(filepath)

        # Here be very careful with F --> C
        #complex(dp) U_matrix(mwan,mwan,nkpt,nsppol)
        #complex(dp) U_matrix_opt(mband,mwan,nkpt,nsppol)
        # complex(dpc) :: hamWR(:,:,:,:)
        # ! hamWR(mwan,mwan,nrpts,nsppol))
        # ! Hamiltonian in k-space (ab-initio grid) in the Wannier gauge.

    @lazy_property
    def nwan_spin(self):
        """Number of Wannier functions for each spin."""
        return self.reader.read_value("nwan")

    @lazy_property
    def mwan(self):
        """Max number of Wannier functions over spins, i.e max(nwan_spin) (used to dimension arrays)."""
        return self.reader.read_dimvalue("mwan")

    #@lazy_property
    #def nntot(self):
    #    """Number of k-point neighbours."""
    #    return self.reader.read_dimvalue("nntot")

    #@lazy_property
    #def ndimwin(self):
    #    """
    #    [nsppol, nkpt] array giving the number of bands inside the outer window for each k-point and spin.
    #    """
    #    return self.reader.read_value("ndimwin")

    #@lazy_property
    #def band_in(self):
    #    """[nsppol, mband] logical array. True if (spin, band) is included in the calculation."""
    #    return self.reader.read_value("band_in")

    #@lazy_property
    #def lwindow(self):
    #    """
    #    [nsppol, mband, nkpt] array. Only if disentanglement
    #    True if this band at this k-point lies within the outer window
    #    """
    #    return self.reader.read_value("lwindow")

    #@lazy_property
    #def have_disentangled(self):
    #    """[nsppol] array. Whether a disentanglement has been performed."""
    #    return self.reader.read_value("have_disentangled")

    @lazy_property
    def wann_centers(self):
        """Numpy array of shape (nsppol, mwan, 3). with Wannier centers in Ang."""
        return self.reader.read_value("wann_centres")

    @lazy_property
    def wann_spreads(self):
        """Numpy array of shape (nsppol, mwan) with spreads."""
        return self.reader.read_value("wann_spreads")

    @lazy_property
    def irvec(self):
        """
        [nrpts, 3] array with the lattice vectors in the WS cell
        in the basis of the lattice vectors defining the unit cell
        """
        return self.reader.read_value("irvec")

    @lazy_property
    def ndegen(self):
        """
        [nrpts] array with the degeneracy of each point.
        It will be weighted using 1 / ndegen[ir]
        """
        return self.reader.read_value("ndegen")

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        # TODO
        return od

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level verbose."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        app(marquee("Wannier centers and spreads", mark="="))
        for spin in range(self.nsppol):
            for iwan in range(self.nwan_spin[spin]):
                app("Center %.6f %.6f %.6f, Spread: %.6f" % (
                    self.wann_centers[spin, iwan, 0],
                    self.wann_centers[spin, iwan, 1],
                    self.wann_centers[spin, iwan, 2],
                    self.wann_spreads[spin, iwan]))

        #self.have_disentangled:
        #self.band_in
        #self.lwindow

        #if verbose > 1:
        #    app("")
        #    app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    def close(self):
        self.reader.close()

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    def get_hwr(self):
        """
        Construct the matrix elements of the KS Hamiltonian in real space
        """
        self.hwr_spin = [None] * self.nsppol
        nr = len(self.irvec)
        #nkpt = ??
        for spin in range(self.nsppol):
            nwan = self.nwan_spin[spin]
            hr = self.hwr_spin[spin] = np.empty((nr, nwan, nwan), dtype=np.complex)
            ham_k = np.zeros((nkpt, nwan, nwan), dtype=np.complex)
            eigval2 = np.zeros((nkpt, nwan), dtype=np.double)

            #! Real-space Hamiltonian H(R) is calculated by Fourier
            #! transforming H(q) defined on the ab-initio reciprocal mesh
            #!
            #allocate(HH_q(num_wann,num_wann,num_kpts))
            #allocate(num_states(num_kpts))

            #HH_q=cmplx_0
            #do ik=1,num_kpts
            #   if(have_disentangled) then
            #      num_states(ik)=ndimwin(ik)
            #   else
            #      num_states(ik)=num_wann
            #   endif
            #   call get_win_min(ik,winmin_q)
            #   do m=1,num_wann
            #      do n=1,m
            #         do i=1,num_states(ik)
            #            ii=winmin_q+i-1
            #            HH_q(n,m,ik)=HH_q(n,m,ik)&
            #             +conjg(v_matrix(i,n,ik))*eigval(ii,ik)&
            #             *v_matrix(i,m,ik)
            #         enddo
            #         HH_q(m,n,ik)=conjg(HH_q(n,m,ik))
            #      enddo
            #   enddo
            #enddo
            #call fourier_q_to_R(HH_q,HH_R)

    def interpolate_ebands(self, vertices_names=None, line_density=20):
        """

        Args:
            vertices_names: Used to specify the k-path for the interpolated QP band structure
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.

        Returns: |ElectronBands| object with interpolated energies.
        """
        if vertices_names is None:
            vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
        kpoints = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)
        kfrac_coords, knames = kpoints.frac_coords, kpoints.names

        # NB: May have different number of wannier functions if nsppol == 2
        nk = len(kpoints)
        eigens = np.zeros((self.nsppol, nk, self.mwan)
        occfacts = np.zeros(eigens.shape)

        # Interpolate the Hamiltonian for each kpoint and spin
        for spin in range(self.nsppol):
            for ik, kpt in enumerate(kpoints):
               oeigs = self.hwr.eval_sk(spin, kpt.frac_coords, der1=None, der2=None)
               eigens[spin, ik, :] = oeigs

        return ElectronBands(self.structure, kpoints, eigens, self.ebands.fermie,
                             occfacts, self.ebands.nelect, self.nspinor, self.nspden)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.ebands.plot(show=False)
        if kwargs.get("verbose"):
            ebw = self.interpolate_ebands()
            yield ebw.plot(show=False)
        else:
            cprint("Use verbose option to plot interpolated band and H(R)", "yellow")

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abiwan = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abiwan)"),
            nbv.new_code_cell("abiwan.ebands.plot();"),
            nbv.new_code_cell("abiwan.ebands.kpoints.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


# TODO: Interface with ElectronsInterpolator
class HWannierR(object)

    def __init__(self, irvec, ndegen):
        self.irvec = irvec
        self.ndegen = ndegen

    def eval_sk(self, spin, kpt, der1=None, der2=None):
        """
        Interpolate eigenvalues for all bands at a given (spin, k-point).
        Optionally compute gradients and Hessian matrices.

        Args:
            spin: Spin index.
            kpt: K-point in reduced coordinates.
            der1: If not None, ouput gradient is stored in der1[nband, 3].
            der2: If not None, output Hessian is der2[nband, 3, 3].

        Return:
            oeigs[nband]
        """
        #call pw90common_fourier_R_to_k(kpt,HH_R,HH,0)
        #call utility_diagonalize(HH,num_wann,my_eig(:,loop_kpt),UU)

        #  !! For alpha=0:
        #  !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
        #  !!
        #  !! For alpha=1,2,3:
        #  !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
        #  !! where R_alpha is a Cartesian component of R
        #  !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)

        #    real(kind=dp)                                   :: kpt(3)
        #    complex(kind=dp), dimension(:,:,:), intent(in)  :: OO_R
        #    complex(kind=dp), dimension(:,:), intent(out)   :: OO
        #    integer                                         :: alpha

        #    integer          :: ir, i,j,ideg
        #    real(kind=dp)    :: rdotk
        #    complex(kind=dp) :: phase_fac

        #    OO(:,:)=cmplx_0
        #    do ir=1,nrpts
        #        ! [lp] Original code, without IJ-dependent shift:
        #         rdotk=twopi*dot_product(kpt(:),irvec(:,ir))
        #         phase_fac=cmplx(cos(rdotk),sin(rdotk),dp)/real(ndegen(ir),dp)
        #         if(alpha==0) then
        #            OO(:,:)=OO(:,:)+phase_fac*OO_R(:,:,ir)
        #         elseif(alpha==1.or.alpha==2.or.alpha==3) then
        #            OO(:,:)=OO(:,:)+&
        #                  cmplx_i*crvec(alpha,ir)*phase_fac*OO_R(:,:,ir)
        #         else
        #            stop 'wrong value of alpha in pw90common_fourier_R_to_k'
        #         endif
        #    enddo

        #  end subroutine pw90common_fourier_R_to_k

    @add_fig_kwargs
    def plot_hwr(self, ax=None, **kwargs):
        """
        Plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Sort points by length
        # Length(ii), Lattice point Rws(1:3,ii), (MAXVAL |hamW_ij(R,isp)|, isp=1,nsppol)'
        for spin in range(self.nsppol):
            nwan = self.nwan_spin[spin]
            hr = self.hwr_spin[spin]

        return fig


class AbiwanReader(ElectronsReader):
    """
    This object reads the results stored in the ABIWAN file produced by ABINIT after having
    called wannier90 in library mode.
    It provides helper function to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbiwanReader
    """
